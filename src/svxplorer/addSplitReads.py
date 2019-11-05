#!/usr/bin/env python

# Add Split Reads to support existing PE variants and create new SR variants

import pysam
import sys
import argparse
import logging
import gc
from shared import formExcludeHash, ignoreRead, readChromosomeLengths, countLines

#global variables
detectIntraChrCopyInv = 0

class newSRVar(object):
    def __init__(self):
        self.l_orient = -1
        self.r_orient = -1
        self.swapped = -1
        self.bp2 = -1
        self.bp3 = -1
        self.count = 1
        self.support = []
        self.bp3tid = -1
        self.typeSV = -1
        self.tag = -1
        self.neighbor_tags = []
        self.hash_pair_tag = -1
        self.n_changes = 0
        self.write = -1
        self.isOriginal = -1
        self.orientStatus = -1
    def __str__(self):
       return "%s\t%s\t%s\t%s\t%s" %(self.bp2, self.bp3, self.count, self.isOriginal, self.typeSV)

class PEVarDetails(object):
    def __init__(self):
        self.bp1_2 = -1
        self.bp2_1 = -1
        self.bp2_2 = -1
        self.bp3_1 = -1
        self.bp3_2 = -1
        self.typeSV = -1
        self.num = -1
        self.orient = "22"

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" %(self.bp2_1, self.bp2_2, self.bp3_1, self.bp3_2, self.typeSV, self.num)

def transferSupport(variant1, variant2):
    # 1st fragment is shared
    for elem in variant2.support[1:]:
        variant1.support.append(elem)
        variant1.count+=1

def mapSVtoNum(SV_type):
    if SV_type== "DEL":
        return 0
    elif SV_type == "TD":
        return 1
    elif SV_type == "INV":
        return 2
    elif SV_type in ["INS","INS_C","INS_C_P"]:
        if SV_type == "INS_C":
            return 5
        return 3
    elif SV_type in ["INS_I", "INS_C_I", "INS_C_I_P"]:
        if SV_type == "INS_C_I":
            return 6
        return 4
    else:
        return -1

def formPEHash(fAV, iObjects, slop):
    logging.info('Started reading the PE variants')
    global SVHashPE
    for line_num, line in enumerate(fAV):
        line_s = line.split()
        SV_specsPE = PEVarDetails()
        SV_specsPE.num = int(line_s[0])
        SV_specsPE.typeSV= mapSVtoNum(line_s[1])
        if SV_specsPE.typeSV == -1:
            continue
        if SV_specsPE.typeSV in [3,4,5,6]:
            SV_specsPE.bp3_1 = int(line_s[9]) - slop
            SV_specsPE.bp3_2 = int(line_s[10]) + slop
        SV_specsPE.bp1_2 = int(line_s[4])
        SV_specsPE.bp2_1 = int(line_s[6]) - slop
        SV_specsPE.bp2_2 = int(line_s[7]) + slop
        SV_specsPE.orient = line_s[15]
        # hash all values within bp margin
        bp = int(line_s[3])
        tid1 = line_s[2]
        tid2 = line_s[5]
        almt = (tid1, tid2, bp)
        #immutable hash objects -- preserve so memory doesn't get overwritten
        iObjects[line_num] = almt
        if almt not in SVHashPE:
            SVHashPE[almt] = SV_specsPE
    logging.info('Finished reading the PE variants')
    return SV_specsPE.num

def addSplitReads(workDir, variantMapFilePE, allVariantFilePE, bamFileSR,
                  slop, refRate, min_vs, mapThresh, ignoreChr, minSizeINS,
                  minSRtoPEsupport, ignoreBED, noCCleanUp, maxClusterMargin):
    fAV = open(workDir+"/allVariants.pe.txt","r")
    fVM = open(workDir+"/variantMap.pe.txt","r")
    fAVN = open(workDir+"/allVariants.pe_sr.txt","w")
    fVMN = open(workDir+"/variantMap.pe_sr.txt","w")
    global SVHashPE
    SVHashPE = {}
    SRVarHash = {}
    # preserve list of complex hash objects
    MAX_ARRAY_SIZE = countLines(workDir+"/allVariants.pe.txt")
    immutable_objects = [None]*MAX_ARRAY_SIZE

    # save the PE variants
    headerAV = fAV.readline()
    nSVsPE = formPEHash(fAV, immutable_objects, slop)

    SRFrag = 0
    SRtoPESuppFrags = [[] for _ in range(1+nSVsPE)]
    SRtoPESuppBPs = {}
    newSRList = []
    bamfile = pysam.Samfile(bamFileSR,"rb")
    # if subsampling: shouldn't be required`
    bp1Prev = -1
    bp1TID = -1
    bp2Prev = -1
    bp2TID = -1
    ins_slop = 30

    ignoreList = set()
    ignoreTIDAll = set()
    if ignoreChr is not None:
        with open(ignoreChr, "r") as fIC:
            for line in fIC:
                chrI = line.strip().split()[0]
                if not chrI.startswith("*"):
                    ignoreList.add(chrI)
                    logging.info("Adding SRs: Chromosome %s will be ignored.", chrI)
                else:
                    ignoreTIDAll.add(chrI[1:])
                    logging.info("Adding SRs: Chr names starting with %s will be ignored", chrI[1:])

    chromosome_lengths = readChromosomeLengths(bamFileSR)

    chrHash = {}
    if ignoreBED is not None:
        logging.info("Regions in %s will be ignored", ignoreBED)
        chrHash = formExcludeHash(chrHash, 0, ignoreBED, chromosome_lengths)
    # all split reads should be mapped, unique alignments
    counterSR = 0
    while True:
        try:
            sr1 = bamfile.next()
            sr2 = bamfile.next()
        except StopIteration:
            break
        varType = -1
        if sr1.qname == sr2.qname:
            SRFrag-=1
            sr_bp1 = sr1.reference_start
            sr_bp2 = sr2.reference_start
            sr_bp1_tid = sr1.reference_name
            sr_bp2_tid = sr2.reference_name
            counterSR+=1
            if counterSR % 10000 == 0:
                logging.debug("Processed %s SRs", counterSR)

            # ignore marked chromosomes
            if sr1.reference_name in ignoreList or \
               sr2.reference_name in ignoreList or \
               sr1.mapping_quality < mapThresh or \
               sr2.mapping_quality < mapThresh or \
               ignoreRead(sr_bp1_tid, sr_bp1, sr_bp2_tid, sr_bp2, chrHash) or \
               (sr_bp1_tid == bp1TID and sr_bp2_tid == bp2TID and abs(sr_bp1 - bp1Prev) < refRate and abs(sr_bp2 - bp2Prev) < refRate):
                    continue

            skip = 0
            for chrI in ignoreTIDAll:
                if sr1.reference_name.startswith(chrI) or sr2.reference_name.startswith(chrI):
                    skip = 1
                    break
            if skip:
                logging.debug("Ignoring SR almt %s and %s, as occurs as * entry in ignoreTIDs", sr1, sr2)
                continue
            bp1Prev = sr_bp1
            bp1TID = sr_bp1_tid
            bp2Prev = sr_bp2
            bp2TID = sr_bp2_tid

            ## SET SWAP AND RISK
            if sr_bp1 < sr_bp2:
                minsr = sr1
                maxsr = sr2
            else:
                maxsr = sr1
                minsr = sr2
                sr_bp2, sr_bp1 = sr_bp1, sr_bp2
            sr_bp1_tid = minsr.reference_name
            sr_bp2_tid = maxsr.reference_name
            # QAS below refers to the alignment position of split read relative to whole read
            # swap value accurate for 75% of inversion reads but even if incorrect, unused
            if sr_bp1_tid == sr_bp2_tid and minsr.query_alignment_start > maxsr.query_alignment_start:
                swap = 1
                if minsr.is_reverse == maxsr.is_reverse:
                    # unless both reads are inverted and swapped -- less likely
                    sr_bp1 = minsr.reference_start
                    sr_bp2 = maxsr.reference_end
                # "risk" not catching some inverted copy-paste insertions but these are few
                else:
                    sr_bp1 = minsr.reference_end
                    sr_bp2 = maxsr.reference_end
            # note: query start positions (QAS) of both reads can be 0, which happens in 50% of SRs from INVs
            else:
                swap = 0
                if minsr.is_reverse == maxsr.is_reverse:
                    sr_bp1 = minsr.reference_end
                    sr_bp2 = maxsr.reference_start
                else:
                    if minsr.is_reverse:
                        sr_bp1 = minsr.reference_start
                        sr_bp2 = maxsr.reference_start
                    else:
                        sr_bp1 = minsr.reference_end
                        sr_bp2 = maxsr.reference_end
        else:
            sys.stderr.write("Please check if split-read file is name-sorted. Exactly two entries per read name.")
            exit(1)

        ## CHECK CURRENT SR ALMT AGAINST EXISTING PE VARIANTS FOR MATCH
        match = 0
        peFound = False
        for x in range(sr_bp1 + slop, sr_bp1 - maxClusterMargin - slop,-1):
            searchAlmt = (sr_bp1_tid, sr_bp2_tid, x)
            if searchAlmt in SVHashPE and sr_bp1 < SVHashPE[searchAlmt].bp1_2 and \
                ((SVHashPE[searchAlmt].bp2_1 - slop < sr_bp2 < SVHashPE[searchAlmt].bp2_2 + slop) or \
                (SVHashPE[searchAlmt].bp3_1 - slop < sr_bp2 < SVHashPE[searchAlmt].bp3_2 + slop)):
                peFound = True
                break
        if peFound and (SVHashPE[searchAlmt].bp2_1 - slop < sr_bp2 < SVHashPE[searchAlmt].bp2_2 + slop): 
            varNumPE = SVHashPE[searchAlmt].num
            varType = SVHashPE[searchAlmt].typeSV
            if varNumPE not in SRtoPESuppBPs:
                newBp = [sr_bp1, sr_bp2, -1, -1]
        elif peFound and (SVHashPE[searchAlmt].bp3_1 - slop < sr_bp2 < SVHashPE[searchAlmt].bp3_2 + slop):

            varNumPE = SVHashPE[searchAlmt].num
            varType = SVHashPE[searchAlmt].typeSV
            if varNumPE not in SRtoPESuppBPs:
                newBp = [sr_bp1, -1, sr_bp2, -1]
        # if not found, try other side of SR
        else:
            #should be unset anyway
            peFound = False
            for x in range(sr_bp2 + slop, sr_bp2 - maxClusterMargin - slop,-1):
                searchAlmt = (sr_bp2_tid, sr_bp1_tid, x)
                if searchAlmt in SVHashPE and sr_bp2 < SVHashPE[searchAlmt].bp1_2 and \
                    ((SVHashPE[searchAlmt].bp2_1 - slop < sr_bp1 < SVHashPE[searchAlmt].bp2_2 + slop) or \
                    (SVHashPE[searchAlmt].bp3_1 - slop < sr_bp1 < SVHashPE[searchAlmt].bp3_2 + slop)):
                    peFound = True
                    break

            if peFound and (SVHashPE[searchAlmt].bp2_1 - slop < sr_bp1 < SVHashPE[searchAlmt].bp2_2 + slop):
                varNumPE = SVHashPE[searchAlmt].num
                varType = SVHashPE[searchAlmt].typeSV
                if varNumPE not in SRtoPESuppBPs:
                    newBp = [sr_bp2, sr_bp1, -1, -1]
            elif peFound and (SVHashPE[searchAlmt].bp3_1 - slop < sr_bp1 < SVHashPE[searchAlmt].bp3_2 + slop):
                varNumPE = SVHashPE[searchAlmt].num
                varType = SVHashPE[searchAlmt].typeSV
                if varNumPE not in SRtoPESuppBPs:
                    newBp = [sr_bp2, -1, sr_bp1, -1]

        #check DEL
        if varType == 0 and swap==0 and minsr.is_reverse == maxsr.is_reverse:
            match = 1
        #check TD
        elif varType == 1 and swap==1 and minsr.is_reverse == maxsr.is_reverse:
            match = 1
        #check INV
        elif varType == 2 and swap==0 and minsr.is_reverse != maxsr.is_reverse:
            match = 1
            # if both ends of inversion confirmed as contiguous with reference
            if (SVHashPE[searchAlmt].orient == "00" and minsr.reference_end > sr_bp1 and maxsr.reference_end > sr_bp2) or \
               (SVHashPE[searchAlmt].orient == "11" and minsr.reference_start < sr_bp1 and maxsr.reference_start < sr_bp2):
                newBp[3] = 1                               
        #check INS
        elif (varType in [3,5] and minsr.is_reverse == maxsr.is_reverse) or \
            (varType in [4,6] and minsr.is_reverse != maxsr.is_reverse):
            match = 1
            logging.debug("INS match of frag %d with PE cluster %d", SRFrag, varNumPE)
            # Set new bp3 of insertion if unset
            if varNumPE in SRtoPESuppBPs and (SRtoPESuppBPs[varNumPE][2] == -1 \
                or SRtoPESuppBPs[varNumPE][1] == -1):

                SRtoPESupp_bp = SRtoPESuppBPs[varNumPE]
                if SRtoPESuppBPs[varNumPE][2] == -1:
                    if abs(sr_bp1 - SRtoPESupp_bp[0]) > ins_slop and abs(sr_bp1 - SRtoPESupp_bp[1]) > ins_slop \
                        and SVHashPE[searchAlmt].bp3_1 < sr_bp1 < SVHashPE[searchAlmt].bp3_2:
                        SRtoPESuppBPs[varNumPE][2] = sr_bp1
                    elif abs(sr_bp2 - SRtoPESupp_bp[0]) > ins_slop and abs(sr_bp2 - SRtoPESupp_bp[1]) > ins_slop and \
                        SVHashPE[searchAlmt].bp3_1 < sr_bp2 < SVHashPE[searchAlmt].bp3_2:
                        SRtoPESuppBPs[varNumPE][2] = sr_bp2
                elif SRtoPESuppBPs[varNumPE][1] == -1:
                    if abs(sr_bp1 - SRtoPESupp_bp[0]) > ins_slop and abs(sr_bp1 - SRtoPESupp_bp[2]) > ins_slop \
                        and SVHashPE[searchAlmt].bp2_1 < sr_bp1 < SVHashPE[searchAlmt].bp2_2:
                        SRtoPESuppBPs[varNumPE][1] = sr_bp1
                    elif abs(sr_bp2 - SRtoPESupp_bp[0]) > ins_slop and abs(sr_bp2 - SRtoPESupp_bp[2]) > ins_slop and \
                        SVHashPE[searchAlmt].bp2_1 < sr_bp2 < SVHashPE[searchAlmt].bp2_2:
                        SRtoPESuppBPs[varNumPE][1] = sr_bp2

                # insertion bp2 should be < bp3
                if varType in [3,4] and SRtoPESuppBPs[varNumPE][1] != -1 \
                    and SRtoPESuppBPs[varNumPE][2] != -1 and SRtoPESuppBPs[varNumPE][2] < SRtoPESuppBPs[varNumPE][1]:
                    SRtoPESuppBPs[varNumPE][1], SRtoPESuppBPs[varNumPE][2] =\
                        SRtoPESuppBPs[varNumPE][2], SRtoPESuppBPs[varNumPE][1]
        # if matches existing PE SV
        if match:
            if varNumPE not in SRtoPESuppBPs:
                SRtoPESuppBPs[varNumPE] = newBp
                SRtoPESuppFrags[varNumPE].append(SRFrag)
            else:
                SRtoPESuppFrags[varNumPE].append(SRFrag)
        # if SR fragment did not match existing PE variant, check existing SR variant list
        else:
            almtMatchesSVbps = 0
            almtSupportsSV = 0
            listBPSearch = range(int(sr_bp1-slop/2),int(sr_bp1+slop/2) + 1)
            for x in listBPSearch:
                newAlmt = (sr_bp1_tid, sr_bp2_tid, x)

                if newAlmt in SRVarHash:
                    almtMatchesSVbps = 1
                    other_bp = sr_bp2
                    other_bp_tid = sr_bp2_tid
                    break
            if almtMatchesSVbps == 0:
                listBPSearch = range(int(sr_bp2-slop/2),int(sr_bp2+slop/2) + 1)
                for x in listBPSearch:
                    newAlmt = (sr_bp2_tid, sr_bp1_tid, x) 
                    if newAlmt in SRVarHash:
                        almtMatchesSVbps = 1
                        other_bp = sr_bp1
                        other_bp_tid = sr_bp1_tid
                        break

            # bp1 is left bp, bp2 is right bp if same chr
            l_orient = minsr.is_reverse
            r_orient = maxsr.is_reverse

            ## IF FRAGMENT SUPPORTS EXISTING SR VARIANT
            if almtMatchesSVbps == 1:
                # do not look for inverted insertions due to ambiguity of swap parameter
                if SRVarHash[newAlmt].typeSV == "DEL_INS" or SRVarHash[newAlmt].typeSV == "DEL":
                    if l_orient == r_orient and swap==1 and not (SRVarHash[newAlmt].bp2-slop/2 <= other_bp \
                        <= SRVarHash[newAlmt].bp2+slop/2) and not (SRVarHash[newAlmt].bp2 < newAlmt[2] < other_bp or \
                        SRVarHash[newAlmt].bp2 > newAlmt[2] > other_bp): 

                        SRVarHash[newAlmt].typeSV = "INS"
                        SRVarHash[newAlmt].bp3 = other_bp
                        SRVarHash[newAlmt].bp3tid = other_bp_tid
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                        SRVarHash[newAlmt].n_changes+=1
                        almtSupportsSV = 1
                        #print "DEL_INS"
                    # cannot update type based on this info
                    elif swap==0 and SRVarHash[newAlmt].bp2-slop/2 <= other_bp <= SRVarHash[newAlmt].bp2+slop/2 \
                        and l_orient == r_orient == SRVarHash[newAlmt].l_orient:
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                    elif swap==0 and SRVarHash[newAlmt].bp2-slop/2 <= other_bp <= SRVarHash[newAlmt].bp2+slop/2 \
                        and l_orient == r_orient and l_orient != SRVarHash[newAlmt].l_orient:
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                #TD_I
                elif SRVarHash[newAlmt].typeSV == "TD_I":
                #cannot update type here
                    if swap==1 and l_orient == r_orient and SRVarHash[newAlmt].bp2-slop/2 <= other_bp \
                        <= SRVarHash[newAlmt].bp2+slop/2:
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                        #print "TD 1", SRVarHash[newAlmt].count

                    elif l_orient == r_orient and swap==0 and not (SRVarHash[newAlmt].bp2-slop/2 <= other_bp \
                        <= SRVarHash[newAlmt].bp2+slop/2) and not (other_bp < newAlmt[2] < SRVarHash[newAlmt].bp2 or \
                        other_bp > newAlmt[2] > SRVarHash[newAlmt].bp2): 

                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                        SRVarHash[newAlmt].typeSV = "INS"
                        SRVarHash[newAlmt].bp3 = other_bp
                        SRVarHash[newAlmt].bp3tid = other_bp_tid
                        SRVarHash[newAlmt].n_changes+=1
                elif SRVarHash[newAlmt].typeSV[:3] == "INV" and l_orient != r_orient:
                    # if both ends of inversion contiguous with reference, confirmed
                    if (SRVarHash[newAlmt].orientStatus == 1 and \
                        minsr.reference_start < sr_bp1 and maxsr.reference_start < sr_bp2) or \
                        (SRVarHash[newAlmt].orientStatus == 0 and \
                        minsr.reference_end > sr_bp1 and maxsr.reference_end > sr_bp2):
                        SRVarHash[newAlmt].typeSV = "INV_B"
                        SRVarHash[newAlmt].orientStatus = 2
                    elif SRVarHash[newAlmt].orientStatus == -1:
                        if minsr.reference_start < sr_bp1 and maxsr.reference_start < sr_bp2:
                            SRVarHash[newAlmt].orientStatus = 0
                        elif minsr.reference_end > sr_bp1 and maxsr.reference_end > sr_bp2:
                            SRVarHash[newAlmt].orientStatus = 1

                    if SRVarHash[newAlmt].bp2-slop/2 <= other_bp <= SRVarHash[newAlmt].bp2+slop/2:
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                    elif detectIntraChrCopyInv and (SRVarHash[newAlmt].bp2-slop/2 <= other_bp <= SRVarHash[newAlmt].bp2+slop/2 or \
                        SRVarHash[newAlmt].bp2 < newAlmt[2] < other_bp or other_bp < newAlmt[2] < SRVarHash[newAlmt].bp2):

                        SRVarHash[newAlmt].typeSV = "INS_I"
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                        SRVarHash[newAlmt].bp3 = other_bp
                        SRVarHash[newAlmt].bp3tid = other_bp_tid
                        SRVarHash[newAlmt].n_changes+=1
                        #print "INV INS"
                elif SRVarHash[newAlmt].typeSV == "INS":
                    #print "INS"
                    if l_orient == r_orient and (SRVarHash[newAlmt].bp3 == -1 or SRVarHash[newAlmt].bp2-slop/2 \
                        <= other_bp <= SRVarHash[newAlmt].bp2+slop/2 or SRVarHash[newAlmt].bp3-slop/2 \
                        <= other_bp <= SRVarHash[newAlmt].bp3+slop/2):

                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                        if SRVarHash[newAlmt].bp3 == -1 and not (SRVarHash[newAlmt].bp2-slop/2 <= other_bp \
                            <= SRVarHash[newAlmt].bp2+slop/2 or SRVarHash[newAlmt].bp2 < newAlmt[2] < other_bp \
                            or other_bp < newAlmt[2] < SRVarHash[newAlmt].bp2):

                            SRVarHash[newAlmt].bp3 = other_bp
                            SRVarHash[newAlmt].bp3tid = other_bp_tid
                            SRVarHash[newAlmt].n_changes+=1
                #$ if 3rd bp absent write as inversion
                elif SRVarHash[newAlmt].typeSV == "INS_I" and l_orient != r_orient:
                    #print "INS_I"
                    if (SRVarHash[newAlmt].bp2-slop/2 <= other_bp <= SRVarHash[newAlmt].bp2+slop/2 or \
                        SRVarHash[newAlmt].bp3-slop/2 <= other_bp <= SRVarHash[newAlmt].bp3+slop/2):
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)

                    elif SRVarHash[newAlmt].bp3 == -1 and not (SRVarHash[newAlmt].bp2 < newAlmt[2] < other_bp \
                        or other_bp < newAlmt[2] < SRVarHash[newAlmt].bp2):
                        SRVarHash[newAlmt].count+=1
                        SRVarHash[newAlmt].support.append(SRFrag)
                        SRVarHash[newAlmt].bp3 = other_bp
                        SRVarHash[newAlmt].bp3tid = other_bp_tid
                        SRVarHash[newAlmt].n_changes+=1
            ## FORM NEW SR VARIANT
            else:
                newAlmt = (sr_bp1_tid, sr_bp2_tid, sr_bp1)
                newVariant = newSRVar()
                newVariant.bp2 = sr_bp2
                newVariant.support.append(SRFrag)
                newVariant.l_orient = l_orient
                newVariant.r_orient = r_orient
                if newAlmt[0] == newAlmt[1] and swap==0:
                    newVariant.swapped = 0
                    newVariant.typeSV = "DEL_INS"
                    if l_orient != r_orient:
                        newVariant.typeSV = "INV"
                        if minsr.reference_start < sr_bp1 and maxsr.reference_start < sr_bp2:
                            newVariant.orientStatus = 0
                        elif minsr.reference_end > sr_bp1 and maxsr.reference_end > sr_bp2:
                            newVariant.orientStatus = 1
                elif newAlmt[0] == newAlmt[1] and swap==1:
                    newVariant.typeSV = "TD_I"
                    newVariant.swapped = 1
                    #$handle this case in INS_I matches
                    if l_orient != r_orient:
                        newVariant.typeSV = "INV"
                        if minsr.reference_start < sr_bp1 and maxsr.reference_start < sr_bp2:
                            newVariant.orientStatus = 0
                        elif minsr.reference_end > sr_bp1 and maxsr.reference_end > sr_bp2:
                            newVariant.orientStatus = 1
                elif newAlmt[0] != newAlmt[1]:
                    newVariant.typeSV = "INS"
                    newVariant.swapped = 0
                    if l_orient != r_orient:
                        newVariant.typeSV = "INS_I"
                newVariant.isOriginal = 1
                if newAlmt not in SRVarHash:
                    SRVarHash[newAlmt] = newVariant

                newVariant2 = newSRVar()
                newVariant2.bp2 = sr_bp1
                newVariant2.support.append(SRFrag)
                newVariant2.l_orient = newVariant.r_orient
                newVariant2.r_orient = newVariant.l_orient
                newVariant2.typeSV = newVariant.typeSV
                newVariant2.swapped = newVariant.swapped
                newAlmtExt = (newAlmt[1], newAlmt[0], sr_bp2) 
                if newAlmtExt not in SRVarHash:
                    SRVarHash[newAlmtExt] = newVariant2

    ## WRITE REVISED PE VARIANTS TO FILE
    fAVN.write("VariantNum\tType\tchr1\tstart1\tstop1\tchr2\tstart2\tstop2\tchr3\tstart3\tstop3\tSupportBy\tNPEClusterSupp\tNFragPESupp\tNFragSRSupp\n")
    fAV.seek(0)
    header = fAV.readline()
    for lineVM in fVM:
        lineVM_split = lineVM.split()
        varNumPE = int(lineVM_split[0])
        for lineAV in fAV:
            lineAV_split = lineAV.split()
            if varNumPE in SRtoPESuppBPs:
                if len(SRtoPESuppFrags[varNumPE]) >= minSRtoPEsupport:
                    lineAV_split[11] = lineAV_split[11] + "_SR"
                    lineAV_split[3] = str(SRtoPESuppBPs[varNumPE][0])
                    lineAV_split[4] = str(SRtoPESuppBPs[varNumPE][0] + 1)
                if len(SRtoPESuppFrags[varNumPE]) >= minSRtoPEsupport and SRtoPESuppBPs[varNumPE][2] == -1:
                    if int(lineAV_split[6])  - slop < SRtoPESuppBPs[varNumPE][1] < int(lineAV_split[7]) + slop:
                        lineAV_split[6] = str(SRtoPESuppBPs[varNumPE][1])
                        lineAV_split[7] = str(SRtoPESuppBPs[varNumPE][1] + 1)
                    elif int(lineAV_split[9]) - slop < SRtoPESuppBPs[varNumPE][1] < int(lineAV_split[10]) + slop:
                        lineAV_split[9] = str(SRtoPESuppBPs[varNumPE][1])
                        lineAV_split[10] = str(SRtoPESuppBPs[varNumPE][1] + 1)
                    if (lineAV_split[1] == "DEL" or lineAV_split[1] == "TD" or lineAV_split[1] == "INV") \
                        and int(lineAV_split[3]) > int(lineAV_split[7]):
                        lineAV_split[3], lineAV_split[6] = lineAV_split[6], lineAV_split[3]
                        lineAV_split[4], lineAV_split[7] = lineAV_split[7], lineAV_split[4]
                    if lineAV_split[1] == "INV" and SRtoPESuppBPs[varNumPE][3] == 1:
                        lineAV_split[1] = "INV_B"
                # insertion matches
                elif len(SRtoPESuppFrags[varNumPE]) >= minSRtoPEsupport:
                    if SRtoPESuppBPs[varNumPE][1] != -1:
                        lineAV_split[6] = str(SRtoPESuppBPs[varNumPE][1])
                        lineAV_split[7] = str(SRtoPESuppBPs[varNumPE][1] + 1)
                    lineAV_split[9] = str(SRtoPESuppBPs[varNumPE][2])
                    lineAV_split[10] = str(SRtoPESuppBPs[varNumPE][2] + 1)

            break
        lineVMJ = "\t".join(lineVM_split)
        fVMN.write("%s" %lineVMJ)
        suppCount = 0
        if varNumPE in SRtoPESuppBPs:
            for SRFrag in SRtoPESuppFrags[varNumPE]:
                suppCount+=1
                fVMN.write("\t%s" %SRFrag)
        fVMN.write("\n")
        lineAV_split[14] = str(suppCount)
        lineAVJ = "\t".join(lineAV_split[:15])
        fAVN.write("%s\n" %lineAVJ)

    ## POSTPROCESS DE NOVO SR VARIANTS AND WRITE TO FILE
    k = 0
    for SRVar in SRVarHash:
        #print SRVar, SRVarHash[SRVar],SRVarHash[SRVar].typeSV, SRVarHash[SRVar].count
        if SRVarHash[SRVar].count > 0:
            bpMate = SRVarHash[SRVar].bp2
            chosenVar = (-1,-1,-1) 
            neighborSupport = []
            origSV = SRVar
            #of hash entry and mate entry for given SR, pick 1 and transfer support to 
            #one more "developed"
            if SRVarHash[SRVar].write == -1:

                    if SRVarHash[SRVar].n_changes > 0:
                        chosenVar = SRVar
                        SRVarHash[SRVar].write = 1
                    else:
                        SRVarHash[SRVar].write = 0
                        for ns in SRVarHash[SRVar].support[1:]:
                            neighborSupport.append(ns)
                    if SRVarHash[SRVar].isOriginal == 1:
                        origSV = SRVar

                    SRVarMate = (SRVar[1], SRVar[0], bpMate) 
                    if SRVarMate in SRVarHash and \
                        SRVarHash[SRVarMate].support[0] == SRVarHash[SRVar].support[0]:
                        if SRVarHash[SRVarMate].n_changes > 0 and chosenVar[2] == -1:
                            chosenVar = SRVarMate
                            SRVarHash[SRVarMate].write = 1
                        else:
                            SRVarHash[SRVarMate].write = 0
                            for ns in SRVarHash[SRVarMate].support[1:]:
                                neighborSupport.append(ns)
                        if SRVarHash[SRVarMate].isOriginal == 1:
                            origSV = SRVarMate

                    # if not chosen still, pick the original one among all neighbors
                    if chosenVar[2] == -1:
                        chosenVar = origSV
                        SRVarHash[chosenVar].write = 1
                    SRVarHash[chosenVar].support = \
                       SRVarHash[chosenVar].support + list(set(neighborSupport) - \
                       set(SRVarHash[chosenVar].support))
                    SRVarHash[chosenVar].count = len(SRVarHash[chosenVar].support)

            if SRVarHash[SRVar].write == 1 and SRVarHash[SRVar].count >= min_vs:
                k+=1
                if (SRVarHash[SRVar].typeSV == "INS_I" or SRVarHash[SRVar].typeSV == "INS") \
                    and SRVar[1] == SRVarHash[SRVar].bp3tid and abs(SRVarHash[SRVar].bp2 - \
                    SRVarHash[SRVar].bp3) < minSizeINS:
                    SRVarHash[SRVar].typeSV = "INS_POSS"
                elif (SRVarHash[SRVar].typeSV == "INS_I" or SRVarHash[SRVar].typeSV == "INS") \
                    and SRVar[1] == SRVarHash[SRVar].bp3tid and SRVarHash[SRVar].bp2 > \
                    SRVarHash[SRVar].bp3:
                    SRVarHash[SRVar].bp2, SRVarHash[SRVar].bp3 = SRVarHash[SRVar].bp3,\
                        SRVarHash[SRVar].bp2

                if (SRVarHash[SRVar].typeSV == "DEL_INS" or SRVarHash[SRVar].typeSV == "TD_I" \
                    or SRVarHash[SRVar].typeSV.startswith("INV")) and SRVar[2] > SRVarHash[SRVar].bp2:

                        output = [k+varNumPE+1, SRVarHash[SRVar].typeSV, SRVar[1], 
                                 SRVarHash[SRVar].bp2, SRVarHash[SRVar].bp2+1, SRVar[0], 
                                 SRVar[2], SRVar[2] + 1, SRVarHash[SRVar].bp3tid, SRVarHash[SRVar].bp3, 
                                 SRVarHash[SRVar].bp3 + 1, "SR", "0", "0", SRVarHash[SRVar].count]
                        outputN = map(str, output)
                        fAVN.write("%s\n" %("\t".join(outputN)))
                else:
                        output = [k+varNumPE+1, SRVarHash[SRVar].typeSV, SRVar[0], SRVar[2], 
                                 SRVar[2] + 1, SRVar[1], SRVarHash[SRVar].bp2, SRVarHash[SRVar].bp2 + 1,
                                 SRVarHash[SRVar].bp3tid, SRVarHash[SRVar].bp3, SRVarHash[SRVar].bp3
                                 + 1, "SR", "0", "0", SRVarHash[SRVar].count]
                        outputN = map(str, output)
                        fAVN.write("%s\n" %("\t".join(outputN)))
                fVMN.write("%s" %(k+varNumPE+1))
                for elem in SRVarHash[SRVar].support:
                    fVMN.write("\t%s" %elem)
                fVMN.write("\n")

    bamfile.close()

    #free memory
    SVHashPE.clear()
    del SVHashPE
    SRVarHash.clear()
    del SRVarHash
    chrHash.clear()
    del chrHash
    immutable_objects = []
    gc.collect()

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Add split reads to support existing PE variants and create new SR variants')
    PARSER.add_argument('workDir', help='Work directory')
    PARSER.add_argument('variantMapFilePE', help='File containing PEvariant map, typically variantMap.pe.txt')
    PARSER.add_argument('allVariantFilePE', help='File containing list of PE variants, typically allVariants.pe.txt')
    PARSER.add_argument('bamFileSR', help='File containing all split reads, name-sorted')
    PARSER.add_argument('-d', dest='debug', action='store_true',
                        help='print debug information')
    PARSER.add_argument('-x', action='store_true',
                                    help='no cluster cleanup used')
    PARSER.add_argument('-s', default=8.0, dest='slop', type=float, help='SR breakpoint slop')
    PARSER.add_argument('-f', default=0, dest='refRate', type=int, help='Subsample every so many split reads')
    PARSER.add_argument('-m', default=3, dest='min_vs', type=int,
        help='Minimum support for SR variants')
    PARSER.add_argument('-q', default=10, dest='mapThresh', type=int, help='SR Mapping quality threshold')
    PARSER.add_argument('-i', default=None, dest='ignoreBED',
        help='Exclude-regions file in BED format')
    PARSER.add_argument('-c', default=None, dest='ignoreChr',
        help='File listing chromosomes to exclude from analysis')
    PARSER.add_argument('-n', default=10, dest='minSizeINS', type=int,
        help='Minimum size for SR INS calls')
    PARSER.add_argument('-t', default=1, dest='minSRtoPEsupport', type=int,
        help='Minimum support for PE variants required to update breakpoints')
    ARGS = PARSER.parse_args()

    LEVEL = logging.INFO
    if ARGS.debug:
        LEVEL = logging.DEBUG

    logging.basicConfig(level=LEVEL,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    addSplitReads(ARGS.workDir, ARGS.variantMapFilePE, ARGS.allVariantFilePE,
                  ARGS.bamFileSR, ARGS.slop, ARGS.refRate, ARGS.min_vs,
                  ARGS.mapThresh, ARGS.ignoreChr, ARGS.minSizeINS,
                  ARGS.minSRtoPEsupport, ARGS.ignoreBED, ARGS.x)

    logging.shutdown()
