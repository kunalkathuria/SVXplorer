#!/usr/bin/env python
# Kunal Kathuria 1/18
# Filter variants based on local coverage in variant regions and write final SVs in BED format

import sys
import pysam
from collections import Counter
from bitarray import bitarray
import argparse

def formChrHash():
    print "Forming pile-up hash table..."
    fo=open(args.NH_REGIONS_FILE,"r")
    prev_start = -1
    prev_stop = -1
    for k,line in enumerate(fo):
        line_s = line.split()
        currentTID = line_s[0]
        start = int(line_s[1])
        stop = int(line_s[2])+1
        if currentTID not in chrHash:
            prev_stop = -1
            chrHash[currentTID] = bitarray()
            #bed file is 1-based
            for y in range(1,start):
                chrHash[currentTID].append(0)
        for x in range(start, stop):
            chrHash[currentTID].append(1)
        # make hash table of unreliable regions if greater than RDL (almt would be doubtful there)
        if prev_stop != -1 and currentTID == prevTID:
            addBit = 1
            if start - prev_stop > RDL:
                addBit = 0
                for x in range(prev_stop, start):
                    if currentTID not in chrHash:
                        chrHash[currentTID] = bitarray()
                    chrHash[currentTID].append(addBit)
        prev_start = start
        prev_stop = stop
        prevTID = currentTID
    print "Done"

def readBamStats():
    rdl, sd, coverage = -1, -1, -1
    for i,line in enumerate(fStat):
        if i == 0:
            rdl=float(line[:-1])
        elif i==2:
            sd=float(line[:-1])
        elif i==3:
            coverage = float(line[:-1])
            break
    return rdl, sd, coverage

def calculateDELThreshPE():
    PE_DEL_THRESH=PE_DEL_THRESH_S + int((SD-SD_S)*(PE_DEL_THRESH_L-PE_DEL_THRESH_S)/(SD_L-SD_S))
    if PE_DEL_THRESH > PE_DEL_THRESH_S:
        PE_DEL_THRESH = PE_DEL_THRESH_S
    elif PE_DEL_THRESH < PE_DEL_THRESH_L:
        PE_DEL_THRESH = UNIV_VAR_THRESH

    return PE_DEL_THRESH

def calculateLocCovg(chr_n, bpFirst, bpSecond):
    gap = bpSecond - bpFirst
    start = .25*gap + bpFirst
    stop = min(start+.5*gap,start +3*args.PILEUP_THRESH)
    covLoc, counter, confRegion = 0,0,0
    if stop > start:                    
        for pileupcolumn in fBAM.pileup(chr_n, start, stop):
            if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
                covLoc = covLoc + pileupcolumn.n
                counter+=1
                if counter > args.PILEUP_THRESH:
                    break
        if (counter > MIN_PILEUP_THRESH and (counter > GOOD_REG_THRESH*(stop-start) or counter > args.PILEUP_THRESH)):
            confRegion = 1
            covLoc = (1.0*covLoc)/(1.0*counter)
    return covLoc, confRegion

def updateBEDFiles():
    if lineAV_split[1] == "DEL":
        if lineAV_split[11].find("PE") == -1 and lineAV_split[11].find("SR") != -1 and \
            int(lineAV_split[7]) - int(lineAV_split[3]) < SR_DEL_THRESH:
            return 
        elif lineAV_split[11].find("SR") == -1 and lineAV_split[11].find("PE") != -1 and \
            int(lineAV_split[7]) - int(lineAV_split[3]) < PE_DEL_THRESH:
            return
        elif lineAV_split[11].find("SR") != -1 and lineAV_split[11].find("PE") != -1 and \
            int(lineAV_split[7]) - int(lineAV_split[3]) < MIX_DEL_THRESH:
            return
        fDEL.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(lineAV_split[2], lineAV_split[3], 
            lineAV_split[4], lineAV_split[5], lineAV_split[6], lineAV_split[7],"DEL",lineAV_split[11],GT,support))

    elif lineAV_split[1] == "TD" or (lineAV_split[1]=="TD_I" and lineAV_split[11].find("PE") != -1):
        [bp1_s, bp1_e] = min(int(lineAV_split[3]),int(lineAV_split[6])), min(int(lineAV_split[4]), int(lineAV_split[7]))
        [bp2_s, bp2_e] = max(int(lineAV_split[3]),int(lineAV_split[6])), max(int(lineAV_split[4]), int(lineAV_split[7]))

        fTD.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(lineAV_split[2], bp1_s, bp1_e, 
            lineAV_split[5], bp2_s, bp2_e,"TD",lineAV_split[11],GT,support))

    elif lineAV_split[1] == "INS_I" and lineAV_split[2] == lineAV_split[5] and lineAV_split[8] == "-1" \
        and lineAV_split[11][:2] == "PE":
    
        [bp1_s, bp1_e] = min(int(lineAV_split[3]),int(lineAV_split[6])), min(int(lineAV_split[4]), int(lineAV_split[7]))
        [bp2_s, bp2_e] = max(int(lineAV_split[3]),int(lineAV_split[6])), max(int(lineAV_split[4]), int(lineAV_split[7]))

        fTD.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(lineAV_split[2], bp1_s, bp1_e, lineAV_split[5], 
            bp2_s, bp2_e,"TD_INV",lineAV_split[11],GT,support))

    elif lineAV_split[1] == "INV" and ( (args.libINV and lineAV_split[11].find("SR") != -1) or \
        int(lineAV_split[12]) > 1): #or lineAV_split[1]=="INV_POSS":

        fINV.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(lineAV_split[2], lineAV_split[3], 
            lineAV_split[4], lineAV_split[5], lineAV_split[6], lineAV_split[7],lineAV_split[1],lineAV_split[11],GT,support))

    elif lineAV_split[1] == "Unknown" or lineAV_split[1] == "INS_POSS" or lineAV_split[1] == "TD_I" or \
        lineAV_split[1] == "INV_POSS" or ( (lineAV_split[1] == "INS" or lineAV_split[1] == "INS_I" or \
        lineAV_split[1] == "INS_C" or lineAV_split[1] == "INS_C_I" or lineAV_split[1]=="INS_C_P" or \
        lineAV_split[1]=="INS_C_I_P") and (lineAV_split[9] == "-1" or lineAV_split[6] == "-1" or bnd) ) \
        or lineAV_split[1] == "INS_C" or lineAV_split[1] == "INS_C_I":

        if lineAV_split[1] == "INS" or lineAV_split[1] == "INS_I" or lineAV_split[1] == "INS_C" or \
            lineAV_split[1] == "INS_C_I" or lineAV_split[1]=="INS_C_P" or lineAV_split[1]=="INS_C_I_P":

            fUnknown.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(lineAV_split[2], lineAV_split[3], 
                lineAV_split[7], lineAV_split[8], lineAV_split[9], lineAV_split[10],"BND", lineAV_split[1], 
                lineAV_split[11],GT,support))
        else:
            fUnknown.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(lineAV_split[2], lineAV_split[3], 
                lineAV_split[4], lineAV_split[5], lineAV_split[6], lineAV_split[7],"BND", lineAV_split[1], 
                lineAV_split[11],GT,support))

    elif len(lineAV_split[1]) > 2 and lineAV_split[1][0:3] == "INS" and not (lineAV_split[1]== "INS_C" or \
        lineAV_split[1] == "INS_C_I"):

        # two lines for insertion as in bedpe format; bp 1 and bp3 were flanks of bp2 by convention in INS_C classification 
        # unless confirmed further as INS_C_P
        [bp1_s, bp1_e]= int(lineAV_split[3]), int(lineAV_split[4])
        [bp2_s, bp2_e] = int(lineAV_split[6]),int(lineAV_split[7])
        [bp3_s, bp3_e] = int(lineAV_split[9]),int(lineAV_split[10])
        if swap:
            bp1_s, bp3_s = bp3_s, bp1_s
            bp1_e = bp3_e
        if bp3_s < bp2_s:
            bp2_s = bp3_s
            bp3_e = bp2_e
                
        fINS.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(lineAV_split[5], bp2_s, bp3_e, lineAV_split[2], 
            bp1_s, bp1_e, lineAV_split[1],GT,support))

    else:
        fUnknown.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(lineAV_split[2], lineAV_split[3], 
            lineAV_split[4], lineAV_split[5], lineAV_split[6], lineAV_split[7],"BND", lineAV_split[1],lineAV_split[11],GT,support))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Filter variants based on local coverage in variant regions \
        and write final SVs in BED format')
    parser.add_argument('workDir', help='Work directory')
    parser.add_argument('ufFile', help='Uniqueness filter output file, typically variants.uniqueFilter.txt')
    parser.add_argument('avFile', help='All variant file, typically allVariants._.txt')
    parser.add_argument('vmFile', help='Variant map file, typically variantMap._.txt')
    parser.add_argument('bamFile', help='Position-sorted args.BAM file for sample')
    parser.add_argument('statFile', help='File containing BAM statistics, typically bamStats.txt')
    parser.add_argument('-a', default=.6, dest='DEL_THRESH', type=float, 
        help='Local coverage threshold used for calling a deletion')
    parser.add_argument('-b', default=1.4, dest='DUP_THRESH', type=float,
        help='Local coverage threshold used for calling a duplication')
    parser.add_argument('-c', default='none', dest='NH_REGIONS_FILE',
        help='File containing non-homologous regions as chr, start, stop (provided in SVC)')
    parser.add_argument('-d', default=.8, dest='GOOD_REG_THRESH', type=float,
        help='Minimum ratio of unexcluded bases to to total bases covered by variant region')
    parser.add_argument('-e', default=1.4, dest='PILEUP_THRESH', type=float,
        help='Minimum unexcluded bases required to trust/use local coverage value in SV')
    parser.add_argument('-f', default=0, dest='splitINS', type=int, 
        help='Whether to split accidental complex variants into simple diploid variants based on local coverage')
    parser.add_argument('-g', default=0, dest='libINV', type=int,
        help='Liberal INV calling: 1 PE cluster + SR support sufficient')
    parser.add_argument('-v', default=0, dest='verbose', type=int, help='Verbose output')
    args = parser.parse_args()

    workDir= args.workDir
    # global constants
    DEL_THRESH2 = .125
    DEL_THRESH_L = .8
    DUP_THRESH_L = 1.2
    MIN_PILEUP_THRESH = 80
    # to trust pile-up depth in given region, this percentage of bases should return data
    GOOD_REG_THRESH=.8
    fAV = open(args.avFile,"r")
    fVM=open(args.vmFile,"r")
    fUF = open(args.ufFile,"r")
    fDEL = open(workDir+"/deletions.bedpe","w")
    fTD = open(workDir+"/tandemDuplications.bedpe","w")
    fINV = open(workDir+"/inversions.bedpe","w")
    fINS = open(workDir+"/insertions.bedpe","w")
    fUnknown = open(workDir+"/unknowns.bedpe","w")
    print "Writing final bedpe files using coverage information..."
    fStat = open(args.statFile,"r")
    RDL,SD,COVERAGE = readBamStats()
    if args.verbose:
        print "Coverage is:", COVERAGE

    # calculate min PE size based on insert length standard deviation under empirical model
    # aligner tends to mark concordants as discordants when SD is small unless distr v good, seemingly sharp effects    
    # around sig_IL of 20 and lower for Poisson and other related distributions.
    DEL_CONF_THRESH = .7
    UNIV_VAR_THRESH=100
    INS_VAR_THRESH=50
    SR_DEL_THRESH=80
    MIX_DEL_THRESH=50
    PE_DEL_THRESH_S=250
    PE_DEL_THRESH_L=150
    SD_S = 14
    SD_L = 24
    PE_DEL_THRESH = calculateDELThreshPE()

    uniqueFilterSVs = set()
    chrHash = {}
    for line in fUF:
        uniqueFilterSVs.add(int(line))
    fBAM = pysam.AlignmentFile(args.bamFile, "rb" )
    if args.NH_REGIONS_FILE != "none":
        formChrHash()

    for counter, lineAV in enumerate(fAV):
        support = 0
        for entry in fVM:
            support = len(entry.split()) - 1
            break
        counter+=1
        if args.verbose and counter % 100 == 0:
            print "Writing Variant", counter
        lineAV_split = lineAV.split()
        varNum = int(lineAV_split[0])

        if varNum in uniqueFilterSVs:
            ## skip SV if size thresholds not met
            if lineAV_split[1][0:3] == "TD" or (len(lineAV_split) > 2 and (lineAV_split[1][0:3] == "INV" \
                or lineAV_split[1][0:3] == "DEL")):
                if int(lineAV_split[7])-int(lineAV_split[3]) < UNIV_VAR_THRESH:
                    continue
            elif len(lineAV_split[1]) > 2 and lineAV_split[1][0:3] == "INS" and not (lineAV_split[1] == "INS_C" \
                or lineAV_split[1] == "INS_C_I"):
                if (0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH):
                    continue
            elif lineAV_split == "INS_C" or lineAV_split == "INS_C_I":
                if (0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH and 0 < \
                    int(lineAV_split[7])-int(lineAV_split[3]) < INS_VAR_THRESH) or (0 < int(lineAV_split[7])-\
                    int(lineAV_split[9]) < INS_VAR_THRESH and 0 < int(lineAV_split[4])-int(lineAV_split[6]) < INS_VAR_THRESH):
                    continue
  
            ## RUN PILEUP FILTER
            if lineAV_split[11].find("RD") == -1: 
                bnd=0
                swap = 0
                GT=""
                if (lineAV_split[1] == "DEL_INS" or lineAV_split[1] == "DEL" or lineAV_split[1][0:2]== "TD") and \
                    int(lineAV_split[4]) + MIN_PILEUP_THRESH < int(lineAV_split[6]):
                    confL, confR, covLocL, covLocR = 0,0,0,0
                    covLocM, confMiddle = calculateLocCovg(lineAV_split[2], int(lineAV_split[4]), int(lineAV_split[6]))

                    if lineAV_split[1][0:2] != "TD": 
                        verifL_start = int(lineAV_split[4])    
                        verifR_stop = int(lineAV_split[6])
                        VER_BUFFER = max(1.5*MIN_PILEUP_THRESH,.25*(verifR_stop-verifL_start))
                        verifL_stop = verifL_start + VER_BUFFER
                        verifR_start = verifR_stop - VER_BUFFER
                        if verifR_stop - verifL_start > VER_BUFFER:
                            covLocL, confL = calculateLocCovg(lineAV_split[2], verifL_start, verifL_stop)
                            covLocR, confR = calculateLocCovg(lineAV_split[2], verifR_start, verifR_stop)
                    if confMiddle == 1 or (confL == 1 and confR == 1):
                        if lineAV_split[1][0:2] == "TD" and ((confMiddle == 1 and covLocM/COVERAGE < args.DEL_THRESH) or \
                            (confL == 1 and covLocL/COVERAGE < DEL_CONF_THRESH and confR == 1 and \
                            covLocR/COVERAGE < DEL_CONF_THRESH)):
                       
                            if args.verbose:
                                print "TD confirmed (pileup)"   
                            lineAV_split[1] = "TD"

                        elif lineAV_split[1][0:2] == "TD" and confMiddle == 1:
                            if covLocM/COVERAGE < DEL_THRESH_L: #1.0:
                                lineAV_split[1] = "BND"
                
                        elif lineAV_split[1][0:3] == "DEL" and ((confMiddle == 1 and covLocM/COVERAGE < args.DEL_THRESH) \
                            or (confL == 1 and covLocL/COVERAGE < DEL_CONF_THRESH and confR == 1 and \
                            covLocR/COVERAGE < DEL_CONF_THRESH)):
                       
                            if args.verbose:
                                print "DEL confirmed (pileup)"
                            lineAV_split[1] = "DEL"
                            if covLocM/COVERAGE < DEL_THRESH2:
                                GT="GT:1/1"
                            elif covLocM/COVERAGE > 3*DEL_THRESH2:   
                                GT="GT:0/1"
                    
                        elif lineAV_split[1][0:3] == "DEL":
                            if confMiddle and covLocM/COVERAGE > DUP_THRESH_L: #1.0:
                                # since bp3 = -1, this will be written as a BND event
                                lineAV_split[1] = "INS"
                    
                elif len(lineAV_split[1]) > 2 and lineAV_split[11].find("PE") != -1 and (lineAV_split[1] == "INS" \
                    or lineAV_split[1] == "INS_I") and int(lineAV_split[7]) + MIN_PILEUP_THRESH < \
                    int(lineAV_split[9]) and lineAV_split[8] != "-1":

                    #bp2-3
                    covLoc, confReg = calculateLocCovg(lineAV_split[5], int(lineAV_split[7]), int(lineAV_split[9]))
                    if confReg and covLoc/COVERAGE < DUP_THRESH_L:
                        lineAV_split[1] = "INS_C_P"
                        if lineAV_split[1] == "INS_I":
                            lineAV_split[1] = "INS_C_I_P"
                    elif confReg and covLoc/COVERAGE < DEL_THRESH_L:
                        bnd = 1

                elif len(lineAV_split[1]) > 4 and lineAV_split[11].find("PE") != -1 and \
                    lineAV_split[1][0:5] == "INS_C" and lineAV_split[8] != "-1":
                    del_23 = 0
                    del_21 = 0
                    dup_23 = 0
                    dup_21 = 0
                    #bp2-3
                    start = int(lineAV_split[7])
                    stop = int(lineAV_split[9])
                    if start > stop:
                        start, stop = stop, start
                    covLoc_23, conf_23 = calculateLocCovg(lineAV_split[5], start, stop)
                    #bp1-2
                    start = int(lineAV_split[4])
                    stop = int(lineAV_split[6])
                    if start > stop:
                        start, stop = stop, start
                    if lineAV_split[2] == lineAV_split[5]:
                        covLoc_12, conf_12 = calculateLocCovg(lineAV_split[5], start, stop)
                    else:
                        convLoc_12, conf_12 = 0,0

                    if conf_23: 
                        if covLoc_23/COVERAGE > args.DUP_THRESH:
                            dup_23 = 1
                        elif covLoc_23/COVERAGE < args.DEL_THRESH:
                            del_23 =1
                    if conf_12: 
                        if covLoc_12/COVERAGE > args.DUP_THRESH:
                            dup_21 = 1
                        elif covLoc_12/COVERAGE < args.DEL_THRESH:
                            del_21 =1

                    confINSBP = 0
                    if (lineAV_split[1] == "INS_C" or lineAV_split[1] == "INS_C_I"): 
                        if dup_23:
                            lineAV_split[1]+="_P"
                            confINSBP = 1
                        elif dup_21:
                            lineAV_split[1]+="_P"
                            confINSBP = 1  
                            swap = 1
                        elif (del_23 or del_21) and not args.splitINS:
                            bnd = 1

            ## write in BED files
            updateBEDFiles()

    fAV.close()
    fVM.close()
    fUF.close()
    fDEL.close()
    fTD.close()
    fINV.close()
    fINS.close()
    fUnknown.close()
    fBAM.close()


