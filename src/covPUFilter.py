#!/usr/bin/env python

# Filter variants based on local coverage in variant regions and write final SVs in BED format
#PILEUP THRESH should be > 800 and < 10000

from sys import stderr
import pysam
from collections import Counter
import numpy as np
import argparse
import logging
from bitarray import bitarray
from shared import readChromosomeLengths, readBamStats

# global variables
DEL_THRESH_GT = .125
DUP_THRESH_S = 1.15
DEL_THRESH_S = .85
INS_COPY_THRESH = 1.1
INS_CUT_THRESH = .7
CALC_THRESH = 2000000
MQT_COV = 5
MQ0_Set = set()
# empirical calculation of DEL_THRESH b/w 15 and 25 stdev of IL
SD_S = 14
SD_L = 24
MIN_SPLIT_INS_COV = 7

def countCvg(fBAM, start, stop, chr_n):
    covOut = fBAM.count_coverage(chr_n, start, stop, read_callback="all", quality_threshold = MQT_COV)
    #covTotal = sum(map(covOut[0][1][:][,int))

def formChrHash(NH_REGIONS_FILE, RDL, chrLengths):
    global chrHash
    logging.info("Forming PU hash table...")
    fo=open(NH_REGIONS_FILE,"r")
    prev_start = -1
    prev_stop = -1
    for k,line in enumerate(fo):
        line_s = line.split()
        currentTID = line_s[0]
        start = int(line_s[1])
        stop = int(line_s[2])+1
        if currentTID not in chrHash:
            prev_stop = -1
            chrHash[currentTID] = bitarray(chrLengths[currentTID])
            chrHash[currentTID].setall(0)
            #bed file is 1-based
        # mark unreliable regions as 0 if greater than RDL (almt would be doubtful there)
        if prev_stop != -1 and currentTID == prevTID:
            if 0 < start - prev_stop <= RDL:
                chrHash[currentTID][prev_stop:start] = 1

        chrHash[currentTID][start:stop] = 1
        prev_start = start
        prev_stop = stop
        prevTID = currentTID
    logging.info("Done forming PU hash table")

def calculateLocCovg(NH_REGIONS_FILE,chr_n, bpFirst, bpSecond, PILEUP_THRESH, fBAM, chrHash, 
        GOOD_REG_THRESH, outerBPs, MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH, isTD):
    global covHash
    COV_REJ_FACTOR = 20
    MIN_BP_SPAN = 80
    bin_size_loc, largeDupBinThresh, TD_SIZE_SUSPECT_BOUND = 20,.5, 100000
    MIN_NBINS_LOC=40 #min((PILEUP_THRESH/bin_size_loc) - 10,40)
    bin_size = 100
    if chr_n not in covHash:
        logging.debug("Calculating coverage for %s", chr_n)
        counterBase, refLoop, cov_100bp, totalCov = 0,0,0,0
        MAX_ARRAY_SIZE = int(1.1*CALC_THRESH/bin_size)
        covList = np.empty((MAX_ARRAY_SIZE,))
        covListCounter = 0
        for pileupcolumn in fBAM.pileup(chr_n, stepper="all"):
            if NH_REGIONS_FILE is None or \
            (chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]):
                cov_100bp += pileupcolumn.n
                totalCov += pileupcolumn.n
                counterBase += 1
                refLoop += 1
                if refLoop == bin_size:
                    #logging.debug("Covg for chr %s in bin %s: %s", chr_n, 1 + len(covList), cov_100bp)  
                    covList[covListCounter] = 1.0*cov_100bp/refLoop
                    covListCounter+=1
                    cov_100bp, refLoop = 0,0
                if counterBase > CALC_THRESH:
                    break
        if counterBase > 0:
            #logging.debug("100bp cvg list: %s", covList)
            avgCov = 1.0*totalCov/counterBase
            covList = covList[0:covListCounter]
            covList.sort()
            covHash[chr_n] = covList[covList.size/2] #avgCov
            #change to debug when test done
            logging.debug("Median coverage of Chr %s written as %f; average was %f",
                          chr_n, covHash[chr_n], avgCov)
        else:
            logging.debug("Unable to calculate coverage reliably in chromosome %s", chr_n)
            print >> stderr, "Note: unable to calculate coverage reliably in chromosome", chr_n
            covHash[chr_n] = 0

    if bpSecond - bpFirst < MIN_BP_SPAN:
        bpFirstL = outerBPs[0]
        bpSecondL = outerBPs[1]
    else:
        bpFirstL = bpFirst
        bpSecondL = bpSecond

    covList = np.empty((1,))
    gap = bpSecondL - bpFirstL
    start = .25*gap + bpFirstL
    stop = min(start+.5*gap,start +3*PILEUP_THRESH)
    covLoc,covLocNH, counter, counterNH, confRegion, gbCount, bCountNH, covLocG, \
        covBinLoc, refLoopLoc, covListLocCounter = 0,0,0,0,0,0,0,0,0,0,0
    MAX_TD_SIZE = 1500000
    MAX_PU_SIZE = 100000
    MAX_ARRAY_SIZE = int(1.1*MAX_PU_SIZE/bin_size_loc)
    covListLoc = np.empty((MAX_ARRAY_SIZE,))
    rejRegion = 0
    if stop > start and ((not isTD) or (stop - start < MAX_TD_SIZE)):
        for pileupcolumn in fBAM.pileup(chr_n, start, stop, stepper="all", truncate=True):
            puVal = pileupcolumn.n
            #if region is too busy, do not trust it at all
            if puVal > COV_REJ_FACTOR * covHash[chr_n]:
                covLoc, covLocNH, counter, counterNH = 0,0,0,0
                rejRegion = 1
                break
            covLocNH+= puVal
            counterNH+=1
            
            badReadCount, allBadReads = 0,0
            for read in pileupcolumn.pileups:
                if read.alignment.mapping_quality == 0:
                    covLocNH-= 1
                    badReadCount+=1
            if badReadCount == len(pileupcolumn.pileups) and badReadCount != 0:
                allBadReads = 1
                counterNH-=1
                MQ0_Set.add((chr_n,pileupcolumn.pos))

            if NH_REGIONS_FILE is None and counterNH > PILEUP_THRESH:
                break
            logging.debug("PU pos, val: %s %s", pileupcolumn.pos, puVal)

            if (NH_REGIONS_FILE is not None) and \
            (chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]):
                logging.debug("PUval, pu.n: %s, %s", puVal, pileupcolumn.n)
                covLoc+= puVal
                counter+=1
                covLoc-= badReadCount
                if allBadReads:
                    counter-= 1
                    
                covBinLoc += puVal
                covBinLoc-= badReadCount
                refLoopLoc += 1
                if allBadReads:
                    refLoopLoc-=1

                if isTD == 1 and refLoopLoc == bin_size_loc and \
                    bpSecond - bpFirst > TD_SIZE_SUSPECT_BOUND:
                    #logging.debug("Covg for chr %s in bin %s: %s", chr_n, 1 + len(covListLoc), cov_100bp)  
                    logging.debug("Appending %s to covBinLoc", 1.0*covBinLoc/refLoopLoc)
                    covListLoc[covListLocCounter] = 1.0*covBinLoc/refLoopLoc
                    covListLocCounter+=1
                    covBinLoc, refLoopLoc = 0,0

                if counter > PILEUP_THRESH:
                    break
        logging.debug("CounterNH, covLocNH, rejRegion 1: %s, %s, %s", counterNH, covLocNH, rejRegion)

        # pile-up does not step through loop if no reads lie at any point in variant region
        # we need to account for these bases nevertheless in calculating true coverage
        counterNH = 0
        counter = 0
        if rejRegion == 0:
            for x in range(int(start), int(stop)):
                if (NH_REGIONS_FILE is not None and \
                    (chr_n in chrHash and x < len(chrHash[chr_n]) and chrHash[chr_n][x])) and \
                    (chr_n,x) not in MQ0_Set:
                    counter+=1
                counterNH+=1
                if NH_REGIONS_FILE is None and counterNH > PILEUP_THRESH:
                    break
                if counter > PILEUP_THRESH:
                    break
        rejRegion = 0

        # add sides also 
        gbCount, bCountNH = counter, counterNH
        if bpSecond - bpFirst >= MIN_BP_SPAN:
            counter, counterNH,covBinLoc, refLoopLoc = 0,0,0,0
            start, stop = bpFirstL, bpFirstL + .25*gap
            logging.debug("Start, stop 2:%s, %s", start, stop)
            for pileupcolumn in fBAM.pileup(chr_n, start, stop, stepper="all", truncate=True):
                # we wish to avoid adding coverage from boundaries, especially if it's not a good region 
                # but also need more stats/evidence if counterNH doesn't have many bases from good region
                puVal = pileupcolumn.n
                #if region is too busy, do not trust it at all
                if puVal > COV_REJ_FACTOR * covHash[chr_n]:
                    covLoc, covLocNH, counter, counterNH = 0,0,0,0
                    rejRegion = 1
                    break
                covLocNH+= puVal
                counterNH+=1
                
                badReadCount, allBadReads = 0,0 
                for read in pileupcolumn.pileups:
                    if read.alignment.mapping_quality == 0:
                        covLocNH-= 1
                        badReadCount+=1
                if badReadCount == len(pileupcolumn.pileups) and badReadCount != 0:
                    allBadReads = 1
                    counterNH-=1
                    MQ0_Set.add((chr_n,pileupcolumn.pos))

                if NH_REGIONS_FILE is None and counterNH > PILEUP_THRESH:
                    break

                if NH_REGIONS_FILE is not None and \
                (chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]):
                    covLoc+= puVal
                    counter+=1
                    covLoc-= badReadCount
                    if allBadReads:
                        counter-= 1
                    
                    covBinLoc += puVal
                    covBinLoc-= badReadCount
                    refLoopLoc += 1
                    if allBadReads:
                        refLoopLoc-=1

                    if isTD == 1 and refLoopLoc == bin_size_loc and \
                        bpSecond - bpFirst > TD_SIZE_SUSPECT_BOUND:
                        logging.debug("Appending %s to covBinLoc", 1.0*covBinLoc/refLoopLoc)
                        covListLoc[covListLocCounter] = 1.0*covBinLoc/refLoopLoc
                        covListLocCounter+=1
                        covBinLoc, refLoopLoc = 0,0
                    if counter > PILEUP_THRESH:
                        break

            logging.debug("counterNH, covLocNH 2a: %s, %s", counterNH, covLocNH)
            # pile-up does not step through loop if no portion of reads lies in variant region
            counter, counterNH = 0,0
            if rejRegion == 0:
                for x in range(int(start), int(stop)):
                    if (NH_REGIONS_FILE is not None and \
                        (chr_n in chrHash and x < len(chrHash[chr_n]) and chrHash[chr_n][x])) and \
                        (chr_n,x) not in MQ0_Set:
                        counter+=1         
                    counterNH+=1
                    if NH_REGIONS_FILE is None and counterNH > PILEUP_THRESH:
                        break
                    if counter > PILEUP_THRESH:
                        break
            rejRegion = 0

            logging.debug("counterNH, covLocNH 2b: %s, %s", counterNH, covLocNH)
            gbCount+= counter
            bCountNH+= counterNH

            counter, counterNH, covBinLoc, refLoopLoc = 0,0,0,0
            start, stop = bpSecondL - .25*gap, bpSecondL
            logging.debug("Start, stop 3:%s, %s", start, stop)
            for pileupcolumn in fBAM.pileup(chr_n, start, stop, stepper="all", truncate=True):
                puVal = pileupcolumn.n
                #if region is too busy, do not trust it at all
                if puVal > COV_REJ_FACTOR * covHash[chr_n]:
                    covLoc, covLocNH, counter, counterNH = 0,0,0,0
                    rejRegion = 1
                    break
                covLocNH+=puVal
                counterNH+=1
                
                badReadCount, allBadReads = 0,0
                for read in pileupcolumn.pileups:
                    if read.alignment.mapping_quality == 0:
                        covLocNH-= 1
                        badReadCount+=1
                if badReadCount == len(pileupcolumn.pileups) and badReadCount != 0:
                    allBadReads = 1
                    counterNH-=1
                    MQ0_Set.add((chr_n,pileupcolumn.pos))

                if NH_REGIONS_FILE is None and counterNH > PILEUP_THRESH:
                    break

                if NH_REGIONS_FILE is not None and \
                (chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]):
                    covLoc+=puVal
                    counter+=1
                    covLoc-= badReadCount
                    if allBadReads:
                        counter-= 1

                    covBinLoc += puVal
                    covBinLoc-= badReadCount
                    refLoopLoc += 1
                    if allBadReads:
                        refLoopLoc-=1
                
                    if isTD == 1 and refLoopLoc == bin_size_loc and \
                        bpSecond - bpFirst > TD_SIZE_SUSPECT_BOUND:
                        logging.debug("Appending %s to covBinLoc", 1.0*covBinLoc/refLoopLoc)
                        covListLoc[covListLocCounter] = 1.0*covBinLoc/refLoopLoc
                        covListLocCounter+=1
                        covBinLoc, refLoopLoc = 0,0
                    if counter > PILEUP_THRESH:
                        break
            logging.debug("counterNH, covLocNH 3a: %s, %s", counterNH, covLocNH)
            # pile-up does not step through loop if no portion of reads lies in variant region
            counter, counterNH = 0,0
            if rejRegion == 0:
                for x in range(int(start), int(stop)):
                    if (NH_REGIONS_FILE is not None and \
                        (chr_n in chrHash and x < len(chrHash[chr_n]) and chrHash[chr_n][x])) and \
                        (chr_n,x) not in MQ0_Set:
                        counter+=1         
                    counterNH+=1
                    if NH_REGIONS_FILE is None and counterNH > PILEUP_THRESH:
                        break
                    if counter > PILEUP_THRESH:
                        break
            gbCount+= counter
            bCountNH+= counterNH
            logging.debug("counterNH, covLocNH 3b: %s, %s", counterNH, covLocNH)

        logging.debug("bCountNH final: %s", bCountNH)
        covLocG = covLoc
        if NH_REGIONS_FILE is not None and gbCount >= MIN_PILEUP_THRESH:
            confRegion = 1
            covLoc = (1.0*covLoc)/(1.0*gbCount)
            logging.debug("Good Reg used. confRegion, covLoc: %s, %s", confRegion, covLoc)
        elif bCountNH >= MIN_PILEUP_THRESH_NH:
            confRegion = 1
            covLoc = (1.0*covLocNH)/(1.0*bCountNH)
            logging.debug("NH used. confNH, covLocNH: %s, %s", confRegion, covLoc)

    logging.debug("NH bases found, covNH final: %s, %s", bCountNH, covLocNH/(.01+bCountNH))
    logging.debug("Good bases found, confirmation, covGood final, chr cov: %s, %s, %s, %s", gbCount, confRegion, covLocG/(.01+gbCount), covHash[chr_n])
    covListLoc = covListLoc[0:covListLocCounter]
    if covHash[chr_n] != 0:
        largeDupRet, threshIndex = 0,0
        if isTD == 1 and bpSecond - bpFirst > TD_SIZE_SUSPECT_BOUND and \
            covListLoc.size > MIN_NBINS_LOC:
            covListLoc.sort()
            for k,val in enumerate(covListLoc):
                if 1.*val/covHash[chr_n] > DUP_THRESH_S:
                    break
            threshIndex = k
            largeDupRet = 1
            if 1.*threshIndex/covListLoc.size < largeDupBinThresh:
                largeDupRet = 2
        logging.debug("largeDupRet, threshIndex, lenCovListLoc: %s, %s, %s", largeDupRet, threshIndex, covListLoc.size)
        return 1.0*covLoc/covHash[chr_n], confRegion, largeDupRet
    else:
        return 0, 0,0

def writeVariants(lineAV_split, swap, bnd, support, GT, fAVN, 
                  SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo):
   
    svtype = lineAV_split[1]
    lineAV_split_T = list(lineAV_split)
    if svtype.startswith("TD") or svtype.startswith("INV") or svtype.startswith("DEL") or\
        (svtype.startswith("INS_half") and lineAV_split[2] == lineAV_split[5]):
        if int(lineAV_split[7])-int(lineAV_split[3]) < UNIV_VAR_THRESH:
            return

    if svtype == "DEL":
        if lineAV_split[11].find("PE") == -1 and lineAV_split[11].find("SR") != -1 and \
            int(lineAV_split[7]) - int(lineAV_split[3]) < SR_DEL_THRESH:
            return
        elif lineAV_split[11].find("SR") != -1 and lineAV_split[11].find("PE") != -1 and \
            int(lineAV_split[7]) - int(lineAV_split[3]) < MIX_DEL_THRESH:
            return

    # correct SPLIT_INS fields just in case
    if  svtype.startswith("TD") or  svtype.startswith("DEL") or  svtype.startswith("BND"):
        lineAV_split_T[8:11] = ["-1", "-1", "-1"]
    lineAV_split_T.extend([str(swap),str(bnd),str(support),GT, str(covInfo)])
    fAVN.write("%s\n" %("\t".join(lineAV_split_T)))

def covPUFilter(workDir, avFile, vmFile, ufFile, statFile, bamFile,
                NH_REGIONS_FILE, DEL_THRESH, DUP_THRESH, splitINS, 
                PILEUP_THRESH, GOOD_REG_THRESH, minVariantSize):

    #statistical constants 
    MIN_PILEUP_THRESH = 80
    MIN_PILEUP_THRESH_NH = 80
    UNIV_VAR_THRESH = minVariantSize
    INS_VAR_THRESH = minVariantSize
    SR_DEL_THRESH = max(80, minVariantSize)
    MIX_DEL_THRESH = max(50, minVariantSize)
    fAV = open(avFile,"r")
    fVM=open(vmFile,"r")
    fUF = open(ufFile,"r")
    fAVN = open(workDir+"/allVariants.pu.txt","w")
    fAVN.write("VariantNum\tType\tchr1\tstart1\tstop1\tchr2\tstart2\tstop2\tchr3\tstart3\tstop3\t SupportBy\tNPEClusterSupp\tNFragPESupp\tNFragSRSupp\tSwapBP\tBNDFlag\tSupport\tGT\n")
    logging.info("Writing final bedpe files using coverage information")
    RDL,SD,COV = readBamStats(statFile)
    logging.info("Some stats from BAM. RDL: %d, Sd: %f", 
                  RDL, SD)

    # calculate min PE size based on insert length standard deviation under empirical model
    # aligner tends to mark concordants as discordants when SD is small unless distr v good, seemingly sharp effects
    # around sig_IL of 20 and lower for generally well-behaved distributions. Without rigor, though.
    logging.info("DEL_THRESH, DEL_THRESH_S, DUP_THRESH, DUP_THRESH_S, MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH: %s, %s, %s, %s, %s, %s", DEL_THRESH, DEL_THRESH_S, DUP_THRESH, DUP_THRESH_S, MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH)

    uniqueFilterSVs = set()
    global chrHash
    global covHash
    chrHash = {}
    covHash = {}

    for line in fUF:
        uniqueFilterSVs.add(int(line))

    fBAM = pysam.AlignmentFile(bamFile, "rb" )
    if NH_REGIONS_FILE is not None:
        logging.info("Using good-regions BED file %s in cov PU", NH_REGIONS_FILE)
        chrLengths = readChromosomeLengths(bamFile)
        formChrHash(NH_REGIONS_FILE, RDL, chrLengths)
    else:
        print >> stderr, "Warning! Not using a good regions file for pile-up filter! This can affect some coverage-based results adversely."

    svtype = None
    header = fAV.readline()

    for counter, lineAV in enumerate(fAV):
        if counter % 250 == 0:
            logging.info("Variant %s written.", counter)
        counter+=1
        lineAV_split = lineAV.split()
        varNum = int(lineAV_split[0])
        support, isTD = 0,0
        if lineAV_split[13] != "." and lineAV_split[14] != ".":
            support = int(lineAV_split[13]) + int(lineAV_split[14])
        elif lineAV_split[13] != ".":
            support = int(lineAV_split[13])
        elif lineAV_split[14] != ".":
            support = int(lineAV_split[14])
        if varNum in uniqueFilterSVs:
            svtype = lineAV_split[1]

            ## skip SV if size thresholds not met
            if svtype.startswith("TD") or svtype.startswith("INV") or svtype.startswith("DEL") or\
                (svtype.startswith("INS_half") and lineAV_split[2] == lineAV_split[5]):
                if int(lineAV_split[7])-int(lineAV_split[3]) < UNIV_VAR_THRESH:
                    continue

            ## RUN PILEUP FILTER
            if lineAV_split[11].find("RD") == -1:
                bnd=0
                swap = 0
                GT="."
                covInfo = -1
                if (svtype == "DEL_INS" or svtype == "DEL" or svtype[0:2]== "TD") and \
                    int(lineAV_split[3]) + min(MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH) < int(lineAV_split[7]):

                    if svtype.startswith("TD"):
                        isTD = 1
                    covLocM, confMiddle, largeDupRet = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[2],
                            int(lineAV_split[4]), int(lineAV_split[6]),
                            PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, [int(lineAV_split[3]), int(lineAV_split[7])],
                            MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH, isTD)
                    covInfo = covLocM
                    if svtype.startswith("TD") and lineAV_split[11].find("PE") == -1 and \
                        confMiddle == 1 and ((largeDupRet == 0 and covLocM >= DUP_THRESH_S) \
                        or largeDupRet == 2): 

                        logging.debug("TD confirmed using pileup")
                        svtype = "TD"

                    elif svtype.startswith("TD") and (confMiddle == 0 or \
                        (largeDupRet == 0 and covLocM < DUP_THRESH_S) or largeDupRet == 1):
                        logging.debug("Rejected %s due to low cvg or insufficient evidence of coverage due to small variant region", lineAV)
                        # since bp3 = -1, this will be written as a BND event
                        svtype = "INS_halfRF"

                    elif svtype.startswith("DEL") and confMiddle == 1 and covLocM <= DEL_THRESH: 

                        logging.debug("DEL confirmed using pileup")
                        svtype = "DEL"
                        if covLocM < DEL_THRESH_GT:
                            GT="GT:1/1"
                        elif covLocM > 3*DEL_THRESH_GT:
                            GT="GT:0/1"

                    elif svtype.startswith("DEL") and ((lineAV_split[11].find("PE") != -1 and \
                        (confMiddle == 0 or covLocM > DEL_THRESH_S)) or \
                        (lineAV_split[11].find("PE") == -1 and (confMiddle == 0 or covLocM > DEL_THRESH))):
                        logging.debug("Rejected %s due to high cvg or insufficient evidence of coverage due to small variant region", lineAV)
                        # since bp3 = -1, this will be written as a BND event
                        svtype = "INS_halfFR"
                        
                elif (svtype in ["INS", "INS_I"] or svtype.startswith("INS_C")) and lineAV_split[8] != "-1":

                    del_23, del_12, del_13, dup_23, dup_12, dup_13 = 0, 0, 0, 0, 0, 0
                    swap_12, swap_13 = 0,0
                    #bp2-3
                    #use inner bounds for all
                    start = int(lineAV_split[7])
                    stop = int(lineAV_split[9])
                    covLoc_23, conf_23, _ = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                            start, stop, PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, 
                            [int(lineAV_split[6]), int(lineAV_split[10])], MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH, isTD)
                    covInfo = covLoc_23

                    #bp1-2 ("12" will here refer to 1, the paste bp, and the closest of the other 2)
                    start = int(lineAV_split[4])
                    stop = int(lineAV_split[6])
                    outerBPs = [int(lineAV_split[3]), int(lineAV_split[7])]
                    if start > stop:
                        swap_12 = 1
                        start = int(lineAV_split[10])
                        stop = int(lineAV_split[3])
                        outerBPs = [int(lineAV_split[9]), int(lineAV_split[4])]
                    if lineAV_split[2] == lineAV_split[5]:
                        covLoc_12, conf_12, _ = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                                start, stop, PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, outerBPs, 
                                MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH, isTD)
                    else:
                        covLoc_12, conf_12 = 0,0
                   
                    #bp1-3 ("13" will here refer to 1, the paste bp, and the farther of the other 2)
                    covLoc_13, conf_13 = 0,0
                    if splitINS == True:
                        start = int(lineAV_split[4])
                        stop = int(lineAV_split[9])
                        outerBPs = [int(lineAV_split[3]), int(lineAV_split[10])]
                        if start > stop:
                            swap_13 = 1
                            start = int(lineAV_split[7])
                            stop = int(lineAV_split[3])
                            outerBPs = [int(lineAV_split[6]), int(lineAV_split[4])]
                        if lineAV_split[2] == lineAV_split[5]:
                            covLoc_13, conf_13, _ = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                                    start, stop, PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, outerBPs, 
                                    MIN_PILEUP_THRESH, MIN_PILEUP_THRESH_NH, isTD)

                    if conf_23:
                        if covLoc_23 > DUP_THRESH:
                            dup_23 = 1
                        elif covLoc_23 < DEL_THRESH:
                            del_23 =1
                    if conf_12:
                        if covLoc_12 > DUP_THRESH:
                            dup_12 = 1
                        elif covLoc_12 < DEL_THRESH:
                            del_12 =1
                    if conf_13:
                        if covLoc_13 > DUP_THRESH:
                            dup_13 = 1
                        elif covLoc_13 < DEL_THRESH:
                            del_13 =1

                    lineAV_split1 = list(lineAV_split)
                    if svtype in ["INS", "INS_I"] and lineAV_split[8] != "-1":

                        # bp2-3
                        if (0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH) or \
                            (conf_23 and covLoc_23 < INS_COPY_THRESH): 
                            bnd = 1
                            
                    elif svtype.startswith("INS_C") and lineAV_split[11].find("PE") != -1 and \
                        lineAV_split[8] != "-1":

                        if splitINS == True and COV > MIN_SPLIT_INS_COV:    
                            logging.debug("Split INS is true -- INS_C: %s", lineAV)
                            
                            if lineAV_split[2] == lineAV_split[5] and del_12 and del_23:
                                logging.debug("Split INS_C 1")
                                #1-2 is del
                                lineAV_split1[1] = "DEL"
                                if swap_12:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[9], lineAV_split[10], \
                                                        lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)
                                #2-3 is DEL
                                lineAV_split1[1] = "DEL"
                                lineAV_split1[2:8] = lineAV_split[5], lineAV_split[6], lineAV_split[7], \
                                                    lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)

                                #1-3 is bnd
                                if int(lineAV_split[12]) >= 3:
                                    lineAV_split1[1] = "BND"
                                    if swap_13:
                                        lineAV_split1[2:8] = lineAV_split[2], lineAV_split[6], lineAV_split[7], \
                                                            lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                    else:
                                        lineAV_split1[2:8] = lineAV_split[2], lineAV_split[3], lineAV_split[4], \
                                                            lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                    writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                                  SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)
                                continue
                            
                            elif lineAV_split[2] == lineAV_split[5] and del_23 and dup_13:
                                logging.debug("Split INS_C 3")
                                #1-2 is bnd
                                lineAV_split1[1] = "BND"
                                if swap_12:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[9], lineAV_split[10], \
                                                        lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)
                                #1-3 is TD
                                lineAV_split[1] = "TD"
                                if swap_13:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[6], lineAV_split[7], \
                                                        lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                else:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[3], lineAV_split[4], \
                                                        lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)

                                #2-3 is del
                                lineAV_split1[1] = "DEL"
                                lineAV_split1[2:8] = lineAV_split[5], lineAV_split[6], lineAV_split[7], \
                                                    lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)
                                continue
                            
                            elif del_23:
                                logging.debug("Split INS_C 5")
                                #1-2 is bnd
                                lineAV_split1[1] = "BND"
                                if swap_12:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[9], lineAV_split[10], \
                                                        lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)
                                #1-3 is bnd
                                if int(lineAV_split[12]) >= 3:
                                    lineAV_split1[1] = "BND"
                                    if swap_13:
                                        lineAV_split1[2:8] = lineAV_split[2], lineAV_split[6], lineAV_split[7], \
                                                            lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                    else:
                                        lineAV_split1[2:8] = lineAV_split[2], lineAV_split[3], lineAV_split[4], \
                                                            lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                    writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                                  SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)

                                #2-3 is del
                                lineAV_split1[1] = "DEL"
                                lineAV_split1[2:8] = lineAV_split[5], lineAV_split[6], lineAV_split[7], \
                                                    lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, 
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)
                                continue

                        if (svtype == "INS_C" or svtype == "INS_C_I"):
                            if dup_23 and not dup_12:
                                svtype+="_P"
                            elif dup_12 and not dup_23:
                                svtype+="_P"
                                swap = 1
                
                        #NPE_CLUSTERS_SUPP = int(lineAV_split[12])
                        if (svtype == "INS_C_P" or svtype == "INS_C_I_P") and \
                            ((0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH) or \
                            (conf_23 and covLoc_23 < INS_CUT_THRESH)):
                            # or (conf_12 and not (DEL_THRESH < covLoc_12 < DUP_THRESH))
                            # or NPE_CLUSTERS < 3
                            bnd = 1

            ## write in BED files
            lineAV_split[1] = svtype
            writeVariants(lineAV_split, swap, bnd, support, GT, fAVN,
                    SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH, covInfo)

    fAV.close()
    fAVN.close()
    fVM.close()
    fUF.close()
    fBAM.close()
    chrHash.clear()
    del chrHash
    covHash.clear()
    del covHash

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Filter variants based on local coverage in variant regions and write final SVs in BED format')
    PARSER.add_argument('workDir', help='Work directory')
    PARSER.add_argument('ufFile', help='Uniqueness filter output file, typically variants.uniqueFilter.txt')
    PARSER.add_argument('avFile', help='All variant file, typically allVariants._.txt')
    PARSER.add_argument('vmFile', help='Variant map file, typically variantMap._.txt')
    PARSER.add_argument('bamFile', help='Position-sorted args.BAM file for sample')
    PARSER.add_argument('statFile', help='File containing BAM statistics, typically bamStats.txt')
    PARSER.add_argument('-a', default=.6, dest='DEL_THRESH', type=float,
        help='Local coverage threshold used for calling a deletion')
    PARSER.add_argument('-b', default=1.4, dest='DUP_THRESH', type=float,
        help='Local coverage threshold used for calling a duplication')
    PARSER.add_argument('-c', default=None, dest='NH_REGIONS_FILE',
        help='File containing non-homologous regions as chr, start, stop')
    PARSER.add_argument('-d', action='store_true', dest='debug', 
        help='print debug information')
    PARSER.add_argument('-i', default=.8, dest='GOOD_REG_THRESH', type=float,
        help='Minimum ratio of unexcluded bases to to total bases covered by variant region')
    PARSER.add_argument('-e', default=500.0, dest='PILEUP_THRESH', type=float,
        help='Minimum unexcluded bases required to trust/use local coverage value in SV')
    PARSER.add_argument('-f', default=0, dest='splitINS', type=int,
        help='Whether to split accidental complex variants into simple diploid variants based on local coverage')
    PARSER.add_argument('-g', default=1, dest='libINV', type=int,
        help='Liberal INV calling: 1 PE cluster + SR support sufficient')
    PARSER.add_argument('-s', default=100, dest='minVarSize', type=int, help='minimum size of variants called')
    PARSER.add_argument('-v', default=0, dest='verbose', type=int, help='Verbose output')
    ARGS = PARSER.parse_args()

    LEVEL = logging.INFO
    if ARGS.debug:
        LEVEL = logging.DEBUG

    logging.basicConfig(level=LEVEL,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    covPUFilter(ARGS.workDir, ARGS.avFile, ARGS.vmFile, ARGS.ufFile, 
                ARGS.statFile, ARGS.bamFile, ARGS.NH_REGIONS_FILE,
                ARGS.DEL_THRESH, ARGS.DUP_THRESH, ARGS.splitINS,
                ARGS.PILEUP_THRESH, ARGS.GOOD_REG_THRESH, ARGS.minVarSize)

    logging.shutdown()
