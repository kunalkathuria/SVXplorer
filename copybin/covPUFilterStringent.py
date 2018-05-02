#!/usr/bin/env python

# Filter variants based on local coverage in variant regions and write final SVs in BED format

from sys import stderr
import pysam
from collections import Counter
from bitarray import bitarray
import argparse
import logging

# global variables
DEL_THRESH2 = .125
DEL_THRESH_L = .8
DUP_THRESH_L = 1.2
DEL_THRESH_H = .25
DUP_THRESH_H = 1.75
MIN_PILEUP_THRESH = 80
CALC_THRESH = 1000000
PE_DEL_THRESH_S = 102
PE_DEL_THRESH_L = 100
chrHash = {}
covHash = {}
# empirical calculation of DEL_THRESH b/w 15 and 25 stdev of IL
SD_S = 14
SD_L = 24
MIN_SPLIT_INS_COV = 7

def formChrHash(NH_REGIONS_FILE, RDL):
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
            chrHash[currentTID] = bitarray()
            #bed file is 1-based
            for y in range(1,start):
                chrHash[currentTID].append(0)
        # make hash table of unreliable regions if greater than RDL (almt would be doubtful there)
        if prev_stop != -1 and currentTID == prevTID:
            addBit = 1
            if start - prev_stop > RDL:
                addBit = 0
            for x in range(prev_stop, start):
                if currentTID not in chrHash:
                    chrHash[currentTID] = bitarray()
                chrHash[currentTID].append(addBit)
        for x in range(start, stop):
            chrHash[currentTID].append(1)
        prev_start = start
        prev_stop = stop
        prevTID = currentTID
    logging.info("Done forming PU hash table")

def readBamStats(statFile):
    rdl, sd, coverage = -1, -1, -1
    with open(statFile, 'r') as fStat:
        for i,line in enumerate(fStat):
            if i == 0:
                rdl = float(line[:-1])
            elif i == 2:
                sd = float(line[:-1])
            elif i == 3:
                coverage = float(line[:-1])
                break
    return rdl, sd, coverage

def calculateDELThreshPE(SD):
    PE_DEL_LOWEST = 100
    PE_DEL_THRESH=PE_DEL_THRESH_S + int((SD-SD_S)*(PE_DEL_THRESH_L-PE_DEL_THRESH_S)/(SD_L-SD_S))
    if PE_DEL_THRESH > PE_DEL_THRESH_S:
        PE_DEL_THRESH = PE_DEL_THRESH_S
    elif PE_DEL_THRESH < PE_DEL_THRESH_L:
        PE_DEL_THRESH = PE_DEL_LOWEST

    return PE_DEL_THRESH

def calculateLocCovg(NH_REGIONS_FILE,chr_n, bpFirst, bpSecond, PILEUP_THRESH, fBAM, chrHash, 
                    GOOD_REG_THRESH, outerBPs):
    global covHash
    bin_size = 100
    if chr_n not in covHash:
        logging.debug("Calculating coverage for %s", chr_n)
        counterBase, refLoop, cov_100bp, totalCov = 0,0,0,0
        covList = []
        for pileupcolumn in fBAM.pileup(chr_n):
            cov_100bp += pileupcolumn.n
            totalCov += pileupcolumn.n
            counterBase += 1
            refLoop += 1
            if refLoop == bin_size:
                covList.append(1.0*cov_100bp/refLoop)
                cov_100bp, refLoop = 0,0
            if counterBase > CALC_THRESH:
                break
        if len(covList) > 0:
            avgCov = 1.0*totalCov/counterBase
            covHash[chr_n] = avgCov #covList[len(covList)/2]
            #change to debug when test done
            logging.debug("Median coverage of Chr %s written as %f; average was %f",
                          chr_n, covHash[chr_n], avgCov)
        else:
            logging.debug("Unable to calculate coverage in chromosome %s", chr_n)
            print >> stderr, ("Note: unable to calculate coverage in chromosome %s", chr_n)
            covHash[chr_n] = 0

    if bpSecond - bpFirst < 1.25*MIN_PILEUP_THRESH:
        bpFirstL = outerBPs[0]
        bpSecondL = outerBPs[1]
    else:
        bpFirstL = bpFirst
        bpSecondL = bpSecond

    gap = bpSecondL - bpFirstL
    start = .25*gap + bpFirstL
    stop = min(start+.5*gap,start +3*PILEUP_THRESH)
    covLoc, counter, confRegion = 0,0,0
    if stop > start:
        for pileupcolumn in fBAM.pileup(chr_n, start, stop):
            if NH_REGIONS_FILE is None or \
            (chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]):
                covLoc = covLoc + pileupcolumn.n
                counter+=1
                if counter > PILEUP_THRESH:
                    break
        if (counter > MIN_PILEUP_THRESH and (counter > GOOD_REG_THRESH*(stop-start) or counter > PILEUP_THRESH)):
            confRegion = 1
            covLoc = (1.0*covLoc)/(1.0*counter)

    if covHash[chr_n] != 0:
        return 1.0*covLoc/covHash[chr_n], confRegion
    else:
        return 0, 0

def writeVariants(lineAV_split, swap, bnd, support, GT, fAVN, PE_DEL_THRESH, 
                  SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH):
   
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
        elif lineAV_split[11].find("SR") == -1 and lineAV_split[11].find("PE") != -1 and \
            int(lineAV_split[7]) - int(lineAV_split[3]) < PE_DEL_THRESH:
            return
        elif lineAV_split[11].find("SR") != -1 and lineAV_split[11].find("PE") != -1 and \
            int(lineAV_split[7]) - int(lineAV_split[3]) < MIX_DEL_THRESH:
            return

    # correct SPLIT_INS fields just in case
    if  svtype.startswith("TD") or  svtype.startswith("DEL") or  svtype.startswith("BND"):
        lineAV_split_T[8:11] = ["-1", "-1", "-1"]
    lineAV_split_T.extend([str(swap),str(bnd),str(support),GT])
    fAVN.write("%s\n" %("\t".join(lineAV_split_T)))

def covPUFilter(workDir, avFile, vmFile, ufFile, statFile, bamFile,
                NH_REGIONS_FILE, DEL_THRESH, DUP_THRESH, splitINS, 
                PILEUP_THRESH, GOOD_REG_THRESH, minVariantSize):
    
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
    logging.info("DEL_H, DEL_L, DEL, DUP_H, DUP_L, DUP thresholds: %s, %s, %s, %s, %s; %s", DEL_THRESH_H, DEL_THRESH_L, DEL_THRESH, DUP_THRESH_H, DUP_THRESH_L, DUP_THRESH)
    # calculate min PE size based on insert length standard deviation under empirical model
    # aligner tends to mark concordants as discordants when SD is small unless distr v good, seemingly sharp effects
    # around sig_IL of 20 and lower for generally well-behaved distributions. Without rigor, though.
    PE_DEL_THRESH = max(calculateDELThreshPE(SD),UNIV_VAR_THRESH)
    logging.info("PE deletion threshold is %d", PE_DEL_THRESH)

    uniqueFilterSVs = set()
    global chrHash
    global covHash

    for line in fUF:
        uniqueFilterSVs.add(int(line))

    fBAM = pysam.AlignmentFile(bamFile, "rb" )
    if NH_REGIONS_FILE is not None:
        logging.info("Using good regions BED file: %s", NH_REGIONS_FILE)
        formChrHash(NH_REGIONS_FILE, RDL)
    else:
        print >> stderr, "Warning! Not using a good regions file for pile-up filter! This can affect some coverage-based results adversely."
    header = fAV.readline()
    for counter, lineAV in enumerate(fAV):
        counter+=1
        lineAV_split = lineAV.split()
        varNum = int(lineAV_split[0])
        support = 0
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
                (svtype.startswith("INS_half") and lineAV_split[2] == lineAV_split[5]) or\
                svtype == "DN_INS":
                if int(lineAV_split[7])-int(lineAV_split[3]) < UNIV_VAR_THRESH:
                    continue

            ## RUN PILEUP FILTER
            if lineAV_split[11].find("RD") == -1:
                bnd=0
                swap = 0
                GT="."
                if (svtype == "DEL_INS" or svtype == "DEL" or svtype[0:2]== "TD") and \
                    int(lineAV_split[4]) + MIN_PILEUP_THRESH < int(lineAV_split[6]):
                    covLocM, confMiddle = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[2],
                            int(lineAV_split[4]), int(lineAV_split[6]),
                            PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, [int(lineAV_split[3]), int(lineAV_split[7])])

                    if confMiddle == 1:
                        if svtype[0:2] == "TD" and covLocM > DUP_THRESH: 

                            logging.debug("TD confirmed using pileup")
                            svtype = "TD"

                        elif svtype[0:2] == "TD" and covLocM < DUP_THRESH_L:
                            # since bp3 = -1, this will be written as a BND event
                            svtype = "INS_halfRF"

                        elif svtype[0:3] == "DEL" and covLocM < .7: 

                            logging.debug("DEL confirmed using pileup")
                            svtype = "DEL"
                            if covLocM < DEL_THRESH2:
                                GT="GT:1/1"
                            elif covLocM > 3*DEL_THRESH2:
                                GT="GT:0/1"

                        #elif svtype[0:3] == "DEL" and covLocM > DEL_THRESH_L: #1.0:
                            # since bp3 = -1, this will be written as a BND event
                            #svtype = "INS_halfFR"

                elif svtype in ["INS", "INS_I"] or svtype.startswith("INS_C"):

                    del_23, del_12, del_13, dup_23, dup_12, dup_13 = 0, 0, 0, 0, 0, 0
                    swap_12, swap_13 = 0,0
                    # bp2-3
                    # use inner bounds for all
                    start = int(lineAV_split[7])
                    stop = int(lineAV_split[9])
                    covLoc_23, conf_23 = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                            start, stop, PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, [int(lineAV_split[6]), int(lineAV_split[10])])

                    # bp1-2 ("12" will here refer to 1, the paste bp, and the closest of the other 2)
                    start = int(lineAV_split[4])
                    stop = int(lineAV_split[6])
                    outerBPs = [int(lineAV_split[3]), int(lineAV_split[7])]
                    if start > stop:
                        swap_12 = 1
                        start = int(lineAV_split[10])
                        stop = int(lineAV_split[3])
                        outerBPs = [int(lineAV_split[9]), int(lineAV_split[4])]
                    if lineAV_split[2] == lineAV_split[5]:
                        covLoc_12, conf_12 = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                                start, stop, PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, outerBPs)
                    else:
                        convLoc_12, conf_12 = 0,0
                   
                    # bp1-3 ("13" will here refer to 1, the paste bp, and the farther of the other 2)
                    start = int(lineAV_split[4])
                    stop = int(lineAV_split[9])
                    outerBPs = [int(lineAV_split[3]), int(lineAV_split[10])]
                    if start > stop:
                        swap_13 = 1
                        start = int(lineAV_split[7])
                        stop = int(lineAV_split[3])
                        outerBPs = [int(lineAV_split[6]), int(lineAV_split[4])]
                    if lineAV_split[2] == lineAV_split[5]:
                        covLoc_13, conf_13 = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                                start, stop, PILEUP_THRESH, fBAM, chrHash, GOOD_REG_THRESH, outerBPs)
                    else:
                        convLoc_13, conf_13 = 0,0

                    DEL_THRESH_SPLIT = DEL_THRESH

                    if conf_23:
                        if covLoc_23 > DUP_THRESH_L:
                            dup_23 = 1
                        elif covLoc_23 < DEL_THRESH_SPLIT:
                            del_23 = 1
                    if conf_12:
                        if covLoc_12 > DUP_THRESH_L:
                            dup_12 = 1
                        elif covLoc_12 < DEL_THRESH_SPLIT:
                            del_12 = 1
                    if conf_13:
                        if covLoc_13 > DUP_THRESH_L:
                            dup_13 = 1
                        elif covLoc_13 < DEL_THRESH_SPLIT:
                            del_13 = 1

                    lineAV_split1 = list(lineAV_split)
                    if svtype in ["INS", "INS_I"] and lineAV_split[8] != "-1":
                        if splitINS == True and COV > MIN_SPLIT_INS_COV:    
                            logging.debug("Split INS is true: %s", lineAV)
                            # $$club all these if's into 1 if writing as BND
                            if lineAV_split[2] == lineAV_split[5] and (del_12 or del_23 \
                                or (conf_23 and not dup_23)):
                                #1-2 is del
                                logging.debug("Split INS 1")
                                lineAV_split1[1] = "BND"
                                if swap_12:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[9], lineAV_split[10], \
                                                        lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, PE_DEL_THRESH,
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH)
                                #1-3 is TD
                                lineAV_split1[1] = "BND"
                                if swap_13:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[6], lineAV_split[7], \
                                                        lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                else:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[3], lineAV_split[4], \
                                                        lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, PE_DEL_THRESH,
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH)
                                continue

                        # if split above, writing as BND, so don't worry about min size of variants
                        if (conf_12 and covLoc_12 > DUP_THRESH_L) or\
                           (conf_23 and covLoc_23 < DUP_THRESH_L) or\
                            (0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH):
                            bnd = 1

                    elif svtype.startswith("INS_C") and lineAV_split[11].find("PE") != -1 and \
                        lineAV_split[8] != "-1":

                        if splitINS == True and COV > MIN_SPLIT_INS_COV:    
                            logging.debug("Split INS is true -- INS_C: %s", lineAV)
                            if lineAV_split[2] == lineAV_split[5] and (del_12 or del_23):
                                logging.debug("Split INS_C 1")
                                #1-2 is del
                                lineAV_split1[1] = "BND"
                                if swap_12:
                                    lineAV_split1[2:8] = lineAV_split[2], lineAV_split[9], lineAV_split[10], \
                                                        lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, PE_DEL_THRESH,
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH)
                                #2-3 is DEL
                                lineAV_split1[1] = "BND"
                                lineAV_split1[2:8] = lineAV_split[5], lineAV_split[6], lineAV_split[7], \
                                                    lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, PE_DEL_THRESH,
                                              SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH)

                                #1-3 is bnd
                                if int(lineAV_split[12]) >= 3:
                                    lineAV_split1[1] = "BND"
                                    if swap_13:
                                        lineAV_split1[2:8] = lineAV_split[2], lineAV_split[6], lineAV_split[7], \
                                                            lineAV_split[2],  lineAV_split[3],  lineAV_split[4]
                                    else:
                                        lineAV_split1[2:8] = lineAV_split[2], lineAV_split[3], lineAV_split[4], \
                                                            lineAV_split[8],  lineAV_split[9],  lineAV_split[10]
                                    writeVariants(lineAV_split1, swap, bnd, support, GT, fAVN, PE_DEL_THRESH,
                                                  SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH)
                                continue
                        
                        NPE_CLUSTERS_SUPP = int(lineAV_split[12])
                        if (svtype == "INS_C" or svtype == "INS_C_I"):
                            if dup_23 and not dup_12:
                                svtype+="_P"
                            elif dup_12 and not dup_23:
                                svtype+="_P"
                                swap = 1
                        elif (svtype == "INS_C_P" or svtype == "INS_C_I_P") and \
                            ((0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH) or\
                            (conf_12 and covLoc_12 > DUP_THRESH) or \
                            NPE_CLUSTERS_SUPP < 3 or not (DEL_THRESH < covLoc_23 < DUP_THRESH)):
                            bnd = 1

            ## write in BED files
            lineAV_split[1] = svtype
            writeVariants(lineAV_split, swap, bnd, support, GT, fAVN,
                    PE_DEL_THRESH, SR_DEL_THRESH, MIX_DEL_THRESH, UNIV_VAR_THRESH)

    fAV.close()
    fAVN.close()
    fVM.close()
    fUF.close()
    fBAM.close()

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
