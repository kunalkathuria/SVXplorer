#!/usr/bin/env python

# Filter variants based on local coverage in variant regions and write final SVs in BED format

from sys import stderr
import pysam
from collections import Counter
from bitarray import bitarray
import argparse
import logging

# global constants
DEL_THRESH2 = .125
DEL_THRESH_L = .8
DUP_THRESH_L = 1.2
MIN_PILEUP_THRESH = 80
CALC_THRESH = 1000000

DEL_CONF_THRESH = .7
UNIV_VAR_THRESH = 100
INS_VAR_THRESH = 100
INS_VAR_THRESH_PE = 50
SR_DEL_THRESH = 80
MIX_DEL_THRESH = 50
PE_DEL_THRESH_S = 250
PE_DEL_THRESH_L = 150
SD_S = 14
SD_L = 24

def formChrHash(NH_REGIONS_FILE):
    print "Forming pile-up hash table..."
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

def readBamStats(statFile):
    rdl, sd, coverage = -1, -1, -1
    with open(statFile, 'r') as fStat:
        for i,line in enumerate(fStat):
            if i == 0:
                rdl=float(line[:-1])
            elif i==2:
                sd=float(line[:-1])
                break
    return rdl, sd

def calculateDELThreshPE(SD):
    PE_DEL_THRESH=PE_DEL_THRESH_S + int((SD-SD_S)*(PE_DEL_THRESH_L-PE_DEL_THRESH_S)/(SD_L-SD_S))
    if PE_DEL_THRESH > PE_DEL_THRESH_S:
        PE_DEL_THRESH = PE_DEL_THRESH_S
    elif PE_DEL_THRESH < PE_DEL_THRESH_L:
        PE_DEL_THRESH = UNIV_VAR_THRESH

    return PE_DEL_THRESH

def calculateLocCovg(NH_REGIONS_FILE,chr_n, bpFirst, bpSecond, PILEUP_THRESH, fBAM, chrHash, covHash, 
                    GOOD_REG_THRESH):
    bin_size = 100
    if chr_n not in covHash:
        logging.info("Calculating coverage for %s", chr_n)
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
            covHash[chr_n] = covList[len(covList)/2] 
            avgCov = 1.0*totalCov/counterBase
            #change to debug when test done
            logging.info("Median coverage of Chr %s written as %f; average was %f",
                          chr_n, covHash[chr_n], avgCov)
        else:
            print >> stderr, ("Warning! No good bases in chromosome %s", chr_n)
            covHash[chr_n] = 0

    gap = bpSecond - bpFirst
    start = .25*gap + bpFirst
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

def writeVariants(lineAV_split, swap, bnd, support, GT, fAVN, PE_DEL_THRESH):
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

    lineAV_split.extend([str(swap),str(bnd),str(support),GT])
    fAVN.write("%s\n" %("\t".join(lineAV_split)))

def covPUFilter(workDir, avFile, vmFile, ufFile, statFile, bamFile,
                NH_REGIONS_FILE, DEL_THRESH, DUP_THRESH, splitINS, 
                PILEUP_THRESH, GOOD_REG_THRESH):
    fAV = open(avFile,"r")
    fVM=open(vmFile,"r")
    fUF = open(ufFile,"r")
    fAVN = open(workDir+"/allVariants.pu.txt","w")
    logging.info("Writing final bedpe files using coverage information")
    RDL,SD = readBamStats(statFile)
    logging.info("Some stats from BAM. RDL: %d, Sd: %f", 
                  RDL, SD)

    # calculate min PE size based on insert length standard deviation under empirical model
    # aligner tends to mark concordants as discordants when SD is small unless distr v good, seemingly sharp effects
    # around sig_IL of 20 and lower for Poisson and other related distributions.
    PE_DEL_THRESH = calculateDELThreshPE(SD)
    logging.info("PE deletion threshold is %d", PE_DEL_THRESH)

    covHash = {}
    uniqueFilterSVs = set()
    chrHash = {}
    for line in fUF:
        uniqueFilterSVs.add(int(line))

    fBAM = pysam.AlignmentFile(bamFile, "rb" )
    if NH_REGIONS_FILE is not None:
        formChrHash(NH_REGIONS_FILE)
    else:
        print >> stderr, "Warning! Not using a good regions file for pile-up filter! This can affect some coverage-based results adversely."
    for counter, lineAV in enumerate(fAV):
        support = 0
        for entry in fVM:
            support = len(entry.split()) - 1
            break
        counter+=1
        lineAV_split = lineAV.split()
        varNum = int(lineAV_split[0])

        if varNum in uniqueFilterSVs:
            svtype = lineAV_split[1]

            ## skip SV if size thresholds not met
            if svtype.startswith("TD") or svtype.startswith("INV") or svtype.startswith("DEL"):
                if int(lineAV_split[7])-int(lineAV_split[3]) < UNIV_VAR_THRESH:
                    continue
            elif svtype.startswith("INS") and svtype not in ["INS_C","INS_C_I"]:
                if lineAV_split[11].find("PE") != -1 and \
                0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH_PE:
                    continue
                elif lineAV_split[11].find("PE") == -1 and \
                0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH:
                    continue
            elif lineAV_split == "INS_C" or lineAV_split == "INS_C_I":
                if (0 < int(lineAV_split[10])-int(lineAV_split[6]) < INS_VAR_THRESH and \
                0 < int(lineAV_split[7])-int(lineAV_split[3]) < INS_VAR_THRESH) or \
                (0 < int(lineAV_split[7])-int(lineAV_split[9]) < INS_VAR_THRESH and \
                0 < int(lineAV_split[4])-int(lineAV_split[6]) < INS_VAR_THRESH):
                    continue

            ## RUN PILEUP FILTER
            if lineAV_split[11].find("RD") == -1:
                bnd=0
                swap = 0
                GT=""
                if (svtype == "DEL_INS" or svtype == "DEL" or svtype[0:2]== "TD") and int(lineAV_split[4]) + MIN_PILEUP_THRESH < int(lineAV_split[6]):
                    covLocM, confMiddle = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[2],
                            int(lineAV_split[4]), int(lineAV_split[6]),
                            PILEUP_THRESH, fBAM, chrHash, covHash, GOOD_REG_THRESH)

                    if confMiddle == 1:
                        if svtype[0:2] == "TD" and covLocM < DEL_THRESH: 

                            logging.info("TD confirmed using pileup")
                            svtype = "TD"

                        elif svtype[0:2] == "TD" and covLocM < 1.0:
                            svtype = "BND"

                        elif svtype[0:3] == "DEL" and covLocM < DEL_THRESH: 

                            logging.info("DEL confirmed using pileup")
                            svtype = "DEL"
                            if covLocM < DEL_THRESH2:
                                GT="GT:1/1"
                            elif covLocM > 3*DEL_THRESH2:
                                GT="GT:0/1"

                        elif svtype[0:3] == "DEL" and covLocM > DUP_THRESH_L: #1.0:
                            # since bp3 = -1, this will be written as a BND event
                            svtype = "INS"

                elif len(svtype) > 2 and (svtype == "INS" \
                    or svtype == "INS_I") and int(lineAV_split[7]) + MIN_PILEUP_THRESH < \
                    int(lineAV_split[9]) and lineAV_split[8] != "-1":

                    #bp2-3
                    covLoc, confReg = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                        int(lineAV_split[7]), int(lineAV_split[9]),
                        PILEUP_THRESH, fBAM, chrHash, covHash, GOOD_REG_THRESH)
                    if confReg and covLoc < DUP_THRESH_L:
                        bnd = 1

                elif len(svtype) > 4 and lineAV_split[11].find("PE") != -1 and \
                    svtype[0:5] == "INS_C" and lineAV_split[8] != "-1":
                    del_23 = 0
                    del_12 = 0
                    dup_23 = 0
                    dup_12 = 0
                    #bp2-3
                    start = int(lineAV_split[7])
                    stop = int(lineAV_split[9])
                    if start > stop:
                        start, stop = stop, start
                    covLoc_23, conf_23 = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                            start, stop, PILEUP_THRESH, fBAM, chrHash, covHash, GOOD_REG_THRESH)
                    #bp1-2
                    start = int(lineAV_split[4])
                    stop = int(lineAV_split[6])
                    if start > stop:
                        start, stop = stop, start
                    if lineAV_split[2] == lineAV_split[5]:
                        covLoc_12, conf_12 = calculateLocCovg(NH_REGIONS_FILE,lineAV_split[5],
                                start, stop, PILEUP_THRESH, fBAM, chrHash, covHash, GOOD_REG_THRESH)
                    else:
                        convLoc_12, conf_12 = 0,0

                    if conf_23:
                        if covLoc_23 > DUP_THRESH_L:
                            dup_23 = 1
                        elif covLoc_23 < DEL_THRESH:
                            del_23 =1
                    if conf_12:
                        if covLoc_12 > DUP_THRESH_L:
                            dup_12 = 1
                        elif covLoc_12 < DEL_THRESH:
                            del_12 =1

                    confINSBP = 0
                    if (svtype == "INS_C" or svtype == "INS_C_I"):
                        if dup_23 and not dup_12:
                            svtype+="_P"
                            confINSBP = 1
                        elif dup_12 and not dup_23:
                            svtype+="_P"
                            confINSBP = 1
                            swap = 1

            ## write in BED files
            lineAV_split[1] = svtype
            writeVariants(lineAV_split, swap, bnd, support, GT, fAVN,
                    PE_DEL_THRESH)

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
                ARGS.PILEUP_THRESH, ARGS.GOOD_REG_THRESH)

    logging.shutdown()
