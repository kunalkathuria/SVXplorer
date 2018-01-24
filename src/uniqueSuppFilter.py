#!/usr/bin/env python
# Kunal Kathuria 1/18
# Pick variants having minimum unique support for each support category (SR, PE, MIXED)

import sys
import math
from collections import Counter
import argparse

def readBamStats(statFile):
    f=open(statFile,"r")
    for i,line in enumerate(f):
        if i==3:
            break
        elif i==2:
            sigIL = float(line.split()[0])
    return float(line.split()[0]), sigIL

def formMQSet(mapThresh, mqSet, allDiscordantsFile):
    f=open(allDiscordantsFile,"r")
    throw = f.readline() #absorb header
    for line in f:
        line_split = line.split()
        frag = int(line_split[0])
        mq = int(line_split[7])
        # mapping qual threshold for uniquely supporting fragments
        if frag not in mqSet and mq >= mapThresh:
            mqSet.add(frag)
    f.close()

def calculateSVThresh(SVType, SVSupp):
    if SVType == "INV" or SVType.find("INS") != -1:
        disjThresh = complex_thresh
    elif SVSupp.find("PE") == -1 and SVSupp.find("SR") != -1:
        disjThresh = sr_thresh
    elif SVSupp.find("PE") != -1 and SVSupp.find("SR") == -1:
        disjThresh = pe_thresh
    elif SVSupp.find("PE") != -1 and SVSupp.find("SR") != -1:
        disjThresh = mix_thresh

    return disjThresh

def uniquenessFilter(fragmentList, nInputVariants, mqSet, allDiscordantsFile, mapThresh,\
    variantMapFile, allVariantFile, rdFragIndex, rdVarIndex):

    nSVs = 0
    disjointness = [0]*nInputVariants
    pickV = [0]*nInputVariants
    varNum = []
    fVM=open(variantMapFile,"r")
    fAV=open(allVariantFile,"r")
    fUS=open(workDir+"/variants.uniqueSupport.txt","w")
    fUF=open(workDir+"/variants.uniqueFilter.txt","w")
    formMQSet(mapThresh, mqSet, allDiscordantsFile)

    print "Applying uniqueness filter."
    nFragOccrns = Counter(fragmentList)
    # obtain disjointness count
    # assumed any given fragment only appears once in set list
    for counter, line in enumerate(fVM):
        currentSet = map(int, line.split())
        varNum.append(currentSet[0])
        currentSet = currentSet[1:] 
        disjThresh = -1

        for line in fAV:
            SVType = line.split()[1]
            SVSupp = line.split()[11]
            disjThresh = calculateSVThresh(SVType, SVSupp)
            break

        for elem in currentSet:
            #pick those supported by RD automatically       
            if elem >= rdFragIndex:
                if nFragOccrns[elem] == 1:
                    disjointness[counter]+=1
                pickV[counter] = 2
            #if SR, no secondaries so is above MQ_THRESH automatically
            elif nFragOccrns[elem] == 1 and (elem in mqSet or elem < 0):
                disjointness[counter]+=1
            if disjointness[counter] == disjThresh:
                break
    fAV.seek(0)
    for g,item in enumerate(disjointness):
        # not stored in memory
        for line in fAV:
            SVType = line.split()[1]
            SVSupp = line.split()[11]
            disjThresh = calculateSVThresh(SVType, SVSupp)
            break

        if (item >= disjThresh) or (item >= 1 and pickV[g] == 2):
           fUF.write("%s\n" %varNum[g]) 
           nSVs+=1
        #make list of variants with unique support till threshold was reached
        if varNum[g] < rdVarIndex:
            fUS.write("%s %s\n" %(varNum[g], disjointness[g]))

    fVM.close()
    fAV.close()
    fUS.close()
    fUF.close()
    return nSVs

def readVariantMap(filename, allFrags):
    f=open(filename, 'r')
    for counter, line in enumerate(f):
        parsed = map(int, line.split())
        for frag in parsed[1:]:
            allFrags.append(frag)
    f.close()
    return counter + 1
    
if __name__ == "__main__":  

    # parse arguments
    #try:
    parser = argparse.ArgumentParser(description='Apply uniquenes Filter: pick variants \
        having minimum unique support for each support category (SR, PE, MIXED)')
    parser.add_argument('workDir', help='Work directory')
    parser.add_argument('statFile', help='File containing BAM statistics, typically\
        bamStats.txt')
    parser.add_argument('variantMapFile', help='File containing variant map, typically variantMap._.txt')
    parser.add_argument('allVariantFile', help='File containing list of variants, typically allVariants._.txt')
    parser.add_argument('allDiscordantsFile', help='File containing all discordants with MQ,\
        typically allDiscordants.txt')
    parser.add_argument('-a', default=6, dest='pe_thresh_max', type=int,\
        help='Maximum allowed support threshold for PE-only variants (dynamic)') 
    parser.add_argument('-b', default=6, dest='sr_thresh_max', type=int,\
        help='Maximum allowed support threshold for SR-only variants')
    parser.add_argument('-c', default=3, dest='pe_thresh_min', type=int,\
        help='Minimum allowed support threshold for PE-only variants')
    parser.add_argument('-d', default=3, dest='sr_thresh_min', type=int,\
        help='Minimum allowed support threshold for SR-only variants')
    parser.add_argument('-e', default=4, dest='mix_thresh', type=int,\
        help='Mixed variants support threshold')
    parser.add_argument('-f', default=4, dest='complex_thresh', type=int,\
        help='Complex variants (INV, INS) support threshold')
    parser.add_argument('-g', default=10, dest='mapThresh', type=int,\
        help='Mapping quality threshold for fragments uniquely supporting variant')
    parser.add_argument('-k', default=4, dest='single_thresh', type=int,\
        help='If using one universal support threshold (will need to look at code to use)')
    parser.add_argument('-l', default=3000000, dest='rdVarIndex', type=int,\
        help=argparse.SUPPRESS)
    parser.add_argument('-i', default=100000000, dest='rdFragIndex', type=int,\
        help=argparse.SUPPRESS)
    parser.add_argument('-j', default=-1, dest='ILBoost', type=int,\
        help=argparse.SUPPRESS)
    args = parser.parse_args()
    #except:
        #parser.print_help()
        #exit(1)
 
    allFrags = []
    mqSet = set()
    single_thresh = args.single_thresh
    workDir = args.workDir
    statFile = args.statFile
    variantMapFile = args.variantMapFile
    allVariantFile = args.allVariantFile
    allDiscordantsFile = args.allDiscordantsFile
    mapThresh = args.mapThresh
    rdVarIndex = args.rdVarIndex
    rdFragIndex = args.rdFragIndex
    #constants
    pe_thresh_max = args.pe_thresh_max
    sr_thresh_max = args.sr_thresh_max
    pe_thresh_min = args.pe_thresh_min
    sr_thresh_min = args.sr_thresh_min
    mix_thresh = args.mix_thresh
    # for complex variants, various criteria are satisfied, so threshold static/relaxed
    complex_thresh = args.complex_thresh
    # linear model to calculate support threshold by category
    # familiar developers may tweak model here directly
    pe_low = 3
    pe_high = 6
    covg_low = 5
    covg_high = 50
    sr_low = 4
    sr_high = 6
    il_low1 = 25
    il_low2 = 35
    ILBoost = 0

    # apply above-mentioned support threshold model
    [covg,sig_il] = readBamStats(statFile)
    if il_low1 <= int(sig_il) <= il_low2:
        ILBoost = args.ILBoost
    sr_thresh = math.floor(sr_low + (covg-covg_low)*1.0*(sr_high - sr_low)/(covg_high - covg_low))
    pe_thresh = round(ILBoost + pe_low + (covg-covg_low)*1.0*(pe_high - pe_low)/(covg_high - covg_low))
    if pe_thresh > pe_thresh_max:
        pe_thresh = pe_thresh_max
    if sr_thresh > sr_thresh_max:
        sr_thresh = sr_thresh_max
    if pe_thresh < pe_thresh_min:
        pe_thresh = pe_thresh_min
    if sr_thresh < sr_thresh_min:
        sr_thresh = sr_thresh_min
    print "sr, pe threshes, covg are:", sr_thresh, pe_thresh, covg

    nInputVariants = readVariantMap(variantMapFile, allFrags)
    nSVs = uniquenessFilter(allFrags, nInputVariants, mqSet, allDiscordantsFile, mapThresh,\
        variantMapFile, allVariantFile, rdFragIndex, rdVarIndex)
    fNSVs= open(workDir+"/NSVs.txt","w")
    fNSVs.write("%s\n" %nSVs)
    fNSVs.close()
