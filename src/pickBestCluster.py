import sys
import pandas as pd
import numpy as np
import pysam as ps
import argparse
from shared import readChromosomeLengths, formExcludeHash

SUPP_PERC=.5
MIN_SUPP=10

def determineMin1(chr1, start1, chrHash):
    bit = 1
    index = start1
    while index >= 0:
        index-=1
        if chrHash[chr1][index] == 0:
            break
    return index

def pickBestCluster(clusterFile, wdir, ignoreRegions, sampleBAM):

    ignBuffer = 0
    fCN = open(wdir + "/allClusters.postClean.txt", "w")
    fSuppHist = open(wdir + "/suppHist.txt", "w")
    fCl = open(clusterFile, "r")
    chrHash = {}
    chrLengths = readChromosomeLengths(sampleBAM)
    chrHash = formExcludeHash(chrHash, ignBuffer, ignoreRegions, chrLengths)
    compHash = {}

    for line in fCl:
        line_split = line.split()
        chrL = line_split[3]
        chrR = line_split[6]
        startL, stopL = int(line_split[4]), int(line_split[5])
        startR, stopR = int(line_split[7]), int(line_split[8])
        condn1 = False
        if (chrL not in chrHash) or (chrL in chrHash and chrHash[chrL][startL] == 0 and \
            (stopL >= len(chrHash[chrL]) or chrHash[chrL][stopL] == 0)):
            if stopL >= len(chrHash[chrL]):
                line_split[5] = str(len(chrHash[chrL]) - 1)
            if (chrR not in chrHash) or (chrR in chrHash and chrHash[chrR][startR] == 0 and \
                (stopR >= len(chrHash[chrR]) or chrHash[chrR][stopR] == 0)): 
                condn1 = True
                if stopR >= len(chrHash[chrR]):
                    line_split[8] = str(len(chrHash[chrR]) - 1)
                fCN.write("%s\n" %"\t".join(line_split))

        if not condn1:
            if chrL in chrHash: 
                minBR = determineMin1(chrL, startL, chrHash)
                compID = (chrL, minBR)
            elif chrR in chrHash:
                minBR = determineMin1(chrR, startR, chrHash)
                compID = (chrR, minBR)

            if compID not in compHash:
                compHash[compID] = []
                compHash[compID].append(line)
            else:
                compHash[compID].append(line)

    for compID in compHash:
        suppMax, maxSuppIndex = -1, -1
        suppSum = 0
        for k,cluster in enumerate(compHash[compID]):
            cluster_split = cluster.split()
            suppC = int(cluster_split[1])
            if suppC > suppMax:
                maxSuppIndex = k
                suppMax = suppC
            suppSum+= suppC
        if maxSuppIndex > -1:
             #print "Sum is", suppSum, .75*suppSum, suppMax, maxSuppIndex, listC[maxSuppIndex]
            if suppMax >= max(MIN_SUPP,SUPP_PERC*suppSum):
                clLine = compHash[compID][maxSuppIndex]
                clLine_split = clLine.split()
                if clLine_split[3] in chrHash and int(clLine_split[5]) >= len(chrHash[clLine_split[3]]):
                    clLine_split[5] = str(len(chrHash[clLine_split[3]]) - 1)
                if clLine_split[6] in chrHash and int(clLine_split[8]) >= len(chrHash[clLine_split[6]]):
                    clLine_split[8] = str(len(chrHash[clLine_split[6]]) - 1)
                fCN.write("%s\n" %"\t".join(clLine_split))
            fSuppHist.write("%s\n" %(1.0*suppMax/suppSum))

    fCl.close()
    fCN.close()


if __name__=="__main__":

    PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description='Pick one cluster from regions marked as poor mappability regions')
    PARSER.add_argument('clusterFile', help='file containing all discordant clusters')
    PARSER.add_argument('wdir', help='working directory')
    PARSER.add_argument('ignoreRegions', help='file containing poor-mappability regions in BED format')
    PARSER.add_argument('sampleBAM', help='BAM file of sample')
    ARGS = PARSER.parse_args()
    pickBestCluster(ARGS.clusterFile, ARGS.wdir, ARGS.ignoreRegions, ARGS.sampleBAM)
