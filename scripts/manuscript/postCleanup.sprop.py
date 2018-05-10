import sys
import pandas as pd
import numpy as np
import pysam as ps

SUPP_PERC=float(sys.argv[6]) #.65
MIN_SUPP=int(sys.argv[7]) #10

def determineMin1(chr1, start1, chrHash):
    bit = 1
    index = start1
    while index >= 0:
        index-=1
        if chrHash[chr1][index] == 0:
            break
    return index

def readChromosomeLengths(bamfile):
    lengths = {}
    bfile = ps.Samfile(bamfile, 'rb')
    for chrominfo in bfile.header['SQ']:
        lengths[chrominfo['SN']] = chrominfo['LN']
    bfile.close()
    return lengths

def formExcludeHash(chrHash, ignoreBuffer, ignoreBED, lengths):
    fo=open(ignoreBED, "r")
    for line in fo:
        line_s = line.split()
        currentTID = line_s[0]
        if currentTID not in chrHash:
            chrHash[currentTID] = np.zeros(lengths[currentTID])
        chrHash[currentTID][int(line_s[1])-ignoreBuffer:int(line_s[2])+ignoreBuffer] = 1
    return chrHash

if __name__=="__main__":

    wdir = sys.argv[3]
    ignoreBed = sys.argv[2]
    ignBuffer = int(sys.argv[5])
    clusterFile = sys.argv[1]
    fCN = open(wdir + "/allClusters.postClean.txt", "w")
    fSuppHist = open(wdir + "/suppHist.txt", "w")
    fCl = open(clusterFile, "r")
    chrHash = {}
    chrLengths = readChromosomeLengths(sys.argv[4])
    chrHash = formExcludeHash(chrHash, ignBuffer, ignoreBed, chrLengths)
    compHash = {}

    for line in fCl:
        line_split = line.split()
        chrL = line_split[3]
        chrR = line_split[6]
        startL, stopL = int(line_split[4]), int(line_split[5])
        startR, stopR = int(line_split[7]), int(line_split[8])
        #if startL == 325201:
            #print chrL, chrL in chrHash, chrR, chrR in chrHash
        condn1 = False
        if (chrL not in chrHash) or (chrL in chrHash and chrHash[chrL][startL] == 0 and \
            chrHash[chrL][stopL] == 0): 
            condn1 = True
            if (chrR not in chrHash) or (chrR in chrHash and chrHash[chrR][startR] == 0 and \
                chrHash[chrR][stopR] == 0): 
                fCN.write("%s" %line)
            else:
                condn1 = False
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
                fCN.write("%s" %compHash[compID][maxSuppIndex])
            fSuppHist.write("%s\n" %(1.0*suppMax/suppSum))

    fCl.close()
    fCN.close()

