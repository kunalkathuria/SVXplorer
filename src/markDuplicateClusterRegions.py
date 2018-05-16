import sys
import pandas as pd
import numpy as np
import argparse

def writeRemainingRegions(overlap, badRegFile, listT):
    if len(listT) >= 2:
        maxStop = -1
        if overlap == "L": clStart = listT[0].split()[4]
        elif overlap == "R": clStart = listT[0].split()[7]

        for cluster in listT:
            clusterSplit = cluster.split()
            if overlap == "L":
                clChr = clusterSplit[3]
                clStop = int(clusterSplit[5])
            elif overlap == "R":
                clChr = clusterSplit[6]
                clStop = int(clusterSplit[8])
            if clStop > maxStop:
                maxStop = clStop
        badRegFile.write("%s\t%s\t%s\n" %(clChr, clStart, maxStop))

def separateClusters(line, overlap, badRegFile, listT, chr_n, start, chr_prev, prev_stop, chrB, chrB_prev):

    if chr_n != chr_prev or start > prev_stop:
        if len(listT) >= 2:
            maxStop = -1
            if overlap == "L": clStart = listT[0].split()[4]
            elif overlap == "R": clStart = listT[0].split()[7]

            for cluster in listT:
                clusterSplit = cluster.split()
                if overlap == "L":
                    clChr = clusterSplit[3]
                    clStop = int(clusterSplit[5])
                elif overlap == "R":
                    clChr = clusterSplit[6]
                    clStop = int(clusterSplit[8])
                if clStop > maxStop:
                    maxStop = clStop
            badRegFile.write("%s\t%s\t%s\n" %(clChr, clStart, maxStop))
        del listT[:]
        listT.append(line) # add for now if want to examine clusters
    elif chr_n == chr_prev and start <= prev_stop:
        listT.append(line)

    return listT

def markDuplicateClusterRegions(clusterFile, wdir):
    prev_stop01, prev_stop10, chr_prev01, chr_prev10 = 0, 0, "*", "*"
    prev_stop00, prev_stop11, chr_prev00, chr_prev11 = 0, 0, "*", "*"
    chrB_prev01, chrB_prev10, chrB_prev00, chrB_prev11 = "*", "*", "*", "*"
    list01, list10, list11, list00 = [], [], [], []
    fCL = open(clusterFile, "r")
    badRegFile = open(wdir + "/badRegions.bed", "w")

    for line in fCL:
        line_split = line.split()
        chr_n, start, chrB = line_split[3], int(line_split[4]), line_split[6]
        orient = line_split[2]
        if orient == "01":
            list01 = separateClusters(line, "L", badRegFile, list01, chr_n, start, chr_prev01, prev_stop01, chrB, chrB_prev01)
            chr_prev01, prev_stop01, chrB_prev01 = chr_n, int(line_split[5]), chrB
        elif orient == "10":
            list10 = separateClusters(line, "L", badRegFile, list10, chr_n, start, chr_prev10, prev_stop10, chrB, chrB_prev10)
            chr_prev10, prev_stop10, chrB_prev10 = chr_n, int(line_split[5]), chrB
        elif orient == "00":
            list00 = separateClusters(line, "L", badRegFile, list00, chr_n, start, chr_prev00, prev_stop00, chrB, chrB_prev00)
            chr_prev00, prev_stop00, chrB_prev00 = chr_n, int(line_split[5]), chrB
        elif orient == "11":
            list11 = separateClusters(line, "L", badRegFile, list11, chr_n, start, chr_prev11, prev_stop11, chrB, chrB_prev11)
            chr_prev11, prev_stop11, chrB_prev11 = chr_n, int(line_split[5]), chrB
    fCL.close()

    writeRemainingRegions("L", badRegFile, list01)
    writeRemainingRegions("L", badRegFile, list10)
    writeRemainingRegions("L", badRegFile, list00)
    writeRemainingRegions("L", badRegFile, list11)

    #sort on right position and repeat
    prev_stop01, prev_stop10, chr_prev01, chr_prev10 = 0, 0, "*", "*"
    prev_stop00, prev_stop11, chr_prev00, chr_prev11 = 0, 0, "*", "*"
    chrB_prev01, chrB_prev10, chrB_prev00, chrB_prev11 = "*", "*", "*", "*"
    list01, list10, list00, list11 = [], [], [], []

    data = pd.read_table(clusterFile,
                         names=['index', 'ns', 'orient', 'lchr', 'lpos', 'lend',
                                 'rchr', 'rpos', 'rend', 'small'],
                         dtype={'lchr':np.str, 'rchr':np.str, 'orient':np.str})
    data = data.sort_values(by = ['rchr', 'rpos'])
    data.to_csv(wdir + "/allClusters.rs.cleanup.txt", header=None, index=None, sep='\t')
    fCNR = open(wdir + "/allClusters.rs.cleanup.txt", "r")

    for line in fCNR:
        line_split = line.split()
        chr_n, start, chrB = line_split[6], int(line_split[7]), line_split[3]
        orient = line_split[2]
        if orient == "01":
            list01 = separateClusters(line, "R", badRegFile, list01, chr_n, start, chr_prev01, prev_stop01, chrB, chrB_prev01)
            chr_prev01, prev_stop01, chrB_prev01 = chr_n, int(line_split[8]), chrB
        elif orient == "10":
            list10 = separateClusters(line, "R", badRegFile, list10, chr_n, start, chr_prev10, prev_stop10, chrB, chrB_prev10)
            chr_prev10, prev_stop10, chrB_prev10 = chr_n, int(line_split[8]), chrB
        elif orient == "00":
            list00 = separateClusters(line, "R", badRegFile, list00, chr_n, start, chr_prev00, prev_stop00, chrB, chrB_prev00)
            chr_prev00, prev_stop00, chrB_prev00 = chr_n, int(line_split[8]), chrB
        elif orient == "11":
            list11 = separateClusters(line, "R", badRegFile, list11, chr_n, start, chr_prev11, prev_stop11, chrB, chrB_prev11)
            chr_prev11, prev_stop11, chrB_prev11 = chr_n, int(line_split[8]), chrB

    writeRemainingRegions("R", badRegFile, list01)
    writeRemainingRegions("R", badRegFile, list10)
    writeRemainingRegions("R", badRegFile, list00)
    writeRemainingRegions("R", badRegFile, list11)
    fCNR.close()
    badRegFile.close()

if __name__=="__main__":

    PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description='Mark poor-mappability regions based on overlapping clusters')
    PARSER.add_argument('clusterFile', help='file containing all discordant clusters')
    PARSER.add_argument('wdir', help='working directory')
    ARGS = PARSER.parse_args()
    markDuplicateClusterRegions(ARGS.clusterFile, ARGS.wdir)

