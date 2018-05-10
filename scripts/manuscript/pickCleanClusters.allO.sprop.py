import sys
import pandas as pd
import numpy as np

def writeRemainingCls(overlap, badRegFile, listT, fCN):
    if len(listT) == 1:
        fCN.write("%s" %listT[0])
    elif len(listT) >= 2:
        for entry in listT:
            fCN.write("%s" %entry)
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

def separateClusters(overlap, badRegFile, listT, fCN, chr_n, start, chr_prev, prev_stop, chrB, chrB_prev):

    if chr_n != chr_prev or start > prev_stop: # start > prev_start + BUFFER
        if len(listT) == 1:
            fCN.write("%s" %listT[0])
        elif len(listT) >= 2:
            for entry in listT:
                fCN.write("%s" %entry)
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

if __name__=="__main__":

    wdir = sys.argv[2]
    BUFFER = 1000 # maxClMarginLength
    prev_stop01, prev_stop10, chr_prev01, chr_prev10 = 0, 0, "*", "*"
    prev_stop00, prev_stop11, chr_prev00, chr_prev11 = 0, 0, "*", "*"
    chrB_prev01, chrB_prev10, chrB_prev00, chrB_prev11 = "*", "*", "*", "*"
    list01, list10, list11, list00 = [], [], [], []
    fCL = open(sys.argv[1], "r")
    fCN = open(wdir + "/allClusters.clean.ls.txt", "w")
    badRegFile = open(wdir + "/badRegions.bed", "w")

    for line in fCL:
        line_split = line.split()
        chr_n, start, chrB = line_split[3], int(line_split[4]), line_split[6]
        orient = line_split[2]
        if orient == "01":
            list01 = separateClusters("L", badRegFile, list01, fCN, chr_n, start, chr_prev01, prev_stop01, chrB, chrB_prev01)
            chr_prev01, prev_stop01, chrB_prev01 = chr_n, int(line_split[5]), chrB
        elif orient == "10":
            list10 = separateClusters("L", badRegFile, list10, fCN, chr_n, start, chr_prev10, prev_stop10, chrB, chrB_prev10)
            chr_prev10, prev_stop10, chrB_prev10 = chr_n, int(line_split[5]), chrB
        elif orient == "00":
            list00 = separateClusters("L", badRegFile, list00, fCN, chr_n, start, chr_prev00, prev_stop00, chrB, chrB_prev00)
            chr_prev00, prev_stop00, chrB_prev00 = chr_n, int(line_split[5]), chrB
        elif orient == "11":
            list11 = separateClusters("L", badRegFile, list11, fCN, chr_n, start, chr_prev11, prev_stop11, chrB, chrB_prev11)
            chr_prev11, prev_stop11, chrB_prev11 = chr_n, int(line_split[5]), chrB
    fCL.close()

    writeRemainingCls("L", badRegFile, list01, fCN)
    writeRemainingCls("L", badRegFile, list10, fCN)
    writeRemainingCls("L", badRegFile, list00, fCN)
    writeRemainingCls("L", badRegFile, list11, fCN)
    fCN.close()

    #sort on right position and repeat
    prev_stop01, prev_stop10, chr_prev01, chr_prev10 = 0, 0, "*", "*"
    prev_stop00, prev_stop11, chr_prev00, chr_prev11 = 0, 0, "*", "*"
    chrB_prev01, chrB_prev10, chrB_prev00, chrB_prev11 = "*", "*", "*", "*"
    list01, list10, list00, list11 = [], [], [], []
    data = pd.read_table("%s/allClusters.clean.ls.txt" % wdir,
                         names=['index', 'ns', 'orient', 'lchr', 'lpos', 'lend',
                                'rchr', 'rpos', 'rend', 'small'],
                         dtype={'lchr':np.str, 'rchr':np.str, 'orient':np.str})
    data = data.sort_values(by = ['rchr', 'rpos'])
    data.to_csv(wdir + "/allClusters.rs.txt", header=None, index=None, sep='\t')
    fCNR = open(wdir + "/allClusters.rs.txt", "r")
    fCNW = open(wdir + "/allClusters.clean.allO.txt", "w")
    for line in fCNR:
        line_split = line.split()
        chr_n, start, chrB = line_split[6], int(line_split[7]), line_split[3]
        orient = line_split[2]
        if orient == "01":
            list01 = separateClusters("R", badRegFile, list01, fCNW, chr_n, start, chr_prev01, prev_stop01, chrB, chrB_prev01)
            chr_prev01, prev_stop01, chrB_prev01 = chr_n, int(line_split[8]), chrB
        elif orient == "10":
            list10 = separateClusters("R", badRegFile, list10, fCNW, chr_n, start, chr_prev10, prev_stop10, chrB, chrB_prev10)
            chr_prev10, prev_stop10, chrB_prev10 = chr_n, int(line_split[8]), chrB
        elif orient == "00":
            list00 = separateClusters("R", badRegFile, list00, fCNW, chr_n, start, chr_prev00, prev_stop00, chrB, chrB_prev00)
            chr_prev00, prev_stop00, chrB_prev00 = chr_n, int(line_split[8]), chrB
        elif orient == "11":
            list11 = separateClusters("R", badRegFile, list11, fCNW, chr_n, start, chr_prev11, prev_stop11, chrB, chrB_prev11)
            chr_prev11, prev_stop11, chrB_prev11 = chr_n, int(line_split[8]), chrB

    writeRemainingCls("R", badRegFile, list01, fCNW)
    writeRemainingCls("R", badRegFile, list10, fCNW)
    writeRemainingCls("R", badRegFile, list00, fCNW)
    writeRemainingCls("R", badRegFile, list11, fCNW)
    fCNW.close()
    fCNR.close()
    badRegFile.close()
