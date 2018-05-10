import sys
import pandas as pd
import numpy as np
OV_MARGIN = int(sys.argv[3])

def writeRemainingCls(overlap, badRegFile, listT, fp):
    if len(listT) == 1:
        #print overlap, "Writing", listT[0]
        fp.write("%s" %listT[0])
    elif len(listT) >= 2:
        for entry in listT:
            #print "l > 1:", overlap, "Writing", entry
            fp.write("%s" %entry)
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

def separateClusters(overlap, badRegFile, listT, fp, chr_n, start, chr_prev, prev_stop, prevOrient, orient):

    condn1 = False
    if chr_n == chr_prev:
        condn1 = True
        #print "Condn 2", prevOrient, orient
        # FR close but not overlapping with non-FR; RR or FF overlapping with FR or RF
        if overlap == "L" and prevOrient != "01" and orient == "01" \
            and start >= prev_stop and abs(start - prev_stop) < OV_MARGIN:
            #print "Condn 2 a"
            listT.append(line)
        elif overlap == "R" and prevOrient == "01" and orient != "01" \
            and start >= prev_stop and abs(start - prev_stop) < OV_MARGIN:
            #print "Condn 2 b"
            listT.append(line)
        elif prevOrient in ["00","11"] and orient not in ["00","11"] and \
            (abs(start - prev_stop) < OV_MARGIN or start < prev_stop):
            #print "Condn 2 c"
            listT.append(line)
        elif prevOrient not in ["00","11"] and orient in ["00","11"] and \
            (abs(start - prev_stop) < OV_MARGIN or start < prev_stop):
            listT.append(line)
            #print "Condn 2 d"
        else:
            condn1 = False
    if (not condn1) or chr_n != chr_prev:
        #print "Condn 1"
        if len(listT) == 1:
            #print "S", overlap, "Writing", listT[0]
            fp.write("%s" %listT[0])
        elif len(listT) >= 2:
            for entry in listT:
                #print "S; l > 1:", overlap, "Writing", entry
                fp.write("%s" %entry)
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
    return listT

if __name__=="__main__":

    wdir = sys.argv[2]
    prev_stop, chr_prev, prevOrient = 0, "*", "22"
    listR = []
    fCL = open(sys.argv[1], "r")
    fCN = open(wdir + "/allClusters.clean.str.ls.txt", "w")
    badRegFile = open(wdir + "/badRegions.bed", "w")

    for line in fCL:
        line_split = line.split()
        chr_n, start = line_split[3], int(line_split[4])
        orient = line_split[2]
        listR = separateClusters("L", badRegFile, listR, fCN, chr_n, start, chr_prev, prev_stop, prevOrient, orient)
        chr_prev, prev_stop, prevOrient = chr_n, int(line_split[5]), orient

    fCL.close()
    writeRemainingCls("L", badRegFile, listR, fCN)
    fCN.close()
    
    listR = []
    prev_stop, chr_prev, prevOrient = 0, "*", "22"
    #sort on right position and repeat
    data = pd.read_table("%s/allClusters.clean.str.ls.txt" % wdir,
                         names=['index', 'ns', 'orient', 'lchr', 'lpos', 'lend',
                                'rchr', 'rpos', 'rend', 'small'],
                         dtype={'lchr':np.str, 'rchr':np.str, 'orient':np.str})
    data = data.sort_values(by = ['rchr', 'rpos'])
    data.to_csv(wdir + "/allClusters.str.rs.txt", header=None, index=None, sep='\t')
    fCNR = open(wdir + "/allClusters.str.rs.txt", "r")
    fCNW = open(wdir + "/allClusters.clean.str.txt", "w")
    for line in fCNR:
        #print line, listR
        line_split = line.split()
        chr_n, start = line_split[6], int(line_split[7])
        orient = line_split[2]
        listR = separateClusters("R", badRegFile, listR, fCNW, chr_n, start, chr_prev, prev_stop, prevOrient, orient)
        chr_prev, prev_stop, prevOrient = chr_n, int(line_split[8]), orient
        #print "List at end:", listR, "\n"
    writeRemainingCls("R", badRegFile, listR, fCNW)
    fCNW.close()
    fCNR.close()
    badRegFile.close()
