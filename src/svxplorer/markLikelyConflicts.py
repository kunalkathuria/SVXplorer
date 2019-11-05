import sys
import argparse

OV_MARGIN = int(sys.argv[3])  # TODO: this should be argparse variable


def writeRemainingCls(overlap, badRegFile, listT):
    if len(listT) >= 2:
        maxStop = -1
        if overlap == "L":
            clStart = listT[0].split()[4]
        elif overlap == "R":
            clStart = listT[0].split()[7]

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
        badRegFile.write("%s\t%s\t%s\n" % (clChr, clStart, maxStop))


def separateClusters(line, overlap, badRegFile, listT, chr_n, start, chr_prev, prev_stop, prevOrient, orient):
    condn1 = False
    if chr_n == chr_prev:
        condn1 = True
        # FR close but not overlapping with non-FR; RR or FF overlapping with FR or RF
        if overlap == "L" and prevOrient != "01" and orient == "01" \
                and start >= prev_stop and abs(start - prev_stop) < OV_MARGIN:
            listT.append(line)
        elif overlap == "R" and prevOrient == "01" and orient != "01" \
                and start >= prev_stop and abs(start - prev_stop) < OV_MARGIN:
            listT.append(line)
        elif prevOrient in ["00", "11"] and orient not in ["00", "11"] and \
                (abs(start - prev_stop) < OV_MARGIN or start < prev_stop):
            listT.append(line)
        elif prevOrient not in ["00", "11"] and orient in ["00", "11"] and \
                (abs(start - prev_stop) < OV_MARGIN or start < prev_stop):
            listT.append(line)
        else:
            condn1 = False
    if (not condn1) or chr_n != chr_prev:
        if len(listT) >= 2:
            maxStop = -1
            if overlap == "L":
                clStart = listT[0].split()[4]
            elif overlap == "R":
                clStart = listT[0].split()[7]

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
            badRegFile.write("%s\t%s\t%s\n" % (clChr, clStart, maxStop))
        del listT[:]
        listT.append(line)  # add for now if want to examine clusters
    return listT


def markUnlikelyClusterRegions(clusterFile, wdir):
    prev_stop, chr_prev, prevOrient = 0, "*", "22"
    listR = []
    fCL = open(clusterFile, "r")
    badRegFile = open(wdir + "/badRegions.bed", "w+")

    for line in fCL:
        line_split = line.split()
        chr_n, start = line_split[3], int(line_split[4])
        orient = line_split[2]
        listR = separateClusters(line, "L", badRegFile, listR, chr_n, start, chr_prev, prev_stop, prevOrient, orient)
        chr_prev, prev_stop, prevOrient = chr_n, int(line_split[5]), orient

    fCL.close()
    writeRemainingCls("L", badRegFile, listR)

    listR = []
    prev_stop, chr_prev, prevOrient = 0, "*", "22"
    fCNR = open(wdir + "/allClusters.rs.cleanup.txt", "r")
    for line in fCNR:
        line_split = line.split()
        chr_n, start = line_split[6], int(line_split[7])
        orient = line_split[2]
        listR = separateClusters(line, "R", badRegFile, listR, chr_n, start, chr_prev, prev_stop, prevOrient, orient)
        chr_prev, prev_stop, prevOrient = chr_n, int(line_split[8]), orient

    writeRemainingCls("R", badRegFile, listR)
    fCNR.close()
    badRegFile.close()


def main():
    PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Mark poor-mappability regions based on unlikely clusters')
    PARSER.add_argument('clusterFile', help='file containing all discordant clusters')
    PARSER.add_argument('wdir', help='working directory')
    ARGS = PARSER.parse_args()
    markUnlikelyClusterRegions(ARGS.clusterFile, ARGS.wdir)


if __name__ == "__main__":
    main()
