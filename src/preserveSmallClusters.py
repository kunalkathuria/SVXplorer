import pysam
import argparse
import logging
import numpy as np
from shared import readChromosomeLengths

#gloabl
SVHashPE = {}

def formExcludeHashVN(clusterFile, lengths):
    logging.info('Started reading PE clusters for small cluster preservation')
    global SVHashPE
    fCl=open(clusterFile, "r")
    for line in fCl:
        line_s = line.split()
        varNum = int(line_s[0])
        currentTID1 = line_s[3]
        currentTID2 = line_s[6]
        if currentTID1 not in SVHashPE:
            SVHashPE[currentTID1] = np.zeros(lengths[currentTID1])
        if currentTID2 not in SVHashPE:
            SVHashPE[currentTID2] = np.zeros(lengths[currentTID2])
        SVHashPE[currentTID1][int(line_s[4]):int(line_s[5])] = varNum
        SVHashPE[currentTID2][int(line_s[7]):int(line_s[8])] = varNum
    logging.info('Finished forming PE cluster hash table')
    fCl.close()

def preserveSmallClusters(bamfileSR, clusterFile, mapThresh, preserveSize, slop, wdir):
  
    chrLengths = readChromosomeLengths(bamfileSR)
    global SVHashPE
    bamfileSRH = pysam.Samfile(bamfileSR,"rb")
    iObjects = []
    formExcludeHashVN(clusterFile, chrLengths)
    sizeAddition = {}
    counter = 0
    while True:
        try:
            sr1 = bamfileSRH.next()
            sr2 = bamfileSRH.next()
        except StopIteration:
            break
        if counter % 500 == 0:
            logging.debug("Done with %s SRs", counter)
        counter+= 1

        if sr1.qname == sr2.qname:
            sr_bp1 = sr1.reference_start
            sr_bp2 = sr2.reference_start
            sr_bp1_tid = sr1.reference_name
            sr_bp2_tid = sr2.reference_name

            # ignore marked chromosomes
            if sr1.mapping_quality < mapThresh or sr2.mapping_quality < mapThresh:
                continue

            match = 0
            bp, bp2 = sr_bp1, sr_bp2
            tid, tid_2 = sr_bp1_tid, sr_bp2_tid

            if (tid in SVHashPE) and (tid_2 in SVHashPE) and \
                SVHashPE[tid][bp] != 0 and SVHashPE[tid][bp] == SVHashPE[tid_2][bp2]:
                match = 1
            
            #store changed size as value of hash
            if match:
                varNum = SVHashPE[tid][bp] 
                if varNum not in sizeAddition:
                    sizeAddition[varNum] = 1
                else:    
                    sizeAddition[varNum]+= 1

    bamfileSRH.close()
    fCl = open(clusterFile,"r")
    fClw = open(clusterFile + ".p", "w")

    for line in fCl:
        line_split = line.split()
        clSize, varNum = int(line_split[1]), int(line_split[0])
        if varNum in sizeAddition:
            clSize+= sizeAddition[varNum]
        if clSize >= preserveSize:
            #write to file
            fClw.write("%s" %line)

if __name__=="__main__":

    PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description='Preserve low-support PE clusters if SRs support them as well')
    PARSER.add_argument('clusterFile', help='file containing all discordant clusters')
    PARSER.add_argument('splitBAM', help='file containing name-sorted split reads')
    PARSER.add_argument('mapThresh', help='mapping quality threshold for split reads')
    PARSER.add_argument('preserveSize', help='minimum preserve size for PE clusters after SR support addition')
    PARSER.add_argument('wdir', help='working directory')
    ARGS = PARSER.parse_args()
    logging.basicConfig(filename='%s/run.log' % ARGS.wdir,
                                level=logging.INFO,
                                format='%(asctime)s %(levelname)s %(message)s',
                                datefmt='%m/%d/%Y %I:%M:%S %p',
                                filemode='w')
    preserveSmallClusters(ARGS.splitBAM, ARGS.clusterFile, int(ARGS.mapThresh), int(ARGS.preserveSize), 0, ARGS.wdir)
