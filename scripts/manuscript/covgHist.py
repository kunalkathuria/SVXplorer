import pysam
import sys

covHash = {}
MIN_PILEUP_THRESH = 80

def calculateLocCovg(chr_n, bpFirst, bpSecond, fBAM, fCov):
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
            covHash[chr_n] = covList[len(covList)/2] 
            avgCov = 1.0*totalCov/counterBase
            #change to debug when test done
            fCov.write("Chr %s (med, avg):%s\t%s\n", chr_n, covHash[chr_n], avgCov)
        else:
            covHash[chr_n] = 0

    gap = bpSecond - bpFirst
    start = .25*gap + bpFirst
    stop = min(start+.5*gap,start +3*PILEUP_THRESH)
    covLoc, counter, confRegion = 0,0,0
    if stop > start:
        for pileupcolumn in fBAM.pileup(chr_n, start, stop):
            covLoc = covLoc + pileupcolumn.n
            counter+=1
            if counter > PILEUP_THRESH:
                break
        if counter > MIN_PILEUP_THRESH: 
            confRegion = 1
            covLoc = (1.0*covLoc)/(1.0*counter)

    if covHash[chr_n] != 0:
        return 1.0*covLoc/covHash[chr_n], confRegion
    else:
        return 0, 0

if __name__=="__main__":

    bamFile = sys.argv[1]
    bedFile = sys.argv[2]
    fBAM=pysam.AlignmentFile(bamFile, "rb")
    f=open(bedFile, "r")
    fw=open(bedFile + ".hist","w")
    fCov=open("chrCovg.txt","w")
    for line in f:
        line_split = line.split()
        chr_n = line_split[0]
        bpFirst = int(line_split[1])
        bpSecond = int(line_split[2])

        cov,_ = calculateLocCovg(chr_n, bpFirst, bpSecond, fBAM, fCov)
        fw.write("%s\t%s\n" %(cov,bpSecond - bpFirst))

    f.close()
    fw.close()
    fCov.close()
