import numpy as np
import pysam as ps
from bitarray import bitarray

def readBamStats(statFile):
    rdl, sd, coverage = -1, -1, -1
    with open(statFile, 'r') as fStat:
        for i,line in enumerate(fStat):
            if i == 0:
                rdl = float(line[:-1])
            elif i == 2:
                sd = float(line[:-1])
            elif i == 3:
                coverage = float(line[:-1])
                break
    return rdl, sd, coverage

def findNumberMatches(cigartups):
    nummatches = 0
    if cigartups == None:
        return 0
    for op,oplen in cigartups:
        if op in [0,7]:
            nummatches += int(oplen)

    return nummatches

def readChromosomeLengths(bamfile):
    lengths = {}
    bfile = ps.Samfile(bamfile, 'rb')
    for chrominfo in bfile.header['SQ']:
        lengths[chrominfo['SN']] = chrominfo['LN']
    bfile.close()
    return lengths

def formExcludeHash(chrHash, ignoreBuffer, ignoreBED, lengths):
    """Form double hash table containing chromosomes/genomic units and all corresponding
    locations where alignments are to be excluded from analysis
    Inputs:
        ignoreBED: BED file --'chr start stop'-- containing regions to exclude
        ignoreBuffer: additional buffer margin added to start and stop if desired
    Outputs:
        None
    """
    fo=open(ignoreBED, "r")
    for line in fo:
        line_s = line.split()
        currentTID = line_s[0]
        if currentTID not in chrHash and currentTID in lengths:
            chrHash[currentTID] = bitarray(lengths[currentTID])
            chrHash[currentTID].setall(0)
        if currentTID in chrHash:
            start = int(line_s[1])-ignoreBuffer
            stop = int(line_s[2])+ignoreBuffer
            if stop > start:
                chrHash[currentTID][start:stop] = 1
    fo.close()
    return chrHash

def ignoreRead(chr_l, loc_l, chr_r, loc_r, chrHash):
    """Check if a fragment aligned in particular location is to be excluded from analysis
    Inputs:
        chr_l: chr name of left alignment
        chr_r: chr name of right alignment
        loc_l: left alignment 'arrow-tip' location
        loc_r: right alignment 'arrow-tip' location
    Outputs:
        boolean indicating whether alignment is to be excluded (True) or not (False)
    """
    if chr_l in chrHash and chrHash[chr_l][loc_l] == 1:
        return True
    if chr_r in chrHash and chrHash[chr_r][loc_r] == 1:
        return True
    return False

def countLines(fileName):
    f= open(fileName)
    line_num = -1
    for line_num, _ in enumerate(f):
        pass
    f.close()
    return line_num + 1

