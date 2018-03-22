import numpy as np
import pysam as ps

def readChromosomeLengths(bamfile):
    lengths = {}
    bfile = ps.Samfile(bamfile, 'rb')
    for chrominfo in bfile.header['SQ']:
        lengths[chrominfo['SN']] = chrominfo['LN']
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
    prevTID = "*" # invalid value
    fo=open(ignoreBED, "r")
    for line in fo:
        line_s = line.split()
        currentTID = line_s[0]
        if currentTID != prevTID:
            chrHash[currentTID] = np.zeros(lengths[currentTID])
        chrHash[currentTID][int(line_s[1])-ignoreBuffer:int(line_s[2])+ignoreBuffer] = 1
        prevTID = currentTID
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


