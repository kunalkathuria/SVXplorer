#!/usr/bin/env python

# Identify all discordant reads, analyze, and write them to allDiscordants.us.txt

from itertools import izip
import argparse as ap
import numpy as np
import pysam as ps
import logging
import sys
from shared import formExcludeHash, ignoreRead, readChromosomeLengths

## we would not recommended changing any of these
# primary alignment score threshold (0 recommended due to split reads etc.)
AS_THRESH = 0
# %ile of insert-length dist. where IL value indicates discordant alignment
# due to being too large
DISC_PERC = .9985
# %ile of insert-length dist. where IL value indicates discordant alignment
DISC_PERC_NEG = .0001
#min almt score required to use concordant almt in calculating BAM statistics
AS_CALC_THRESH = .999
BIG_NUM = 100000
# %ile of (concordant) IL dist. where the IL value is considered too large
PENALTY_PERC = .99995
# %tile of IL dist. to start distance penalty assessment
DIST_END_PERC = .99

def calcMeanSig(bamfile1, workDir, calc_thresh):
    """Calculate mean & std of insert-length dist. (and other stats listed below)
    Inputs:
        bamfile1: position-sorted BAM file containing all alignments
    Outputs:
        meanIL: mean insert length
        stdIL: stdev in IL
        disc_thresh: IL value of the input distribution at DISC_PERC, less meanIL
        disc_thresh_neg: IL value of the input distribution at DISC_PERC_NEG,
                         less meanIL (neg value)
    Notes:
        Quantities not returned explicitly but written to bam_stats.txt file:
            meanQL: mean query length
            maxIL: maximum insert length
            dist_penalty: IL value of the input distribution at PENALTY_PERC
            dist_end: IL value of the input distribution at DISC_PERC
        'Alignment' here always refers to PE alignment unless otherwise specified
    """
    bamfile = ps.Samfile(bamfile1, "rb")
    summedIL = 0
    summedQL = 0
    counterRead = 0
    stdIL = 0
    maxIL = 0
    counterLoop = 0
    IL_list_us = []

    while counterRead < calc_thresh and counterLoop < 2*calc_thresh:
        try:
            m1 = bamfile.next()
        except StopIteration:
            break

        # AS_CALC_THRESH condition ensures split read will not be included in
        # analysis because "next" may not be mate then
        if m1.is_proper_pair and \
           m1.is_secondary == False and \
           m1.is_supplementary == False and \
           m1.get_tag("AS") > AS_CALC_THRESH*m1.infer_query_length() and m1.template_length > 0:

            IL = m1.template_length
            IL_list_us.append(IL)

            if IL > maxIL:
                maxIL = IL
            summedIL += IL
            summedQL += m1.infer_query_length()
            counterRead += 1
        counterLoop += 1

    try:
        meanIL = summedIL/counterRead
        meanQL = summedQL/counterRead
    except ZeroDivisionError:
        sys.stderr.write("Please check order of name-sorted and position-sorted files supplied.")
        exit(1)

    IL_list = sorted(IL_list_us)
    try:
        dist_end = IL_list[int(DIST_END_PERC*len(IL_list)) - 1]
        dist_penalty = IL_list[int(PENALTY_PERC*len(IL_list)) - 1]
        disc_thresh_neg = IL_list[int(DISC_PERC_NEG*len(IL_list)) - 1] - meanIL
    except:
        sys.stderr.write("Please check value of DISC_PERC (discordancy percentile for IL).")
        exit(1)

    disc_thresh = IL_list[int(DISC_PERC*len(IL_list)) - 1] - meanIL
    bamfile.close()

    binSize = 10
    binDist = binSize
    distHash = {}
    minIL = IL_list[0]
    for item in IL_list:
        if item - minIL < binSize:
            if binDist not in distHash:
                distHash[binDist] = 0
            distHash[binDist] += 1
        else:
            binDist += binSize
            minIL = item
    binDistHash = {}
    for item in distHash:
        for item2 in distHash:
            temp = item2 - item
            if temp >= 0:
                if temp not in binDistHash:
                    binDistHash[temp] = 0
                if temp > 0:
                    binDistHash[temp] += distHash[item]*distHash[item2]
                else:
                    binDistHash[temp] += (distHash[item])*(distHash[item] -1)*.5
    fh = open(workDir + "/binDist.txt", "w")
    for item in binDistHash:
        fh.write("%s\t%s\n" %(item, binDistHash[item]))

    bamfile = ps.Samfile(bamfile1, "rb")
    counterLoop = 0
    counterRead = 0
    while counterRead < calc_thresh and counterLoop < 2*calc_thresh:
        try:
            m1 = bamfile.next()
        except StopIteration:
            break

        if m1.is_proper_pair and not m1.is_secondary and \
            not m1.is_supplementary and \
            m1.get_tag("AS") > AS_CALC_THRESH*m1.infer_query_length() and \
            m1.template_length > 0:
            stdIL += (m1.template_length - meanIL)**2
            counterRead += 1
        counterLoop += 1

    bamfile.close()
    stdIL = stdIL/counterRead
    stdIL = stdIL**(.5)

    bamfile = ps.AlignmentFile(bamfile1, "rb")
    cov = 0
    width = 20
    counterBase = 0
    tidCounter = 0
    tidThresh = 22
    subsampleRate = 10.0
    while True:
        try:
            chrom = bamfile.get_reference_name(tidCounter)
            if tidCounter > tidThresh or chrom == None or chrom == -1:
                break
        except:
            break
        tidCounter += 1
        counterLoop = 0
        for pileupcolumn in bamfile.pileup(chrom):
            cov += pileupcolumn.n
            counterBase += 1
            counterLoop += 1
            if counterLoop > calc_thresh/subsampleRate:
                break
    bamfile.close()
    if counterBase > 0:
        cov = cov/counterBase

    if disc_thresh < 0:
        disc_thresh = 3*stdIL

    logging.debug("meanQL, meanIL, stdIL, cov, maxIL, disc_thresh: %f %f %f %f %f %f", meanQL, meanIL, stdIL, cov, maxIL, disc_thresh)
    with open(workDir + "/bamStats.txt", "w") as fp:
        print >> fp, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s" %(meanQL, meanIL, stdIL, cov, maxIL, disc_thresh, dist_penalty, dist_end)

    return meanQL, meanIL, disc_thresh, disc_thresh_neg

class alignedFragment(object):
    def __init__(self):
        # whether lBound is start or end of left-aligned read (in reference)
        # depends on orientation of read in reference
        self.lBound = -1
        self.rBound = -1
        self.cType = None #orientation, e.g. "01" represents "FR" etc.
        self.cMapType = None #0 if both reads map, 1 if one read maps only
        self.lTID = None
        self.rTID = None
        self.mapQual = -1
        self.discSmall = -1

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.lTID, self.lBound, self.rTID,
                                         self.rBound, self.cType, self.discSmall,                                         self.mapQual)

def findTotalNMatches(al):
    """ Calculate how many bases in alignment match reference, given cigar string
    Inputs:
        al: pysam alignment read from bamfile
    Outputs:
        Number of matches, ratio of number of matches to query length
    """
    MDtagPos = str(al).find("MD",10) + 6
    temp = str(al)[MDtag_pos:].find(")")
    MDtag = str(al)[MDtag_pos:MDtag_pos +temp-1]
    # MD tag is set for perfect matches, so if not set, dubious
    if MDtagPos == 5:
        return 0,0
    integers = "0123456789"
    n_matches = 0
    prev = "0"
    n_subs = 0
    for string in MDtag:
        if string in integers:
            prev += string
        else:
            n_matches += int(prev)
            prev = "0"
            n_subs += 1
    n_matches += int(prev)
    if n_matches > 0:
        return n_matches, float(n_matches)/float(n_matches + n_subs)
    else:
        return 0,0

def formDiscordant(aln1s, aln2s, disc_thresh, disc_thresh_neg, mean_IL, chrHash,
                   nMatchPct_thresh, nMatch_relative_thresh, as_relative_thresh,
                   map_thresh, permutation_thresh, ignoreBED, ignoreTIDList, ignoreTIDAll):
    """ Analyze all discordant alignment pairs, filter them and write to file those that pass in alignedFragment() format
    Inputs:
        aln1s: list of left alignments with same query name
        aln2s: list of right alignments with same query name
        matches of given alignment to those of primary alignment of same query
        disc_thresh: IL value of the input distribution at DISC_PERC, less mean insert length
        disc_thresh_neg: IL value of the input distribution at DISC_PERC, less
        mean insert length (neg value)
        mean_IL: mean insert length
        locations where alignments are to be excluded from analysis
    Outputs:
        dList1: list of filtered PE alignments where both reads map
        dList2: list of filtered PE alignments where only 1 read maps
    """
    dList1 = []
    dList2 = []
    counterLoop = 0
    outerIL = mean_IL
    firstPass = 1

    for al1 in aln1s:
        for al2 in aln2s:
            counterLoop += 1
            al1_match = 0
            al2_match = 0
            l_orient = 2
            r_orient = 2

            # if both almts unmapped, return nothing
            if al1.is_unmapped and al2.is_unmapped:
                return dList1, dList2
            # if primary alignment and below mapping quality threshold, return nothing
            if (al1.mapping_quality < map_thresh or al2.mapping_quality < map_thresh) and counterLoop == 1:
                return dList1, dList2
            # if too many mappings of one fragment, return nothing
            if counterLoop > permutation_thresh:
                return dList1, dList2

            if al1.is_unmapped == False:
                al1_reference_name = al1.reference_name
            else:
                al1_reference_name = None

            if al2.is_unmapped == False:
                al2_reference_name = al2.reference_name
            else:
                al2_reference_name = None

            al1_reference_start, al1_reference_end, al1_is_unmapped, al1_mapping_quality, al1_is_reverse, al1_score, al1_infer_query_length =  al1.reference_start, al1.reference_end, al1.is_unmapped, al1.mapping_quality, al1.is_reverse, float(al1.get_tag("AS")), al1.infer_query_length()

            al2_reference_start, al2_reference_end, al2_is_unmapped, al2_mapping_quality, al2_is_reverse, al2_score, al2_infer_query_length =  al2.reference_start, al2.reference_end, al2.is_unmapped, al2.mapping_quality, al2.is_reverse, float(al2.get_tag("AS")), al2.infer_query_length()

            if al1_reference_start == None and al1_reference_end == None and \
               al2_reference_start == None and al2_reference_end == None:
                return dList1, dList2

            if al1_reference_name in ignoreTIDList or al2_reference_name in ignoreTIDList:
                logging.info("Ignoring almt combination %s and %s", al1_reference_name, al2_reference_name)
                continue

            skip = 0
            for chrI in ignoreTIDAll:
                if al1_reference_name.startswith(chrI) or al2_reference_name.startswith(chrI):
                    skip = 1
                    break
            if skip:
                logging.info("Ignoring almt combination %s and %s as occurs as * entry in ignoreTIDs", al1_reference_name, al2_reference_name)
                continue

            if ignoreBED is not None and ignoreRead(al1_reference_name, al1_reference_start, al2_reference_name, al2_reference_start, chrHash):
                if counterLoop == 1:
                    return dList1, dList2
                continue

            if nMatchPct_thresh != 0 or nMatch_relative_thresh != 0:
                al1_nMatches, al1_nMatchRatio = findTotalNMatches(al1)
                al2_nMatches, al2_nMatchRatio = findTotalNMatches(al2)
            else:
                al1_nMatches = 0
                al1_nMatchRatio = 0
                al2_nMatches = 0
                al2_nMatchRatio = 0

            if firstPass == 1:
                al1_nMatchRatio_prim = al1_nMatchRatio
                al2_nMatchRatio_prim = al2_nMatchRatio
                al1_score_prim = al1_score
                al2_score_prim = al2_score
                al1_nMatchRatio = BIG_NUM
                al2_nMatchRatio = BIG_NUM
                al1_score = BIG_NUM
                al2_score = BIG_NUM
                firstPass = 0

            al1ToPrim_nMatchRatio = 0
            al2ToPrim_nMatchRatio = 0
            if al1_nMatchRatio_prim > 0:
                al1ToPrim_nMatchRatio = float(al1_nMatchRatio)/float(al1_nMatchRatio_prim)
            if al2_nMatchRatio_prim > 0:
                al2ToPrim_nMatchRatio = float(al2_nMatchRatio)/float(al2_nMatchRatio_prim)

            try:
                if al1_score_prim > AS_THRESH and al1_nMatchRatio >= nMatchPct_thresh and \
                    al1ToPrim_nMatchRatio >= nMatch_relative_thresh and al1_score >= as_relative_thresh*al1_score_prim:
                    al1_match = 1
            except:
                pass
            try:
                if al2_score_prim > AS_THRESH and al2_nMatchRatio >= nMatchPct_thresh and \
                    al2ToPrim_nMatchRatio >= nMatch_relative_thresh and al2_score >= as_relative_thresh*al2_score_prim:
                    al2_match = 1
            except:
                pass
            if (al1_match == False and al2_match == False):
                continue

            if al1_match and al2_match:
                newAlmt = alignedFragment()
                if counterLoop == 1:
                    newAlmt.mapQual = min(al1_mapping_quality, al2_mapping_quality)
                newAlmt.cMapType = 0

                # all right if TID different or orientation different. Won't make any difference in those situations.
                if (al1_reference_start <= al2_reference_start):
                    outerIL = abs(al1_reference_end - al2_reference_start) + al1_infer_query_length + al2_infer_query_length
                    newAlmt.lTID = al1_reference_name
                    newAlmt.rTID = al2_reference_name
                    l_orient = int(al1_is_reverse)
                    r_orient = int(al2_is_reverse)

                else:
                    outerIL = abs(al2_reference_end - al1_reference_start) + al1_infer_query_length + al2_infer_query_length
                    newAlmt.lTID = al2_reference_name
                    newAlmt.rTID = al1_reference_name
                    l_orient = int(al2_is_reverse)
                    r_orient = int(al1_is_reverse)

                if (al1_reference_start <= al2_reference_start) and (not al1_is_reverse) and (not al2_is_reverse):
                    newAlmt.lBound = al1_reference_end
                    newAlmt.rBound = al2_reference_end
                elif (al1_reference_start <= al2_reference_start) and (not al1_is_reverse) and (al2_is_reverse):
                    newAlmt.lBound = al1_reference_end
                    newAlmt.rBound = al2_reference_start
                elif (al1_reference_start <= al2_reference_start) and (al1_is_reverse) and (not al2_is_reverse):
                    newAlmt.lBound = al1_reference_start
                    newAlmt.rBound = al2_reference_end
                elif (al1_reference_start <= al2_reference_start) and (al1_is_reverse) and (al2_is_reverse):
                    newAlmt.lBound = al1_reference_start
                    newAlmt.rBound = al2_reference_start
                elif (al1_reference_start > al2_reference_start) and (not al1_is_reverse) and (not al2_is_reverse):
                    newAlmt.lBound = al2_reference_end
                    newAlmt.rBound = al1_reference_end
                elif (al1_reference_start > al2_reference_start) and (not al1_is_reverse) and (al2_is_reverse):
                    newAlmt.lBound = al2_reference_start
                    newAlmt.rBound = al1_reference_end
                elif (al1_reference_start > al2_reference_start) and (al1_is_reverse) and (not al2_is_reverse):
                    newAlmt.lBound = al2_reference_end
                    newAlmt.rBound = al1_reference_start
                elif (al1_reference_start > al2_reference_start) and (al1_is_reverse) and (al2_is_reverse):
                    newAlmt.lBound = al2_reference_start
                    newAlmt.rBound = al1_reference_start

                newAlmt.cType = str(l_orient) + str(r_orient)
                if newAlmt.lTID == newAlmt.rTID and newAlmt.cType == "01" and outerIL - mean_IL < disc_thresh_neg:
                    newAlmt.discSmall = 1

            elif al2_match and not al1_match:

                newAlmt = alignedFragment()
                if counterLoop == 1:
                    newAlmt.mapQual = al2_mapping_quality
                newAlmt.cMapType = 1
                # left tid set by convention for one mate mappings
                newAlmt.lTID = al2_reference_name

                if al2_is_reverse:
                    newAlmt.lBound = al2_reference_start
                    l_orient = 1
                    outerIL = -1
                else:
                    newAlmt.lBound = al2_reference_end
                    l_orient = 0
                    outerIL = -1
                newAlmt.cType = str(l_orient) + str(r_orient)

            elif al1_match and not al2_match:
                newAlmt = alignedFragment()
                if counterLoop == 1:
                    newAlmt.mapQual = al1_mapping_quality
                newAlmt.cMapType = 1
                # left tid set by convention for one mate mappings
                newAlmt.lTID = al1_reference_name

                if al1_is_reverse:
                    newAlmt.lBound = al1_reference_start
                    l_orient = 1
                    outerIL = -1
                else:
                    newAlmt.lBound = al1_reference_end
                    l_orient = 0
                    outerIL = -1
                newAlmt.cType = str(l_orient) + str(r_orient)

            else:
                continue

            # check discordancy even though discordant flag set by aligner
            if newAlmt.cMapType == 1 or \
                    (newAlmt.cMapType == 0 and ( abs(outerIL - mean_IL) > disc_thresh or \
                    newAlmt.cType != "01" or newAlmt.lTID != newAlmt.rTID )):

                    if newAlmt.lBound == -1 and newAlmt.rBound == -1:
                        logging.warning("WARNING STORING BAD DATA ", newAlmt)

                    if newAlmt.cMapType == 0:
                        dList1.append(newAlmt)
                    elif newAlmt.cMapType == 1:
                        dList2.append(newAlmt)
            # if even 1 alignment is concordant, don't use any mappings from read: safe
            else:
                dList1 =[]
                dList2 = []
                return dList1, dList2

    return dList1, dList2

def readNextReadAlignments(bamname):
    """Generator of alignments returning 1 list of alignments with same name at a time
    Inputs:
        bamname: BAM file containing all left (or right) PE alignments
        (e.g. can include secondary almts) almts
    Outputs:
        1 list of alignments with same name at a time
    """
    bamfile = ps.AlignmentFile(bamname)
    alignments = []
    qname = None

    for alignment in bamfile:
        # ignore QC failed, duplicate and supplementary alignments
        if alignment.is_qcfail or \
           alignment.is_duplicate or \
           alignment.is_supplementary: 
           continue

        if qname == None:
            qname = alignment.qname
            alignments.append(alignment)
        elif qname == alignment.qname:
            alignments.append(alignment)
        else:
            yield qname, alignments
            alignments = [alignment]
            qname = alignment.qname

    if qname != None:
        yield qname, alignments

def writeDiscordantFragments(workDir, readAlmts1, readAlmts2, bamfile, debug,
                             ignoreBED, ignoreChr, permutation_thresh,
                             calc_thresh, nMatchPct_thresh,
                             nMatch_relative_thresh, as_relative_thresh,
                             map_thresh):
    ignoreTIDs = set()
    ignoreTIDAll = set()
    chrHash = {}

    # calculate some basic stats
    rdl, mean_IL, disc_thresh, disc_thresh_neg = calcMeanSig(bamfile, workDir, calc_thresh)

    # read the lengths of the chromosomes
    chromosome_lengths = readChromosomeLengths(bamfile)

    if ignoreChr is not None:
        with open(ignoreChr, 'r') as f:
            for line in f:
                chrI = line.strip().split()[0]
                if not chrI.startswith("*"):
                    ignoreTIDs.add(chrI)
                    logging.info("Chromosome %s will be ignored.", chrI)
                else:
                    ignoreTIDAll.add(chrI[1:])
                    logging.info("Chr names starting with %s will be ignored", chrI[1:])

    ignoreBuffer = 0*rdl
    if ignoreBED is not None:
        logging.info("Regions in %s will be ignored", ignoreBED)
        chrHash = formExcludeHash(chrHash, ignoreBuffer, ignoreBED, chromosome_lengths)

    # read discordant alignments and write all possible discordant pairs to file.
    with open("%s/allDiscordants.us.txt" % workDir, "w") as almtFile:
        currentFrag = 1
        logging.info('Started reading discordant pairs')
        for (q1, aln1s),(q2, aln2s) in izip(readNextReadAlignments(readAlmts1), readNextReadAlignments(readAlmts2)):

            if (currentFrag % 100000) == 0:
                logging.debug("%d fragments analyzed", currentFrag)
            try:
                assert q1[:-2] == q2[:-2]
            except AssertionError:
                logging.debug("Assertion error on query name being same in writeDiscordants")
                sys.stderr.write("Please check if reads pass vendor checks \
                \n{0}\n{1}\nQuitting.\n" .format(q1,q2))
                exit(1)

            dList1, dList2 = formDiscordant(aln1s, aln2s, disc_thresh,
                    disc_thresh_neg, mean_IL, chrHash, nMatchPct_thresh,
                    nMatch_relative_thresh, as_relative_thresh, map_thresh,
                    permutation_thresh, ignoreBED, ignoreTIDs, ignoreTIDAll)
            for item in dList1:
                print >> almtFile, "%s\t%s" %(currentFrag, item)
            for item in dList2:
                print >> almtFile, "%s\t%s" %(currentFrag, item)
            currentFrag += 1
        logging.info('Finished reading discordant pairs')

if __name__ == "__main__":

    # parse arguments
    PARSER = ap.ArgumentParser(description="""
    Filter, format and write all PE discordants from 2 input discordant BAM
    files, each containing first- or second-in-pair discordant alignments, to
    allDiscordants.us.txt""", formatter_class=ap.ArgumentDefaultsHelpFormatter)
    PARSER.add_argument('workDir', help='Work directory')
    PARSER.add_argument('readAlmts1', help='BAM file containing first-in-pair discodrant PE almts')
    PARSER.add_argument('readAlmts2', help='BAM file containing second-in-pair discordant PE almts')
    PARSER.add_argument('bamfile', help='Position-sorted BAM alignment file')
    PARSER.add_argument('-d', action='store_true', dest='debug',
        help='print debug information')
    PARSER.add_argument('-i', default=None, dest='ignoreBED',
        help='Exclude-regions file in BED format')
    PARSER.add_argument('-c', default=None, dest='ignoreChr',
        help='File listing chromosomes to exclude, one per line')
    PARSER.add_argument('-p', default=20, dest='permutation_thresh', type=int,
        help='If all alignments with same query name exceed this number, do not use any in analysis')
    PARSER.add_argument('-t', default=1000000, dest='calc_thresh', type=int,
        help='Max number of concordant alignments to use in calculating BAM statistics')
    PARSER.add_argument('-n', default=0, dest='nMatchPct_thresh', type=int,
        help='Threshold of matching base pairs in alignment to reference')
    PARSER.add_argument('-r', default=0, dest='nMatch_relative_thresh',
        type=int,
        help='Threshold of ratio of matching base pairs in \
        given alignment to primary alignment of same query')
    PARSER.add_argument('-a', default=2, dest='as_relative_thresh', type=int,
        help='Threshold of relative alignment score of given alignment to primary alignment of same name\
        --set to less than 1 (.95 recommended) if using secondary almts')
    PARSER.add_argument('-m', default=10, dest='map_thresh', type=int,
        help='Mapping quality threshold')
    ARGS = PARSER.parse_args()

    LEVEL = logging.INFO
    if ARGS.debug:
        LEVEL = logging.DEBUG

    logging.basicConfig(level=LEVEL,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    writeDiscordantFragments(ARGS.workDir, ARGS.readAlmts1, ARGS.readAlmts2,
                             ARGS.bamfile, ARGS.debug, ARGS.ignoreBED,
                             ARGS.ignoreChr, ARGS.permutation_thresh,
                             ARGS.calc_thresh, ARGS.nMatchPct_thresh,
                             ARGS.nMatch_relative_thresh,
                             ARGS.as_relative_thresh,
                             ARGS.map_thresh)

    logging.shutdown()
