#!/usr/bin/env python
#Kunal Kathuria 12/17
# Identify all discordant reads, analyze, and write them to allDiscordants.us.txt
from sys import argv, stderr
import pysam
from itertools import izip
import argparse
#from time import clock #uncomment if timing loop in main

def calcMeanSig(bamfile1, penalty_perc, disc_perc, disc_perc_neg, calc_thresh, as_calc_thresh, workDir):
    """Calculate mean and standard deviation of insert 
    length distribution (and other stats listed below)
    Inputs:
        bamfile1: position-sorted BAM file containing all concordant and discordant alignments
        penalty_perc: percentile of (concordant) insert length distribution where the IL value is considered too large
        for 2 alignments with this physical separation to belong to the same cluster (see manuscript for details)
        disc_perc: percentile of insert length distribution where the IL value indicates discordant alignment due
        to being too large
        disc_perc_neg: percentile of insert length distribution where the IL value indicates discordant alignment
        due to being too small
        calc_thresh: max number of concordant alignments to use in calculating BAM statistics
        as_calc_thresh: min almt score required to use concordant almt in calculating BAM statistics
    Outputs:
        meanIL: mean insert length 
        stdIL: stdev in IL
        cov: subsampled average coverage
        disc_thresh: IL value of the input distribution at disc_perc, less meanIL
        disc_thresh_neg: IL value of the input distribution at disc_perc_neg, less meanIL (neg value)
    Notes:
        Quantities not returned explicitly but written to bam_stats.txt file:
            meanQL: mean query length
            maxIL: maximum insert length
            dist_penalty: IL value of the input distribution at penalty_perc
            dist_end: IL value of the input distribution at disc_perc
        'Alignment' here always refers to PE alignment unless otherwise specified
    """
    bamfile = pysam.Samfile(bamfile1,"rb")
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

        # as_calc_thresh condition ensures split read will not be included in analysis \
        # because "next" may not be mate then
        if m1.is_proper_pair and \
           m1.is_secondary == False and \
           m1.is_supplementary == False and \
           m1.get_tag("AS") > as_calc_thresh*m1.infer_query_length() and m1.template_length > 0:

            IL = m1.template_length
            IL_list_us.append(IL) 
             
            if IL > maxIL:
                maxIL = IL
            summedIL+= IL
            summedQL+= m1.infer_query_length()
            counterRead+= 1
        counterLoop+=1

    try:
        meanIL = summedIL/counterRead
        meanQL = summedQL/counterRead
    except ZeroDivisionError:
        sys.stderr.write("Please check order of name-sorted and position-sorted files supplied.")
        exit(1)

    IL_list = sorted(IL_list_us)
    try:
        dist_end = IL_list[int(disc_perc*len(IL_list)) - 1]
        dist_penalty = IL_list[int(penalty_perc*len(IL_list)) - 1]
        disc_thresh_neg = IL_list[int(disc_perc_neg*len(IL_list)) - 1] - meanIL
    except:
        sys.stderr.write("Please check value of disc_perc (discordancy percentile for IL).")
        exit(1)

    disc_thresh = dist_end - meanIL
    bamfile.close()
  
    binSize = 10
    binDist = binSize
    distHash = {} 
    minIL = IL_list[0]
    for item in IL_list:
        if item - minIL < binSize:
            if binDist not in distHash:
                distHash[binDist]=0
            distHash[binDist]+=1
        else:
            binDist+=binSize
            minIL = item
    binDistHash = {}
    for item in distHash:
        for item2 in distHash:
            temp = item2 - item
            if temp >= 0: 
                if temp not in binDistHash:
                    binDistHash[temp] = 0
                if temp > 0:
                    binDistHash[temp]+= distHash[item]*distHash[item2]
                else:
                    binDistHash[temp]+= (distHash[item])*(distHash[item] -1)*.5
    fh = open(workDir + "/binDist.txt","w")
    for item in binDistHash:
        fh.write("%s\t%s\n" %(item,binDistHash[item]))         
        
    bamfile = pysam.Samfile(bamfile1,"rb")
    counterLoop=0
    counterRead=0
    while counterRead < calc_thresh and counterLoop < 2*calc_thresh:
        try:
            m1 = bamfile.next()
        except StopIteration:
            break

        if m1.is_proper_pair and not m1.is_secondary and not m1.is_supplementary \
                and m1.get_tag("AS") > as_calc_thresh*m1.infer_query_length() and m1.template_length > 0:
            stdIL+= (m1.template_length - meanIL)**2
            counterRead+=1
        counterLoop+=1

    bamfile.close()

    bamfile = pysam.AlignmentFile(bamfile1,"rb")
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
        tidCounter+=1
        counterLoop = 0
        for pileupcolumn in bamfile.pileup(chrom):
            cov+=pileupcolumn.n
            counterBase+=1
            counterLoop+=1
            if counterLoop > calc_thresh/subsampleRate:
                break
    bamfile.close()
    stdIL = stdIL/counterRead
    stdIL = stdIL**(.5)
    cov=cov/counterBase

    if disc_thresh < 0:
        disc_thresh = 3*stdIL
    if verbose:
        print "meanQL, meanIL, stdIL, cov, maxIL, disc_thresh:", meanQL, meanIL, stdIL, cov, maxIL, disc_thresh
    fp = open(workDir + "/bamStats.txt","w")
    fp.write("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n" %(meanQL, meanIL, stdIL, cov, maxIL, disc_thresh, dist_penalty, dist_end))
    fp.close()
    return meanQL, meanIL, disc_thresh, disc_thresh_neg

def formExcludeHash(ignoreBuffer, chrHash, ignoreBED):
    """Form double hash table containing chromosomes/genomic units and all corresponding
    locations where alignments are to be excluded from analysis
    Inputs:
        ignoreBED: BED file --'chr start stop'-- containing regions to exclude
        chrHash: empty hash table to be populated and used
        ignoreBuffer: additional buffer margin added to start and stop if desired
    Outputs:
        None
    """
    prevTID = "*" # invalid value
    fo=open(ignoreBED,"r")
    for line in fo:
        line_s = line.split()
        currentTID = line_s[0]
        if currentTID != prevTID:
            chrHash[currentTID] = {}
        for x in range(int(line_s[1])-ignoreBuffer, int(line_s[2])+ignoreBuffer):
            chrHash[currentTID][x] = 1
        prevTID = currentTID

def ignoreRead(chr_l, loc_l, chr_r, loc_r):
    """Check if a fragment aligned in particular location is to be excluded from analysis
    Inputs:
        chr_l: chr name of left alignment
        chr_r: chr name of right alignment
        loc_l: left alignment 'arrow-tip' location
        loc_r: right alignment 'arrow-tip' location
    Outputs:
        boolean indicating whether alignment is to be excluded (True) or not (False) 
    """
    if chr_l in chrHash and loc_l in chrHash[chr_l]:
        return 1
    if chr_r in chrHash and loc_r in chrHash[chr_r]:
        return 1
    return 0

def recordExcludeChr(ignoreChr, ignoreTIDs):
    """ Form list of chromosomes/genomic units all of whose alignments are to be
    excluded from analysis
    Inputs:
        ignoreChr: file listing these, one per line
        ignoreTIDs: empty list to be populated by these
    Outputs:
        None
    """
    f= open(ignoreChr, "r")
    for line in f:
        ignoreTIDs.append(line.split()[0])

class alignedFragment(object):
    def __init__(self):
        # whether lBound is start or end of left-aligned read (in reference)
        # depends on orientation of read in reference
        self.lBound = None
        self.rBound = None
        self.cType = None #orientation, e.g. "01" represents "FR" etc.
        self.cMapType = None #0 if both reads map, 1 if one read maps only
        self.lTID = None
        self.rTID = None
        self.mapQual = -1
        self.discSmall = -1

    def __str__(self):
        return "%s %s %s %s %s %s %s" % (self.lTID, self.lBound, self.rTID,
            self.rBound, self.cType, self.discSmall, self.mapQual)

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
            prev+= string
        else:
            n_matches+= int(prev)
            prev = "0"
            n_subs+=1
    n_matches+=int(prev)
    if n_matches > 0:
        return n_matches, float(n_matches)/float(n_matches + n_subs)
    else:
        return 0,0

def formDiscordant(aln1s, aln2s, permutation_thresh, map_thresh, as_thresh, 
        nMatchPct_thresh, nMatch_relative_thresh, as_relative_thresh, disc_thresh, disc_thresh_neg, big_num, mean_IL, chrHash, ignoreTIDs):
    """ Analyze all discordant alignment pairs, filter them and write to file those that pass in alignedFragment() format
    Inputs:
        aln1s: list of left alignments with same query name
        aln2s: list of right alignments with same query name
        permutation_thresh: if all alignment pairs with same query name exceed this number, 
        do not use any in analysis
        map_thresh: mapping quality threshold
        as_thresh: alignment score threshold for primary almt
        nMatchPct_thresh: threshold of number of matching base pairs in alignment to reference
        nMatch_relative_thresh: threshold of ratio of base-pair
        matches of given alignment to those of primary alignment of same query
        as_relative_thresh: threshold of relative alignment score of given
        alignment to primary alignment of same query
        disc_thresh: IL value of the input distribution at disc_perc, less mean insert length
        disc_thresh_neg: IL value of the input distribution at disc_perc, less 
        mean insert length (neg value)
        big_num: a number greater than 1
        mean_IL: mean insert length
        chrHash: double hash table containing chromosomes/genomic units and all corresponding
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
            counterLoop+= 1
            al1_match = 0
            al2_match = 0
            l_orient = 2
            r_orient = 2
   
            al1_reference_name, al1_reference_start, al1_reference_end, al1_is_unmapped, al1_mapping_quality, al1_is_reverse, al1_score, al1_infer_query_length = \
                    al1.reference_name, al1.reference_start, al1.reference_end, al1.is_unmapped, al1.mapping_quality, al1.is_reverse, float(al1.get_tag("AS")), al1.infer_query_length()

            al2_reference_name, al2_reference_start, al2_reference_end, al2_is_unmapped, al2_mapping_quality, al2_is_reverse, al2_score, al2_infer_query_length = \
                    al2.reference_name, al2.reference_start, al2.reference_end, al2.is_unmapped, al2.mapping_quality, al2.is_reverse, float(al2.get_tag("AS")), al2.infer_query_length()

            if al1_reference_start == None and al1_reference_end == None and \
                    al2_reference_start == None and al2_reference_end == None:
                return dList1, dList2

            try:
                if (len(al1_reference_name) > 1 and al1_reference_name[0:2] == "GL") or \
		    (len(al2_reference_name) > 1 and al2_reference_name[0:2] == "GL") \
                    or al1_reference_name in ignoreTIDs or al2_reference_name in ignoreTIDs:
                    continue

                if ignoreBED != "none" and ignoreRead(al1_reference_name, al1_reference_start, al2_reference_name, al2_reference_start):
                    if counterLoop == 1:
                        return dList1, dList2
                    continue

            except:
                pass

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
                al1_nMatchRatio = big_num
                al2_nMatchRatio = big_num
                al1_score = big_num
                al2_score = big_num
                firstPass = 0

            al1ToPrim_nMatchRatio = 0
            al2ToPrim_nMatchRatio = 0
            if al1_nMatchRatio_prim > 0:
                al1ToPrim_nMatchRatio = float(al1_nMatchRatio)/float(al1_nMatchRatio_prim)
            if al2_nMatchRatio_prim > 0:
                al2ToPrim_nMatchRatio = float(al2_nMatchRatio)/float(al2_nMatchRatio_prim)

            try:
                if al1_score_prim > as_thresh*al1_infer_query_length and al1_nMatchRatio >= nMatchPct_thresh and \
                    al1ToPrim_nMatchRatio >= nMatch_relative_thresh and al1_score >= as_relative_thresh*al1_score_prim:
                    al1_match = 1
            except:
                pass
            try:
                if al2_score_prim > as_thresh*al2_infer_query_length and al2_nMatchRatio >= nMatchPct_thresh and \
                    al2ToPrim_nMatchRatio >= nMatch_relative_thresh and al2_score >= as_relative_thresh*al2_score_prim:
                    al2_match = 1
            except:
                pass
            if (al1_match == False and al2_match == False):
                continue

            # if primary alignment and below mapping quality threshold, return nothing
            if (al1_mapping_quality < map_thresh or al2_mapping_quality < map_thresh) and counterLoop == 1:
                return dList1, dList2

            # if too many mappings of one fragment, return nothing
            if counterLoop > permutation_thresh:
                dList1 = []
                dList2 = []
                return dList1, dList2

            if al1_match and al2_match:
                newAlmt = alignedFragment()
                if counterLoop == 1:
                    newAlmt.mapQual = min(al1_mapping_quality, al2_mapping_quality)
                else:
                    newAlmt.mapQual = -1
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
                else:
                    newAlmt.mapQual = -1
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
                    newAlmt.mapQual = al2_mapping_quality
                else:
                    newAlmt.mapQual = -1
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
                        print "WARNING STORING BAD DATA ", newAlmt

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
    bamfile = pysam.AlignmentFile(bamname)
    alignments = []
    qname = None

    for alignment in bamfile:
        if qname == None:
            qname = alignment.qname
            alignments.append(alignment)
        elif qname == alignment.qname:
            alignments.append(alignment)
        else:
            yield qname,alignments
            alignments = [alignment]
            qname = alignment.qname

    if qname != None:
        yield qname,alignments
    
if __name__ == "__main__":

    # parse arguments
    try:
        parser = argparse.ArgumentParser(description='Filter, format and write all PE discordants from 2 \
            input discordant BAM files, each containing first- or second-in-pair discordant alignments, to allDiscordants.us.txt')
        parser.add_argument('workDir', help='Work directory')
        parser.add_argument('readAlmts1', help='BAM file containing first-in-pair discodrant PE almts')
        parser.add_argument('readAlmts2', help='BAM file containing second-in-pair discordant PE almts')
        parser.add_argument('bamfile', help='Position-sorted alignment file (BAM)')
        parser.add_argument('-v', default=0, dest='verbose', type=int,\
            help='1 for verbose output')
        parser.add_argument('-i', default='none', dest='ignoreBED',\
            help='Exclude-regions file in BED format')
        parser.add_argument('-c', default='none', dest='ignoreChr',\
            help='File listing chromosomes to exclude, one per line')
        parser.add_argument('-p', default=20, dest='permutation_thresh', type=int,\
            help='If all alignments with same query name exceed this number, do not use any in analysis')
        parser.add_argument('-t', default=1000000,dest='calc_thresh', type=int,\
            help='Max number of concordant alignments to use in calculating BAM statistics')
        parser.add_argument('-n', default=0, dest='nMatchPct_thresh', type=int,\
            help='Threshold of matching base pairs in alignment to reference')
        parser.add_argument('-r', default=0, dest='nMatch_relative_thresh', type=int,\
            help='Threshold of ratio of matching base pairs in \
            given alignment to primary alignment of same query')
        parser.add_argument('-a', default=2, dest='as_relative_thresh', type=int,\
            help='Threshold of relative alignment score of given alignment to primary alignment of same name\
            --set to less than 1 (.95 recommended) if using secondary almts')
        parser.add_argument('-m', default=10, dest='map_thresh', type=int, \
            help='Mapping quality threshold')
        parser.add_argument('-s', default=0, dest='as_thresh', type=int, \
            help=argparse.SUPPRESS)
        args = parser.parse_args()

    except:
        parser.print_help()
        exit(1)

    ## statistical constants -- most hidden from SVC front-end user. 
    ## see function defs for documentation of specific vars
    ## default settings generally recommended.
    permutation_thresh = args.permutation_thresh 
    calc_thresh = args.calc_thresh
    nMatchPct_thresh = args.nMatchPct_thresh
    nMatch_relative_thresh = args.nMatch_relative_thresh
    as_relative_thresh = args.as_relative_thresh
    map_thresh = args.map_thresh
    # primary alignment score threshold (0 recommended due to split reads etc.)  
    as_thresh = args.as_thresh
    ## we would not recommended changing any of these
    disc_perc = .9985
    disc_perc_neg = .0001
    as_calc_thresh = .999
    big_num = 100000
    penalty_perc = .999999

    readAlmts1 = args.readAlmts1 #workDir + "/aln1s.bam" #read 1 discordants
    readAlmts2 = args.readAlmts2 #workDir + "/aln2s.bam" #read 2 discordants
    workDir = args.workDir
    verbose = args.verbose
    bamfile = args.bamfile
    ignoreBED = args.ignoreBED
    ignoreChr = args.ignoreChr
    ignoreTIDs = []
    chrHash = {}

    RDL, mean_IL, disc_thresh, disc_thresh_neg = \
        calcMeanSig(bamfile, penalty_perc, disc_perc, disc_perc_neg, calc_thresh, as_calc_thresh, workDir)
    ignoreBuffer = 0*RDL
    print "Ignoring:", ignoreChr, ignoreBED
    if ignoreChr != "none":
        recordExcludeChr(ignoreChr, ignoreTIDs)
        print "Chromosomes", ignoreTIDs, "will be ignored."
    if ignoreBED != "none":
        print "Forming hash table for bad regions..."
        formExcludeHash(ignoreBuffer, chrHash, ignoreBED)
        print "Done."

    # read discordant alignments and write all possible discordant pairs to file.
    f1 = open(workDir + "/allDiscordants.us.txt","w")
    f2 = open(workDir + "/allDiscordants.up.us.txt","w")
    f3 = open(workDir + "/allDiscordants.txt","w")
    f4 = open(workDir + "/allDiscordants.up.txt","w")
    f3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("FragIndex", "Chr_l", "Chr_r", "L_pos", "R_pos", "Orient.","Small","MapQual"))
    f4.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("FragIndex", "Chr_l", "Chr_r", "L_pos", "R_pos", "Orient.","Small","MapQual"))
    currentFrag = 1
    print "Entering writeDiscordantFragments.py read loop..."
    for (q1,aln1s),(q2,aln2s) in izip(readNextReadAlignments(readAlmts1),readNextReadAlignments(readAlmts2)):

        if verbose and (currentFrag % 100000) == 0:
                print "Fragment", currentFrag, "analyzed."
        assert q1[:-2] == q2[:-2]
        #start_fd = clock()
        dList1, dList2 = formDiscordant(aln1s, aln2s, permutation_thresh, 
            map_thresh, as_thresh, nMatchPct_thresh, nMatch_relative_thresh, \
            as_relative_thresh, disc_thresh, disc_thresh_neg, big_num, mean_IL, chrHash, ignoreTIDs)
        for item in dList1:
             f1.write("%s %s\n" %(currentFrag, item))
        for item in dList2:
             f2.write("%s %s\n" %(currentFrag, item))
        #end_fd = clock()
        currentFrag = currentFrag + 1
    print "Done with read loop"

    f1.close()
    f2.close()
    f3.close()
    f4.close()
