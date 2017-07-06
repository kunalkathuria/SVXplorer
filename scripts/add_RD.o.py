# Add RD variants to existing PE variants

import pysam
import sys
from random import randint
import math

NDisjEntries = 20
RAN_MIN = pow(10,8)
RAN_INT = pow(10,6)
RO_THRESH = float(sys.argv[1]) #.7
CNV_FILE = sys.argv[2]

class OverlappingCluster(object):

    # at most 4 clusters will overlap
    # $ move to init
    bp1_start = -1
    bp2_start = -1
    bp3_start = -1
    bp1_end = -1
    bp2_end = -1
    bp3_end = -1
    count = 2
    typeO = None

    # it is an insertion cluster if overlapping point has had both orientation of reads
    bp1f = 0
    bp1r = 0
    bp1tid = -1
    bp2tid = -1
    bp3tid = -1

    def __init__(self):
        self.clusterNums = [] # original cluster numbers that support this variant
        self.support = "PE"
        self.complete = 0 # All clusters that can match with this variant have arrived and matched
        self.varNum = None

    def __str__(self):
        #return "%s %s %s %s %s %s %s %s %s" % (self.typeO, self.bp1, self.bp2,self.bp3, self.bp1tid, self.bp2tid, self.bp3tid, self.count, self.support)
        return "%s %s %s %s %s %s %s %s %s %s %s %s" % (self.typeO, self.bp1tid, self.bp1_start, self.bp1_end, self.bp2tid , self.bp2_start, self.bp2_end, self.bp3tid, self.bp3_start, self.bp3_end, self.support, self.count)

def file_len(f):
    
    for i, l in enumerate(f):
        pass
    return i + 1

def ROCalc(min1,max1,min2,max2):

    if min1 <= min2 and max1 >= min2 and max1 <= max2:

        RO = (max1 - min2)/(max2-min2)

    elif min1 <= min2 and max1 >= max2:

        RO = 1

    elif min1 >= min2 and max1 <= max2:

        RO = (max1 - min1)/(max2-min2)

    elif min1 >= min2 and min1 <= max2 and max1 >= max2:

        RO = (max2 - min1)/(max2-min2)
    else:
        RO = 0

    return RO

if __name__ == "__main__":

    # extract out INV from All Variants and store rest in inDels.txt
    f2 = open("../results/text/inDels_S.txt","r")
    fp2 = open("../results/text/VariantMap.txt","r")
    fp3 = open("../results/text/All_Variants_RD.txt","w")
    fp4 = open("../results/text/VariantMap_RD.txt","w")
    f1 = open(CNV_FILE,"r")
    fp6 = open("../results/text/All_Variants.txt","r")

    
    Claimed = {}
    RDfrag = 100000000
    counter = -1
    RDcounter = -1
    debug = 0

    matchingFrags = {}
    bpSwap = {}

    line2 = 1
    prev_pos = 0
    lastTidStart = 0

    for line_cnv in f1:

        RDfrag+=1
        RDcounter+=1

        #print RDcounter

        counter = -1

        line_cnv_split = line_cnv.split()

        min1 = int(line_cnv_split[1])
        max1 = int(line_cnv_split[2])

        tid1 = line_cnv_split[0]
        CN = int(line_cnv_split[3])

        tidFound = 0
        match = 0

        if line2 != '' and line2 !='\n':
            line2 = 1
        else:
            line2= 0

        curr_pos = f2.tell()

        while line2:

            counter+=1
            
            line2 = f2.readline()
            if len(line2) == 0 or line2 == '\n':
                break

            line2_split = line2.split()
            
            variantNum = int(line2_split[0])
            match = 0

            tid2 = line2_split[5]
            #min2 = int(line2_split[3])

            prev_pos = curr_pos
            curr_pos = f2.tell()

            if not tidFound and tid1 == tid2:
                tidFound = 1
                lastTidStart = prev_pos

            #print "Comparing", line1, "and", line2

            #Compare all if TID same as cannot a priori set a SLOP/margin cutoff. Criterion is reciprocal overlap.
            if tidFound and tid1 != tid2:

                break

            elif tid1 == tid2:

                print "Comparing", line_cnv, "and", line2

                #print "4"
                tidFound = 1
                #tidStartPos = f2.tell()

                ## SEE IF MATCH OCCURS
                if line2_split[1] == "DEL" and CN < 2:

                    min2 = int(line2_split[3])
                    max2 = int(line2_split[6])
                    tid2 = line2_split[2]

                    #print min1, max1, min2, max2

                    # change to within range?
                    #if tid1 == tid2 and (min1 < min2 or abs(min1-min2) < SLOP) and (max1 > max2 or abs(max1-max2) < SLOP):

                    RO = ROCalc(min1,max1,min2,max2)

                    if RO > RO_THRESH:
                        match = 1

                    if RO > RO_THRESH:
                        match = 1
                    #print "deletion match"

                elif line2_split[1] == "TD" and CN > 2:

                    min2 = int(line2_split[3])
                    max2 = int(line2_split[6])
                    tid2 = line2_split[2]

                    RO = ROCalc(min1,max1,min2,max2)
                    
                    if RO > RO_THRESH:
                        match = 1
                        
                elif line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I" and CN > 2:

                    # BP1 and BP3 may be switched because PE cannot tell. All 3 will be on same chromosome in that case and above label implies that.

                    bp2 = int(line2_split[6])
                    bp3 = int(line2_split[9])
                    bp1 = int(line2_split[3])
                    tid2 = line2_split[2]

                    min2 = bp1
                    max2 = bp2

                    RO = ROCalc(min2,max2,min1,max1)

                    if RO > RO_THRESH:

                        bp1_n_2 = 1

                    else:                       
                        bp1_n_2 = 0


                    min2 = bp2
                    max2 = bp3

                    RO = ROCalc(min2,max2,min1,max1)

                    if RO > RO_THRESH:

                        bp2_n_3 = 1

                    else:                       
                        bp2_n_3 = 0

                    if bp1_n_2 and not bp2_n_3:

                        bpSwap[variantNum] = 1
                        match = 1

                    if bp2_n_3 and not bp1_n_2:

                        match = 1

                # Other insertions
                elif line2_split[1][0:2] == "IN" and CN > 2:

                    # order guaranteed
                    min2 = int(line2_split[6])
                    max2 = int(line2_split[9])
                    tid2 = line2_split[5]

                    RO = ROCalc(min2,max2,min1,max1)

                    if RO > RO_THRESH:
                        match = 1

                        #print "INS match"

                if match:
                    #print "Match", counter, line_split[1]
                    print "Match", line_cnv, "and", line2
                    
                    if not RDcounter in Claimed:
                        Claimed[RDcounter] = 1

                    if variantNum in matchingFrags:
                        matchingFrags[variantNum].append(RDfrag)
                    else:
                        matchingFrags[variantNum] = [RDfrag]
                    
                    
                    # With RD more than one match is not indicated for haploid. Remove for diploid or otherwise possible conflicts.
                    #break


        
        line2 = 1
        f2.seek(lastTidStart)


    # print new variant map

    for h,line3 in enumerate(fp2):
        
        line3_split = line3.split()
   
        varMapLine = " ".join(line3_split)
        
        # check if order of variants in both files same: YES
        for line in fp6:
            
            line_split = line.split()
            variantNum = int(line_split[0])

            if variantNum in matchingFrags:
                
                if line_split[1] == "INS_C" or line_split[1] == "INS_C_I":
                    
                    line_split[1]+= "_P"
                    
                    if variantNum in bpSwap:
                        
                        temp1 = line_split[3]
                        temp2 = line_split[4]
                        line_split[3] = line_split[9]
                        line_split[4] = line_split[10]
                        line_split[9] = temp1
                        line_split[10] = temp2
                        
                line_split[11] = line_split[11] + "_RD"

            varLine = " ".join(line_split)
            fp3.write("%s\n" %varLine)
            break


        fp4.write("%s" %varMapLine)

        # This assumes that order of variants is same in All_Variants and VariantMap -- which should be true
        #print line3_split[0]

        if variantNum in matchingFrags:
            for item in matchingFrags[variantNum]:
                fp4.write(" %s" %item)

            # create disjointness so that this is picked among ambiguous mappings if only one with RD support
            # can have the same set of random numbers if variant bp's same from PE mapping BUT v unlikely that both will have RD support
            for _ in range(NDisjEntries):
                
                temp = randint(RAN_MIN, RAN_MIN + RAN_INT)
                #print "Writing", h, temp
                fp4.write(" %s" %temp)

        fp4.write("\n")
            
        
   
    f1.close()
    f2.close()
    fp2.close()
    fp3.close()
    fp4.close()
    fp6.close()
        
