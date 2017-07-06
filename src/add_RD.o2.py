# Add RD variants to existing PE variants

import pysam
import sys
from random import randint
import math

RO_THRESH = float(sys.argv[1]) #.7
CNV_FILE = sys.argv[2]
# depends on tool accuracy, but not of consequence if slop is used in comparison to truth set
RD_MARGIN = int(sys.argv[5])
RD_SLOP = int(sys.argv[6])
MAP_THRESH_RDU = int(sys.argv[7])

class Variant(object):

    #$ add inv, ins also (bp 3)
    def __init__(self):
        self.bp1 = -1
	self.bp2 = -1

        self.tid = -1

    def __hash__(self):
        return hash((self.bp1, self.bp2, self.tid))

    def __eq__(self, other):
        return (self.bp1, self.bp2, self.tid) == (other.bp1, other.bp2, other.tid)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

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
    RDfrag = int(sys.argv[4]) #100000000
    RDcounter = int(sys.argv[3]) #3000000
    debug = 0
    newRDVars = []
    matchingFrags = {}
    bpSwap = {}

    line2 = 1
    prev_pos = 0
    lastTidStart = 0
    chr_tag = 1
    #"chr" prefix added by cnvnator, perhaps others
    line_AV = fp6.readline()
    if line_AV.split()[2].find("chr") == -1:
	chr_tag = 0
    
    for line_cnv in f1:

        RDfrag+=1
        RDcounter+=1

        print "Done with CNV number", RDcounter

        counter = -1

        line_cnv_split = line_cnv.split()

        #min1 = int(line_cnv_split[1][line_cnv_split[1].find(":")+1:line_cnv_split[1].find("-")])
        #max1 = int(line_cnv_split[1][line_cnv_split[1].find("-")+1:])

        #tid1 = line_cnv_split[1][:line_cnv_split[1].find(":")]

	tid1 = line_cnv_split[0]
	min1 = line_cnv_split[1]
	max1 = line_cnv_split[2]

	if tid1.find("chr") != -1 and chr_tag == 0:
		tid1 = tid1[3:]

        CN = float(line_cnv_split[3])

	#print tid1, min1, max1, "CNV listing"
        tidFound = 0
        match = 0
	Anymatch = 1

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

                #print "Comparing", line_cnv, "and", line2

                #print "4"
                tidFound = 1
                #tidStartPos = f2.tell()

                ## SEE IF MATCH OCCURS
                if len(line2_split[1]) >=3 and line2_split[1][0:3] == "DEL" and CN < 1:

                    min2 = int(line2_split[3])
                    max2 = int(line2_split[6])
                    tid2 = line2_split[2]

                    #print min1, max1, min2, max2

                    # change to within range?
                    #if tid1 == tid2 and (min1 < min2 or abs(min1-min2) < SLOP) and (max1 > max2 or abs(max1-max2) < SLOP):

                    RO = ROCalc(min1,max1,min2,max2)

                    if RO > RO_THRESH:
                        match = 1

                    elif RO > 0:
                        Anymatch = 0

                    #print "deletion match"

                elif line2_split[1][0:2] == "TD" and CN > 1:

                    min2 = int(line2_split[3])
                    max2 = int(line2_split[6])
                    tid2 = line2_split[2]

                    RO = ROCalc(min1,max1,min2,max2)
                    
                    if RO > RO_THRESH:
                        match = 1
		    elif RO > 0:
                        Anymatch = 0
                        
                elif line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I" and CN > 1:

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

			if RO > 0:
                        	Anymatch = 0

                    min2 = bp2
                    max2 = bp3

                    RO = ROCalc(min2,max2,min1,max1)

                    if RO > RO_THRESH:

                        bp2_n_3 = 1

                    else:                       
                        bp2_n_3 = 0
			if RO > 0:
                        	Anymatch = 0

                    if bp1_n_2 and not bp2_n_3:

                        bpSwap[variantNum] = 1
                        match = 1

                    if bp2_n_3 and not bp1_n_2:

                        match = 1

                # Other insertions; can't handle SR-called DEL_INS here as only 1 source loc caught
                elif line2_split[1][0:2] == "IN" and CN > 1:

                    # order guaranteed
                    min2 = int(line2_split[6])
                    max2 = int(line2_split[9])
                    tid2 = line2_split[5]

                    RO = ROCalc(min2,max2,min1,max1)

                    if RO > RO_THRESH:
                        match = 1

                        #print "INS match"
		    elif RO > 0:
                        Anymatch = 0	
		
                if match:
                    #print "Match", counter, line_split[1]
                    #print "Match", line_cnv, "and", line2
                    
                    if not RDcounter in Claimed:
                        Claimed[RDcounter] = 1

                    if variantNum in matchingFrags:
                        matchingFrags[variantNum].append(RDfrag)
                    else:
                        matchingFrags[variantNum] = [RDfrag]
                    
                    
                    # With RD more than one match is not indicated for haploid. Remove comment for diploid or otherwise possible conflicts.
                    break

 	if not RDcounter in Claimed and Anymatch:
		#Use only those CNV calls with high nbins used
		#fb=open(str("../results/text/" + tid + ".bin.rdtmp"),"r")
		#for line_t in fb:
			#line_t_split = line_t.split()
			#if line_t_split[1] != "start" and int(line_t_split[1]) > min1 and int(line_t_split[1]) < max1:
				#print "New variant found.", line_cnv
				newRDVars.append([tid1, min1, max1, CN])
				#break
			#elif int(line_t_split[1]) > max1:
				#break
		#fb.close()
       
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

        fp4.write("\n")
    
    f=open("../results/text/All_Discords_P_S.txt","r")
    varHash = {}
    iobjects = []
    print "Forming RD hash..."
	
    for line in f:

                chr1 = line.split()[1]
                chr2 = line.split()[3]
		start = int(line.split()[2])
		stop = int(line.split()[4])
		FR = line.split()[5]
		mq = int(line.split()[-1])

		# hash all values within bp margin
                temp = Variant()
		temp.tid = chr1
		temp.bp1 = start
		temp.bp2 = stop

		#immutable hash objects -- preserve so memory doesn't get overwritten
		iobjects.append(temp)
		if mq > MAP_THRESH_RDU and chr1 == chr2 and FR[0] != FR[1] and temp not in varHash:
			varHash[temp] = FR
			#print "Added", temp

    for k,entry in enumerate(newRDVars):

	chrom = entry[0]
	start = int(entry[1])
	stop = int(entry[2])
	write = 0

	for x in range(start - RD_SLOP, start + RD_SLOP):

		for y in (stop - RD_SLOP, stop + RD_SLOP):

			temp = Variant()
			temp.tid = chrom
        	        temp.bp1 = x
               		temp.bp2 = y
			iobjects.append(temp)

			if temp in varHash:
				print "Supporting read for CNV call found..."
				write = 1
				break
		if temp in varHash:
			break

	if write and float(entry[3]) > 1 and varHash[temp] == "10":
		typeV = "INS"
		fp3.write("%s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(RDcounter+1+k, typeV, "-1", "-1", "-1", chrom, start-RD_MARGIN/2, start+RD_MARGIN/2, chrom, stop-RD_MARGIN/2, stop +RD_MARGIN/2, "RD", "1"))

	elif write and varHash[temp] == "01":
		typeV = "DEL"

		fp3.write("%s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(RDcounter+1+k, typeV, chrom, start-RD_MARGIN/2, start+RD_MARGIN/2, chrom, stop-RD_MARGIN/2, stop +RD_MARGIN/2, "-1", "-1", "-1", "RD", "1"))
	# the 2* below is to ensure that frag numbers supporting existing PE, SR variants are not repeated here. The actual number is inconsequential as long as unique.
	fp4.write("%s %s\n" %(RDcounter+1+k, 2*RDfrag+1+k))
   
    f1.close()
    f2.close()
    fp2.close()
    fp3.close()
    fp4.close()
    fp6.close()
        
