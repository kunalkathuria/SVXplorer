# Add Split Reads to existing PE clusters
# This would have been faster if done in above method but pysam has no direct means of identifying split reads

import pysam
import sys
from random import randint

BAM_FILE = sys.argv[1] # "~/heterozygous/alignments/id.ns.splitters.bam"
RISK = int(sys.argv[2]) # 0 # risk some imprecise breakpoint locations and gain many precise ones
SLOP = int(sys.argv[3]) #20
NDisjEntries = 20 # No point making this a user-input, too specific
RAN_MIN = pow(10,8)
RAN_INT = pow(10,6)
SR_MARGIN = int(sys.argv[4]) #250

class OverlappingCluster(object):

    # at most 4 clusters will overlap
    bp1 = -1
    bp2 = -1
    bp3 = -1
    nclusters = 0
    count = 1
    typeO = None

    # it is an insertion cluster if overlapping point has had both orientation of reads
    bp1f = 0
    bp1r = 0
    bp1tid = -1
    bp2tid = -1
    bp3tid = -1
    support = None

    def __init__(self):
           self.clusterNums = [] # original cluster numbers that support this variant

    def __str__(self):
        return "%s %s %s %s %s %s %s %s %s" % (self.typeO, self.bp1, self.bp2,self.bp3, self.bp1tid, self.bp2tid, self.bp3tid, self.count, self.support)

def file_len(f):
    
    for i, l in enumerate(f):
        pass
    return i + 1

if __name__ == "__main__":

    fp1 = open("../results/text/All_Variants.txt","r")
    fp2 = open("../results/text/VariantMap.txt","r")
    fp3 = open("../results/text/All_Variants_SR.txt","w")
    fp4 = open("../results/text/VariantMap_SR.txt","w")
    fp5 = open("../results/text/updatedBPs_v5.txt","w")

    Claimed = {}
    SRFrag = 0
    matchingFrags = {}
    newbp1 = {}
    newbp2 = {}
    newbp3 = {}

    bamfile = pysam.Samfile(BAM_FILE,"rb")

    bp1_cons = -1 # bp numbers that are not skipped. 
    bp1_cons_tid = -1
    bp2_cons = -1
    bp2_cons_tid = -1  
    # SR-list loop first b/c 1) 1 match saves time as most variants do not occur together. Thre may be many alignments for same bp.
    # Either way, no way to decide which of multiple matches is really better
    # $ Can make it so if more than 1 match occurs for same SR, then do not update that bp location
    
    while True:

        try:
                
            sr1 = bamfile.next()
            sr2 = bamfile.next()
                    
        except StopIteration:
            break
    
        if sr1.qname == sr2.qname:

                SRFrag-=1

                if SRFrag % (-100) == 0:

                    print SRFrag

                sr_bp1 = sr1.reference_start
                sr_bp2 = sr2.reference_start
                sr_bp1_tid = sr1.reference_name
		sr_bp2_tid = sr2.reference_name

                if sr_bp1_tid == bp1_cons_tid and sr_bp2_tid == bp2_cons_tid and abs(sr_bp1 - bp1_cons) < SR_MARGIN and abs(sr_bp2 - bp2_cons) < SR_MARGIN:

                    continue

                bp1_cons = sr_bp1
                bp1_cons_tid = sr_bp1_tid
		bp2_cons = sr_bp2
		bp2_cons_tid = sr_bp2_tid

                if sr_bp1 < sr_bp2:
                    
                    minsr = sr1
                    maxsr = sr2
                    
                else:
                    
                    maxsr = sr1
                    minsr = sr2

		sr_bp1_tid = minsr.reference_name
                sr_bp2_tid = maxsr.reference_name

                if minsr.query_alignment_start > maxsr.query_alignment_start:

                    swap = 1

                    if minsr.is_reverse == maxsr.is_reverse:

                    # unless both reads are inverted and swapped -- unlikely
                        sr_bp1 = minsr.reference_start
                        sr_bp2 = maxsr.reference_end
                    

                    # $ Need to see if want to handle inverted copy-paste (may involve swap). Hard to code this criterion. Can take a risk of bad bp's for some variants if
                    # many will be caught correctly. Can avoid using slop in this particular case, perhaps. Okay, done. Small SLOP. Use "RISK".
                    # If no risk, won't support any inversions either, which is fine.
                    
                    elif RISK:

                        sr_bp1 = minsr.reference_end
                        sr_bp2 = maxsr.reference_end

                        
                else:
 
                    swap = 0

                    if minsr.is_reverse == maxsr.is_reverse:
                        
                        sr_bp1 = minsr.reference_end
                        sr_bp2 = maxsr.reference_start

                    # Taking same risk as above
                    elif RISK:

                        sr_bp1 = minsr.reference_start
                        sr_bp2 = maxsr.reference_start

                
                for line in fp1:

                    line_split = line.split()

                    if line_split[0]=="\n":
                        break

                    varNum = int(line_split[0])
                    tid_1 = line_split[1]
                    tid_2 = line_split[5]
                    tid_3 = line_split[8]
                    
                    pe_b1_m = int(line_split[3])
                    pe_b1_x = int(line_split[4])
                    
                    pe_b2_m = int(line_split[6])
                    pe_b2_x = int(line_split[7])

                    
                    pe_b3_m = int(line_split[9])
                    pe_b3_x = int(line_split[10])

                    match = 0

                    if line_split[1] == "DEL":

                        if (sr_bp1_tid == tid_1 and sr_bp1 > pe_b1_m and sr_bp1 < pe_b1_x) or abs(pe_b1_m - sr_bp1) < SLOP or abs(pe_b1_x - sr_bp2) < SLOP:

                            # 1st condition should be ensured below
                            if (sr_bp2_tid == tid_2 and sr_bp2 > pe_b2_m and sr_bp2 < pe_b2_x) or abs(pe_b2_m - sr_bp2) < SLOP or abs(pe_b2_x - sr_bp2) < SLOP:
                        
                                if not swap and minsr.is_reverse == maxsr.is_reverse:
                    
                                    match = 1
                                    #print "deletion match"

                    elif line_split[1] == "TD":

                        if (sr_bp1_tid == tid_1 and sr_bp1 > pe_b1_m and sr_bp1 < pe_b1_x) or abs(pe_b1_m - sr_bp1) < SLOP or abs(pe_b1_x - sr_bp2) < SLOP:

                            if (sr_bp2_tid == tid_2 and sr_bp2 > pe_b2_m and sr_bp2 < pe_b2_x) or abs(pe_b2_m - sr_bp2) < SLOP or abs(pe_b2_x - sr_bp2) < SLOP:

                                if swap and minsr.is_reverse == maxsr.is_reverse:

                                    match = 1
                                    #print "TD match"
                            
                    elif line_split[1] == "INV" and RISK:

                        if (sr_bp1_tid == tid_1 and sr_bp1 > pe_b1_m and sr_bp1 < pe_b1_x) or abs(pe_b1_m - sr_bp1) < SLOP or abs(pe_b1_x - sr_bp2) < SLOP:

                            if (sr_bp2_tid == tid_2 and sr_bp2 > pe_b2_m and sr_bp2 < pe_b2_x) or abs(pe_b2_m - sr_bp2) < SLOP or abs(pe_b2_x - sr_bp2) < SLOP:

                                if not swap and minsr.is_reverse != maxsr.is_reverse:
                            
                                    match = 1
                                    #print "INV match"
                                    
                    # One match per variant only, no point being redundant
                    if match and varNum not in newbp1:

                        newbp1[varNum] = sr_bp1
                        newbp2[varNum] = sr_bp2

                    # Do unless all 3 bp's were already updated by another fragment
                    # $ Remove this condition if want to count all frags to go in set cover and to set newbp's every time.
                    if not (varNum in newbp1 and varNum in newbp2 and varNum in newbp3):

                        if ((line_split[1] == "INS" or line_split[1] == "INS_C_P") and minsr.is_reverse == maxsr.is_reverse) or ((line_split[1] == "INS_I" or line_split[1] == "INS_C_I_P") and minsr.is_reverse != maxsr.is_reverse) :

                            if (sr_bp1_tid == tid_1 and sr_bp1 > pe_b1_m and sr_bp1 < pe_b1_x) or abs(pe_b1_m - sr_bp1) < SLOP or abs(pe_b1_x - sr_bp1) < SLOP:

                                if (sr_bp2_tid == tid_2 and sr_bp2 > pe_b2_m and sr_bp2 < pe_b2_x) or abs(pe_b2_m - sr_bp2) < SLOP or abs(pe_b2_x - sr_bp2) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp1
                                    newbp2[varNum] = sr_bp2

                                elif (sr_bp2_tid == tid_3 and sr_bp2 > pe_b3_m and sr_bp2 < pe_b3_x) or abs(pe_b3_m - sr_bp2) < SLOP or abs(pe_b3_x - sr_bp2) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp1
                                    newbp3[varNum] = sr_bp2

                            elif (sr_bp2_tid == tid_2 and sr_bp2 > pe_b1_m and sr_bp2 < pe_b1_x) or abs(pe_b1_m - sr_bp2) < SLOP or abs(pe_b1_x - sr_bp2) < SLOP:

                                if (sr_bp1_tid == tid_1 and sr_bp1 > pe_b2_m and sr_bp1 < pe_b2_x) or abs(pe_b2_m - sr_bp1) < SLOP or abs(pe_b2_x - sr_bp1) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp2
                                    newbp2[varNum] = sr_bp1

                                elif (sr_bp1_tid == tid_3 and sr_bp1 > pe_b3_m and sr_bp1 < pe_b3_x) or abs(pe_b3_m - sr_bp1) < SLOP or abs(pe_b3_x - sr_bp1) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp2
                                    newbp3[varNum] = sr_bp1

                        elif (line_split[1] == "INS_C" and minsr.is_reverse == maxsr.is_reverse) or (line_split[1] == "INS_C_I" and minsr.is_reverse != maxsr.is_reverse):

                            if (sr_bp1_tid == tid_1 and sr_bp1 > pe_b1_m and sr_bp1 < pe_b1_x) or abs(pe_b1_m - sr_bp1) < SLOP or abs(pe_b1_x - sr_bp1) < SLOP:

                                if (sr_bp2_tid == tid_2 and sr_bp2 > pe_b2_m and sr_bp2 < pe_b2_x) or abs(pe_b2_m - sr_bp2) < SLOP or abs(pe_b2_x - sr_bp2) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp1
                                    newbp2[varNum] = sr_bp2

                                elif (sr_bp2_tid == tid_3 and sr_bp2 > pe_b3_m and sr_bp2 < pe_b3_x) or abs(pe_b3_m - sr_bp2) < SLOP or abs(pe_b3_x - sr_bp2) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp1
                                    newbp3[varNum] = sr_bp2

                            elif (sr_bp2_tid == tid_2 and sr_bp2 > pe_b1_m and sr_bp2 < pe_b1_x) or abs(pe_b1_m - sr_bp2) < SLOP or abs(pe_b1_x - sr_bp2) < SLOP:

                                if (sr_bp1_tid == tid_1 and sr_bp1 > pe_b2_m and sr_bp1 < pe_b2_x) or abs(pe_b2_m - sr_bp1) < SLOP or abs(pe_b2_x - sr_bp1) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp2
                                    newbp2[varNum] = sr_bp1

                                elif (sr_bp1_tid == tid_3 and sr_bp1 > pe_b3_m and sr_bp1 < pe_b3_x) or abs(pe_b3_m - sr_bp1) < SLOP or abs(pe_b3_x - sr_bp1) < SLOP:

                                    match = 1
                                    newbp1[varNum] =  sr_bp2
                                    newbp3[varNum] = sr_bp1
                                    
                            elif (sr_bp1_tid == tid_3 and sr_bp1 > pe_b3_m and sr_bp1 < pe_b3_x) or abs(pe_b3_m - sr_bp1) < SLOP or abs(pe_b3_x - sr_bp1) < SLOP:

                                if (sr_bp2_tid == tid_2 and sr_bp2 > pe_b2_m and sr_bp2 < pe_b2_x) or abs(pe_b2_m - sr_bp2) < SLOP or abs(pe_b2_x - sr_bp2) < SLOP:

                                    match = 1
                                    newbp3[varNum] =  sr_bp1
                                    newbp2[varNum] = sr_bp2

                                elif (sr_bp2_tid == tid_1 and sr_bp2 > pe_b1_m and sr_bp2 < pe_b1_x) or abs(pe_b1_m - sr_bp2) < SLOP or abs(pe_b1_x - sr_bp2) < SLOP:

                                    match = 1
                                    newbp3[varNum] =  sr_bp1
                                    newbp1[varNum] = sr_bp2

                            elif (sr_bp2_tid == tid_3 and sr_bp2 > pe_b3_m and sr_bp2 < pe_b3_x) or abs(pe_b3_m - sr_bp2) < SLOP or abs(pe_b3_x - sr_bp2) < SLOP:

                                if (sr_bp1_tid == tid_2 and sr_bp1 > pe_b2_m and sr_bp1 < pe_b2_x) or abs(pe_b2_m - sr_bp1) < SLOP or abs(pe_b2_x - sr_bp1) < SLOP:

                                    match = 1
                                    newbp3[varNum] =  sr_bp2
                                    newbp2[varNum] = sr_bp1

                                elif (sr_bp1_tid == tid_1 and sr_bp1 > pe_b1_m and sr_bp1 < pe_b1_x) or abs(pe_b1_m - sr_bp1) < SLOP or abs(pe_b1_x - sr_bp1) < SLOP:

                                    match = 1
                                    newbp3[varNum] =  sr_bp2
                                    newbp1[varNum] = sr_bp1

   
                    if match:
                       
			print "Match", line
			
                        if SRFrag not in Claimed:
                            Claimed[SRFrag] = 1
			
			if varNum not in matchingFrags:
				matchingFrags[varNum] = [SRFrag]
			else:
                        	matchingFrags[varNum].append(SRFrag)

                        # There are no secondary alignments reported for split reads here. As a consequence, no ambiguous PE mappings will be picked/supported.
                        break
                
		#print matchingFrags    
                fp1.seek(0)
            
    # print new variant map
    fp1.seek(0)

    debug = 0
    for line3 in fp2:
        
        line3_split = line3.split()
        varNum = int(line3_split[0])
        #print "In loop", varNum
 
        for line in fp1:
            
            line_split = line.split()

            if varNum in matchingFrags:
		
		#print "Changing variant bps", varNum
		
                if varNum in newbp1:
                    
                    line_split[3] = str(newbp1[varNum])
                    line_split[4] = str(newbp1[varNum]+1)

                if varNum in newbp2:
                    
                    line_split[6] = str(newbp2[varNum])
                    line_split[7] = str(newbp2[varNum]+1)

                if varNum in newbp3:
                    
                    line_split[9] = str(newbp3[varNum])
                    line_split[10] = str(newbp3[varNum]+1)

                line_split[11] = line_split[11] + "_SR"

            
            temp2 = " ".join(line_split)
            fp3.write("%s\n" %temp2)

	    #print "Writing", temp2	
	    break

        temp = " ".join(line3_split)
    
        fp4.write("%s" %temp)

        # This assumes that order of variants is same in All_Variants and VariantMap -- which should be true
        # print line3_split[0]
        if varNum in matchingFrags:
                for item in matchingFrags[varNum]:
                    fp4.write(" %s" %item)
                
                # Add random fragments to boost SR support to set cover. Especially for insertions as done above, as only 1 or 2 fragments will match in total.
                # Create disjointness so that this is picked among ambiguous mappings etc.
                # Remove if secondary/multiple alignments for SRs are included, as then all the ambiguous clusters will be picked.
                # "SRFrag" numbered by alignment, not read.
                for _ in range(NDisjEntries):
                
                    temp = randint(RAN_MIN, RAN_MIN + RAN_INT)
                    #print "Writing", h, temp
                    fp4.write(" %s" %temp)

        fp4.write("\n")

          
                                                                                     
    bamfile.close()
    
