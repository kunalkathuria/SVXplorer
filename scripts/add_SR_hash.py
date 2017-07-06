# Add Split Reads to existing PE clusters
# This would have been faster if done in above method but pysam has no direct means of identifying split reads

import pysam
import sys
from random import randint

BAM_FILE = sys.argv[1] # "~/heterozygous/alignments/id.ns.splitters.bam"
RISK = int(sys.argv[2]) # 0 # risk some imprecise breakpoint locations and gain many precise ones
NDisjEntries = 4 # SET = DTHRESH
RAN_MIN = pow(10,8)
RAN_INT = pow(10,6)
SR_SLOP = int(sys.argv[3])
SR_MARGIN = int(sys.argv[4]) #250
varHash = {}

class SRVar(object):

	def __init__(self):

		self.bp1 = -1
		self.bp2 = -1
		self.bp3 = -1

class Variant(object):

    #$ add inv, ins also (bp 3)
    def __init__(self):
    	self.bp = -1

    	self.tid = -1
    	self.tid_2 = -1

    def __str__(self):

	return "%s %s %s" %(self.bp, self.tid, self.tid_2)

    def __hash__(self):
        return hash((self.bp, self.tid, self.tid_2))

    def __eq__(self, other):
        return (self.bp, self.tid, self.tid_2) == (other.bp, other.tid, other.tid_2)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

class Bound(object):

	def __init__(self):
		
		self.bp2_1 = -1
		self.bp2_2 = -1
		self.bp3_1 = -1
		self.bp3_2 = -1
		self.typeO = -1
		self.num = -1

	def __str__(self):
		
		return "%s %s %s %s %s %s" %(self.bp2_1, self.bp2_2, self.bp3_1, self.bp3_2, self.typeO, self.num)

def mapSVtoNum(SV_type):

	if SV_type== "DEL":
		return 0
	elif SV_type == "TD":
		return 1
	elif SV_type == "INV":
		return 2
	elif SV_type == "INS":
		return 3
		# $ implement INS_C separately from other insertions due to extra DEL cluster (one more for loop) and bp3vs1 uncertainty
	else:
		return -1

def formHash(fp1):

	print "Forming hash..."

        for line in fp1:
                line_s = line.split()
		bounds = Bound()		
		bounds.num = int(line_s[0])			
		bounds.typeO= mapSVtoNum(line_s[1])

		if bounds.typeO == -1:
			continue
		if bounds.typeO == 3:

			bounds.bp3_1 = int(line_s[9]) - SR_SLOP
			bounds.bp3_2 = int(line_s[10]) + SR_SLOP

                bounds.bp2_1 = int(line_s[6]) - SR_SLOP
                bounds.bp2_2 = int(line_s[7]) + SR_SLOP

		# hash all values within bp margin		
		for x in range(int(line_s[3])-SR_SLOP, int(line_s[4])+ SR_SLOP):
			temp = Variant()
			temp.tid = line_s[2]
			if line_s[5] != temp.tid:
                                temp.tid_2 = line_s[5]	
			temp.bp = x
			varHash[temp] = bounds
                	#print "Added", temp
	print "Done."

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

    formHash(fp1)

    #for item in varHash:
	#print "In hash:", item

    Claimed = {}
    SRFrag = 0
    matchingFrags = {}
    newbp1 = {}
    newbp2 = {}
    newbp3 = {}
    newVars=1
    bamfile = pysam.Samfile(BAM_FILE,"rb")

    bp1_cons = -1 # bp numbers that are not skipped. 
    bp1_cons_tid = -1
    bp2_cons = -1
    bp2_cons_tid = -1  
    
    while True:

        try:
                
            sr1 = bamfile.next()
            sr2 = bamfile.next()
                    
        except StopIteration:
            break
   
	varType = -1
 
        if sr1.qname == sr2.qname:

                SRFrag-=1

                if SRFrag % (-1000) == 0:

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

                    # unless both reads are inverted and swapped -- less likely
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

			if minsr.is_reverse:

	                        sr_bp1 = minsr.reference_start
        	                sr_bp2 = maxsr.reference_start
			else:

				sr_bp1 = minsr.reference_end
                                sr_bp2 = maxsr.reference_end

                
	else:
		print "Please check if split-read file is name-sorted."
		continue

	match = 0
	temp = Variant()
	temp.bp = sr_bp1
	temp.tid = sr_bp1_tid
	#print "Formed sr:", temp

	if sr_bp1_tid != sr_bp2_tid:	                                                      

		temp.tid_2 = sr_bp2_tid

	if temp in varHash and (varHash[temp].bp2_1 < sr_bp2 < varHash[temp].bp2_2 or varHash[temp].bp3_1 < sr_bp2 < varHash[temp].bp3_2):

		#print "Found"
		varNum = varHash[temp].num
		varType = varHash[temp].typeO

		if varNum not in matchingFrags:
                        newbp = SRVar()
                        newbp.bp1 = sr_bp1
                        newbp.bp2 = sr_bp2

	# if not found, try other side of sr
	else:

		temp.bp = sr_bp2

		if temp.tid_2 != -1:

			swap = temp.tid_2
			temp.tid_2 = temp.tid
			temp.tid = swap

		if temp in varHash and (varHash[temp].bp2_1 < sr_bp2 < varHash[temp].bp2_2 or varHash[temp].bp3_1 < sr_bp2 < varHash[temp].bp3_2):

			varNum = varHash[temp].num
			varType = varHash[temp].typeO
			
			if varNum not in matchingFrags:
				newbp = SRVar()
				newbp.bp1 = sr_bp2
				newbp.bp2 = sr_bp1

	#check DEL
	if varType == 0 and not swap and minsr.is_reverse == maxsr.is_reverse:

		match = 1
	#check TD
	elif varType == 1 and swap and minsr.is_reverse == maxsr.is_reverse:

		match = 1

	#check INV
	elif varType == 2 and not swap and minsr.is_reverse != maxsr.is_reverse:
	
		match = 1

	#check INS
	elif varType == 3 and minsr.is_reverse == maxsr.is_reverse:

		match = 1
		
		# Update bp3 of insertion if unset
		if varNum in matchingFrags and matchingFrags[varNum][0].bp3 == -1:

			temp = matchingFrags[varNum][0]
			if abs(sr_bp1 - temp.bp1) > SR_SLOP and abs(sr_bp1 - temp.bp2) > SR_SLOP:
			
				matchingFrags[varNum][0].bp3 = sr_bp1

			elif abs(sr_bp2 - temp.bp1) > SR_SLOP and abs(sr_bp2 - temp.bp2) > SR_SLOP:
		
				matchingFrags[varNum][0].bp3 = sr_bp2

			# insertion bp2 should be < bp3
                        if matchingFrags[varNum][0].bp3 < matchingFrags[varNum][0].bp2:
			
				swap = matchingFrags[varNum][0].bp3
				matchingFrags[varNum][0].bp3 = matchingFrags[varNum][0].bp2
				matchingFrags[varNum][0].bp2 = swap     
   
	if match:
                       
		#print "Match"
		#print SRFrag, varHash[temp].num, varHash[temp].type0, temp.tid, temp.tid_2, temp.bp, varHash[temp].bp2_1, varHash[temp].bp3_1
		
		if varNum not in matchingFrags:
			matchingFrags[varNum] = [newbp, SRFrag]
		else:
			matchingFrags[varNum].append(SRFrag)

                
        #else:
		
		#print "New possible variant seen.", sr_bp1_tid, sr_bp1, sr_bp2_tid, sr_bp2 
		# add entry to file at top
		#fp3.write("%s %s %s %s %s %s %s %s %s %s %s %s" %(newVars, "BND", sr_bp1_tid, sr_bp1, sr_bp1 + 1, sr_bp2_tid, sr_bp2, sr_bp2 + 1, -1,-1,-1, "SR",1)		
		#newVars+=1
	
    # print new variant map
    fp1.seek(0)


    for line3 in fp2:
        
        line3_split = line3.split()
        varNum = int(line3_split[0])
        #print "In loop", varNum
 
        for line in fp1:
            
            line_split = line.split()

            if varNum in matchingFrags:
                line_split[11] = line_split[11] + "_SR"
		line_split[3] = str(matchingFrags[varNum][0].bp1)
		line_split[4] = str(matchingFrags[varNum][0].bp1 + 1)
		line_split[6] = str(matchingFrags[varNum][0].bp2)
		line_split[7] = str(matchingFrags[varNum][0].bp2 + 1)
		if matchingFrags[varNum][0].bp3 != -1:
			line_split[9] = str(matchingFrags[varNum][0].bp3)
			line_split[10] = str(matchingFrags[varNum][0].bp3 + 1)
            
            temp2 = " ".join(line_split)
            fp3.write("%s\n" %temp2)

	    #print "Writing", temp2	
	    break

        temp = " ".join(line3_split)
    
        fp4.write("%s" %temp)

        # This assumes that order of variants is same in All_Variants and VariantMap -- which should be true
        # print line3_split[0]
        if varNum in matchingFrags:
                for item in matchingFrags[varNum][1:]:
                    fp4.write(" %s" %item)
                
        fp4.write("\n")

          
                                                                                     
    bamfile.close()
    
