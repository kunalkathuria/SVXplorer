# Add Split Reads to existing PE clusters
# This would have been faster if done in above method but pysam has no direct means of identifying split reads

import pysam
import sys
from random import randint
import copy

BAM_FILE = sys.argv[1] # "~/heterozygous/alignments/id.ns.splitters.bam"
RISK = int(sys.argv[2]) # 0 # risk some imprecise breakpoint locations and gain many precise ones
NDisjEntries = 4 # SET = DTHRESH
RAN_MIN = pow(10,8)
RAN_INT = pow(10,6)
SR_SLOP = int(sys.argv[3])
SR_MARGIN = int(sys.argv[4])
MIN_CS_SR = int(sys.argv[5])
MAP_THRESH = int(sys.argv[6])
ignoreChr = sys.argv[7]
varHash = {}
SRVarHash = {}

class newSRVar(object):
	def __init__(self):

		self.l_orient = -1
		self.r_orient = -1
		self.swapped = -1
		self.bp2 = -1
		self.bp3 = -1
		self.count = 1
		self.support = []
		self.bp3tid = -1	
		self.typeO = -1
		self.tag = -1
		self.neighbor_tags = []
		self.hash_pair_tag = -1
		self.n_changes = 0
		self.write = -1

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

def swap_f(a,b):
	temp = a
	a = b
	b = temp

def transferSupport(variant1, variant2):
	temp = variant2
	# 1st fragment is shared
	for elem in variant2.support[1:]:
		variant1.support.append(elem)
		variant1.count+=1        	

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

def formHash(fp1, iobjects):

	print "Forming PE variant hash table..."

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

			#immutable hash objects -- preserve so memory doesn't get overwritten
			iobjects.append(temp)
			if temp not in varHash:
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

    # ensure hash objects immutable
    immutable_objects = []
    formHash(fp1, immutable_objects)

    #for item in varHash:
	#print "In hash:", item

    SRFrag = 0
    matchingFrags = {}
    newSRList = []
    bamfile = pysam.Samfile(BAM_FILE,"rb")

    bp1_cons = -1 # bp numbers that are not skipped. 
    bp1_cons_tid = -1
    bp2_cons = -1
    bp2_cons_tid = -1  
    counter = 1

    # split reads should be mapped, unique alignments  
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

		ignoreList = []
		fo=open(ignoreChr, "r")
		for entry in fo:
			ignoreList.append(entry)

                if sr1.reference_name in ignoreList or sr2.reference_name in ignoreList or sr1.mapping_quality < MAP_THRESH or sr2.mapping_quality < MAP_THRESH or (sr_bp1_tid == bp1_cons_tid and sr_bp2_tid == bp2_cons_tid and abs(sr_bp1 - bp1_cons) < SR_MARGIN and abs(sr_bp2 - bp2_cons) < SR_MARGIN):

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
		temp = Variant()
		temp.bp = sr_bp2
		temp.tid = sr_bp2_tid
		temp.tid_2 = -1

		if sr_bp1_tid != sr_bp2_tid:

			temp.tid_2 = sr_bp1_tid

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
	#elif varType == 0:
		#print sr_bp1_tid, sr_bp1, sr_bp2_tid, sr_bp2, swap, minsr.is_reverse, maxsr.is_reverse

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
			
				swapper = matchingFrags[varNum][0].bp3
				matchingFrags[varNum][0].bp3 = matchingFrags[varNum][0].bp2
				matchingFrags[varNum][0].bp2 = swapper     
   
	if match:
                       
		#print "Match"
		#print SRFrag, varHash[temp].num, varHash[temp].type0, temp.tid, temp.tid_2, temp.bp, varHash[temp].bp2_1, varHash[temp].bp3_1
		
		if varNum not in matchingFrags:
			matchingFrags[varNum] = [newbp, SRFrag]
		else:
			matchingFrags[varNum].append(SRFrag)

                
        else:
		#Look in SR hash -- store only breakpoints up to +/- SR_SLOP/2 (plenty) and then tally up orientation
		#Remove orientation from temp if streamlining PE hash detection
		SRCheck = 0
		claimed = 0
		temp = Variant()
		temp.bp = sr_bp1
	        temp.tid = sr_bp1_tid
       		#print "Formed sr:", temp

		# bp1 is "left" bp, bp2 is "right" bp
		temp.tid_2 = sr_bp2_tid
		l_orient = minsr.is_reverse
                r_orient = maxsr.is_reverse

		if temp in SRVarHash:
		
			SRCheck = 1
			other_bp = sr_bp2
			other_bp_tid = sr_bp2_tid
		else:
			# hashing -- safe to declare new object
			temp = Variant()
			temp.bp = sr_bp2	
			temp.tid = sr_bp2_tid
			temp.tid_2 = sr_bp1_tid

			if temp in SRVarHash:

	                        SRCheck = 1
				other_bp = sr_bp1
				other_bp_tid = sr_bp1_tid
		if SRCheck:

			#note: INS_I can be off by ~.5 RDL or less due to ambiguity of swapped alignment: whether ref start or end alignment point is breakpoint
                        if SRVarHash[temp].typeO == "DEL_INS":

				if swap and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2) and not ( (l_orient != SRVarHash[temp].l_orient and SRVarHash[temp].r_orient == r_orient) or (l_orient == SRVarHash[temp].l_orient and SRVarHash[temp].r_orient != r_orient) ):
					SRVarHash[temp].typeO = "INS" 
                       	
					SRVarHash[temp].bp3 = other_bp
					SRVarHash[temp].bp3tid = other_bp_tid

					if l_orient != r_orient:
						SRVarHash[temp].typeO = "INS_I"
			
					SRVarHash[temp].count+=1
	                                SRVarHash[temp].support.append(SRFrag)

					SRVarHash[temp].n_changes+=1
					claimed = 1

				# cannot update type based on this info	
				elif not swap and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient == r_orient:

					SRVarHash[temp].count+=1
	                                SRVarHash[temp].support.append(SRFrag)
				
				#INS_C
				#elif not swap and l_orient == r_orient and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2):
        		#TD_I
                        elif SRVarHash[temp].typeO == "TD_I":
				# cannot update type $ DO TD_I_I if inverted (included in this case)
				if swap and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and not ( (l_orient != SRVarHash[temp].l_orient and SRVarHash[temp].r_orient == r_orient) or (l_orient == SRVarHash[temp].l_orient and SRVarHash[temp].r_orient != r_orient) ):

	                                SRVarHash[temp].count+=1
        	                        SRVarHash[temp].support.append(SRFrag)

				elif not swap and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2) and not ( (l_orient != SRVarHash[temp].l_orient and SRVarHash[temp].r_orient == r_orient) or (l_orient == SRVarHash[temp].l_orient and SRVarHash[temp].r_orient != r_orient) ):

					SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
					SRVarHash[temp].typeO = "INS"
					SRVarHash[temp].bp3 = other_bp
					SRVarHash[temp].bp3tid = other_bp_tid

                                        if l_orient != r_orient:
                                                SRVarHash[temp].typeO = "INS_I"

					SRVarHash[temp].n_changes+=1

                        elif SRVarHash[temp].typeO[:3] == "INV" and l_orient != r_orient:
			   if SRVarHash[temp].typeO == "INV_POSS":
				#cannot update type
				if not swap and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient == SRVarHash[temp].l_orient:
					SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
				elif not swap and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient != SRVarHash[temp].l_orient:
					SRVarHash[temp].typeO = "INV"
                                	SRVarHash[temp].count+=1
                                	SRVarHash[temp].support.append(SRFrag)
					SRVarHash[temp].n_changes+=1
				elif swap and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2) and l_orient != SRVarHash[temp].l_orient:
					SRVarHash[temp].typeO = "INS_I"
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
					SRVarHash[temp].bp3 = other_bp
					SRVarHash[temp].bp3tid = other_bp_tid
					SRVarHash[temp].n_changes+=1

			   elif SRVarHash[temp].typeO == "INV":
                                
                                if not swap and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2:
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)

			elif SRVarHash[temp].typeO == "INS":

				if l_orient == r_orient and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or (SRVarHash[temp].bp3 != -1 and SRVarHash[temp].bp3-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp3+SR_SLOP/2):
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
					if SRVarHash[temp].bp3 == -1:
                                                SRVarHash[temp].bp3 = other_bp
						SRVarHash[temp].bp3tid = other_bp_tid

			elif SRVarHash[temp].typeO == "INS_I":

                                if l_orient != r_orient and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or (SRVarHash[temp].bp3 != -1 and SRVarHash[temp].bp3-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp3+SR_SLOP/2):
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
					if SRVarHash[temp].bp3 == -1:
						SRVarHash[temp].bp3 = other_bp
						SRVarHash[temp].bp3tid = other_bp_tid

		# for translocations, need to check bp 2 so check all elements
		elif temp.tid == temp.tid_2 and l_orient == r_orient:

			for sr_var_temp in SRVarHash:
				if (SRVarHash[sr_var_temp].typeO == "INS" or SRVarHash[sr_var_temp].typeO == "INS_I"):
					
					if SRVarHash[sr_var_temp].bp3 == -1 and (SRVarHash[temp].bp2-SR_SLOP/2 <= temp.bp <= SRVarHash[sr_var_temp].bp2 +SR_SLOP/2):

						if SRVarHash[sr_var_temp].typeO == "INS":
							SRVarHash[sr_var_temp].typeO == "INS_C"
						elif SRVarHash[sr_var_temp].typeO == "INS_I":
							SRVarHash[sr_var_temp].typeO == "INS_C_I" 
						
						SRVarHash[temp].bp3 = other_bp

					elif SRVarHash[sr_var_temp].bp3 == -1 and (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[sr_var_temp].bp2 +SR_SLOP/2):
					
						if SRVarHash[sr_var_temp].typeO == "INS":
                                                	SRVarHash[sr_var_temp].typeO == "INS_C"
                                        	elif SRVarHash[sr_var_temp].typeO == "INS_I":
                                                	SRVarHash[sr_var_temp].typeO == "INS_C_I"

                                       		 SRVarHash[temp].bp3 = temp.bp

					elif SRVarHash[temp].bp2-SR_SLOP/2 <= temp.bp <= SRVarHash[sr_var_temp].bp2 +SR_SLOP/2 and SRVarHash[temp].bp3-SR_SLOP/2 <= other_bp <= SRVarHash[sr_var_temp].bp3 +SR_SLOP/2:
						if SRVarHash[sr_var_temp].typeO == "INS":
                                                        SRVarHash[sr_var_temp].typeO == "INS_C"
                                                elif SRVarHash[sr_var_temp].typeO == "INS_I":
                                                        SRVarHash[sr_var_temp].typeO == "INS_C_I"
					
					elif SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[sr_var_temp].bp2 +SR_SLOP/2 and SRVarHash[temp].bp3-SR_SLOP/2 <= temp.bp <= SRVarHash[sr_var_temp].bp3 +SR_SLOP/2:
					
						if SRVarHash[sr_var_temp].typeO == "INS":
                                                        SRVarHash[sr_var_temp].typeO == "INS_C"
                                                elif SRVarHash[sr_var_temp].typeO == "INS_I":
                                                        SRVarHash[sr_var_temp].typeO == "INS_C_I"
		#If unclaimed, add to SR hash
		#$add claimed variable
		else:
			# hash in "left" bp to look-up table
			temp.bp = sr_bp1
			temp.tid = sr_bp1_tid
			temp.tid_2 = sr_bp2_tid

			newVariant = newSRVar()
			newVariant.bp2 = sr_bp2
			newVariant.support.append(SRFrag)
			newVariant.l_orient = l_orient
			newVariant.r_orient = r_orient

			#print swap, "is swap"
			if temp.tid == temp.tid_2 and not swap:

				newVariant.swapped = 0

				if l_orient == r_orient:
					newVariant.typeO = "DEL_INS"
				else:
					newVariant.typeO = "INS_I"
					
				#print "DELINS"
			elif temp.tid == temp.tid_2 and swap:

				newVariant.typeO = "TD_I"
				newVariant.swapped = 1
				#print "TD_I"

			elif temp.tid != temp.tid_2:
				newVariant.typeO = "INS"
                                newVariant.swapped = 1
				if l_orient != r_orient:
					newVariant.typeO = "INS_I"

			elif temp.tid == temp.tid_2 and l_orient != r_orient:
				newVariant.typeO = "INV_POSS" 
				newVariant.swapped = 0
					 
				#print "INV_P"
			list1 = range(int(sr_bp1-SR_SLOP/2),int(sr_bp1+SR_SLOP/2) + 1)
                        list2 = range(int(sr_bp2-SR_SLOP/2),int(sr_bp2+SR_SLOP/2) + 1)

                        for x in list1:
				newVariant2 = newSRVar()
				newVariant2.bp2 = newVariant.bp2
				newVariant2.support.append(SRFrag)
				newVariant2.l_orient = newVariant.l_orient
				newVariant2.r_orient = newVariant.r_orient
				newVariant2.typeO = newVariant.typeO
				newVariant2.swapped = newVariant.swapped
				#temp2 = copy.deepcopy(temp)
				temp2 = Variant()
				temp2.tid = temp.tid
				temp2.tid_2 = temp.tid_2
                                temp2.bp = x
				newSRList.append(temp2)
				if newSRList[-1] not in SRVarHash:
	                                SRVarHash[newSRList[-1]] = newVariant2

                                #print SRVarHash[temp].typeO

                        newVariant.bp2 = sr_bp1
			swap_f(newVariant.l_orient, newVariant.r_orient)
			swap_f(temp.tid, temp.tid_2)

                        for x in list2:
				newVariant2 = newSRVar()
                                newVariant2.bp2 = newVariant.bp2
                                newVariant2.support.append(SRFrag)
                                newVariant2.l_orient = newVariant.l_orient
                                newVariant2.r_orient = newVariant.r_orient
                                newVariant2.typeO = newVariant.typeO
                                newVariant2.swapped = newVariant.swapped
                                #temp2 = copy.deepcopy(temp)
				temp2 = Variant()
                                temp2.tid = temp.tid
                                temp2.tid_2 = temp.tid_2
                                temp2.bp = x
				newSRList.append(temp2)
				if newSRList[-1] not in SRVarHash:
	                                SRVarHash[newSRList[-1]] = newVariant2
		
			#print len(newSRList), newSRList[0],newSRList[-1]
			#print temp, temp in SRVarHash
			#print "Case", temp, type(temp.bp), type(temp.tid), type(temp.tid_2), temp in SRVarHash, SRVarHash[temp], newVariant2.typeO, newVariant.typeO, newVariant2, newVariant	

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

    #post-process new SR variants to merge hashed neighbors and hash pairs
    

    #write new SR variants to file	
    k = 0
    for item in SRVarHash:
	k+=1
	#print item, SRVarHash[item],SRVarHash[item].typeO, SRVarHash[item].count
	if SRVarHash[item].count > 0:
			
		if SRVarHash[item].typeO[:3] == "INS":
			if SRVarHash[item].bp3 > SRVarHash[item].bp2:
				swapper = SRVarHash[item].bp3
				SRVarHash[item].bp3 = SRVarHash[item].bp2
				SRVarHash[item].bp2 = swapper
		elif SRVarHash[item].typeO == "DEL_INS":
			SRVarHash[item].typeO = "DEL"
		#print SRVarHash[item].typeO
		
		bp_temp = item.bp
		bp_temp_mate = SRVarHash[item].bp2
		chosen_var = Variant()
		neighborSupport = []

		#check neighbor hash pairs and transfer to one with change in type (bona fide), else pick 1 randomly
		if SRVarHash[item].write == -1:
		   for x in range(int(bp_temp-SR_SLOP-1), int(bp_temp + SR_SLOP+1)):
			temp = Variant()
        	        temp.tid = item.tid
               		temp.tid_2 = item.tid_2
			temp.bp = x

			#if in same variant group or "symmetry group" 
			if temp in SRVarHash and SRVarHash[temp].support[0] == SRVarHash[item].support[0]:
				
				if SRVarHash[temp].n_changes > 0 and chosen_var.bp == -1:
					chosen_var = temp
					SRVarHash[temp].write = 1
				else:
					SRVarHash[temp].write = 0
					for ns in SRVarHash[temp].support[1:]:
						neighborSupport.append(ns)
			
		   #if none, check mate hash pairs and transfer to one with most changes in type (bona fide), else pick 1 randomly
			
		   for x in range(int(bp_temp_mate-SR_SLOP -1), int(bp_temp_mate + SR_SLOP +1)):
			temp = Variant()
                        temp.tid = item.tid_2
                        temp.tid_2 = item.tid
                        temp.bp = x

                        if temp in SRVarHash and SRVarHash[temp].support[0] == SRVarHash[item].support[0]: 
				if SRVarHash[temp].n_changes > 0 and chosen_var.bp == -1:
					chosen_var = temp
					SRVarHash[temp].write = 1
				else:
					SRVarHash[temp].write = 0
					for ns in SRVarHash[temp].support[1:]:
						neighborSupport.append(ns)        

		   #if all neighbors unchanged, just pick 1 among them
		   if chosen_var.bp == -1:
			chosen_var = item
			SRVarHash[chosen_var].write = 1	
		   #print "Chosen:", chosen_var, item		
		   SRVarHash[chosen_var].support = SRVarHash[chosen_var].support + list(set(neighborSupport) - set(SRVarHash[chosen_var].support))
		   SRVarHash[chosen_var].count = len(SRVarHash[chosen_var].support)
	
		if SRVarHash[item].write == 1 and SRVarHash[item].count >= MIN_CS_SR:
			fp3.write("%s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(k+varNum+1,SRVarHash[item].typeO, item.tid, item.bp-SR_SLOP/2, item.bp+SR_SLOP/2, item.tid_2, SRVarHash[item].bp2-SR_SLOP/2, SRVarHash[item].bp2+SR_SLOP/2, SRVarHash[item].bp3tid, SRVarHash[item].bp3, SRVarHash[item].bp3, "SR", SRVarHash[item].count))
			fp4.write("%s" %(k+varNum+1))
          		for elem in SRVarHash[item].support:
				fp4.write(" %s" %elem)
			fp4.write("\n") 
                                                                                     
    bamfile.close()
    
