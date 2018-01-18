#!/usr/bin/env python
# Kunal Kathuria 1/18
# Add Split Reads to existing PE clusters

import pysam
import sys
from random import randint
import copy

WORK_DIR=sys.argv[10]
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
SR_INS_THRESH = int(sys.argv[8])
SR_support_thresh = int(sys.argv[9])

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
		self.isoriginal = -1
		self.instoinv = -1
	def __str__(self):

 	       return "%s %s %s %s %s" %(self.bp2, self.bp3, self.count, self.isoriginal, self.typeO)

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

def swap_f(a):
	temp = a.l_orient
	a.l_orient = a.r_orient
	a.r_orient = temp

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
	elif SV_type == "INS_I":
		return 4
	else:
		return -1

def formHash(fp1, iobjects, varHash):

	print "Forming PE variant hash table..."

        for line in fp1:
                line_s = line.split()
		bounds = Bound()		
		bounds.num = int(line_s[0])			
		bounds.typeO= mapSVtoNum(line_s[1])

		if bounds.typeO == -1:
			continue
		if bounds.typeO == 3 or bounds.typeO == 4:

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

    fp1 = open(WORK_DIR+"/allVariants.pe.txt","r")
    fp2 = open(WORK_DIR+"/variantMap.pe.txt","r")
    fp3 = open(WORK_DIR+"/allVariants.pe_sr.txt","w")
    fp4 = open(WORK_DIR+"/variantMap.pe_sr.txt","w")
    fp5 = open(WORK_DIR+"/updatedBPs.txt","w")
    varHash = {}
    SRVarHash = {}
    # ensure hash objects immutable
    immutable_objects = []
    formHash(fp1, immutable_objects, varHash)

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

		#print sr1, sr2
                SRFrag-=1

                if SRFrag % (-1000) == 0:

                    print SRFrag

                sr_bp1 = sr1.reference_start
                sr_bp2 = sr2.reference_start
                sr_bp1_tid = sr1.reference_name
		sr_bp2_tid = sr2.reference_name

		ignoreList = []
		if ignoreChr != "none":
			fo=open(ignoreChr, "r")
			for entry in fo:
				ignoreList.append(entry[:-1])

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

		#QAS below refers to the alignment position of the aligned portion of the read relative to the whole read depicted in the reference
		#swap criterion works well for 75% of inversion reads but even if their swap parameter is incorrect, no logical error
		#occurs -- since other reads are simply used to construct the inversions while the rest are simply neglected
                if sr_bp1_tid == sr_bp2_tid and minsr.query_alignment_start > maxsr.query_alignment_start:

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

                # this includes case of query start positions (QAS) of both reads being 0, which happens in 50% of inversion split reads        
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
		#print "Please check if split-read file is name-sorted."
		continue

	#this is not perfect for all inverted types, above is almost perfect: see explanation above
	#if swap == -1 and ((minsr.cigarstring.find("S") !=-1 and minsr.cigarstring.find("S") < minsr.cigarstring.find("M")) or  (minsr.cigarstring.find("H") != -1 and minsr.cigarstring.find("H") < minsr.cigarstring.find("M"))):
		#swap = 1
	#elif swap == -1:
		#swap = 0

	#if minsr.query_alignment_start == maxsr.query_alignment_start:
		#print swap, minsr.query_alignment_start, minsr.query_alignment_end, maxsr.query_alignment_start, maxsr.query_alignment_end, minsr.cigarstring, maxsr.cigarstring

	match = 0
	temp = Variant()
	temp.bp = sr_bp1
	temp.tid = sr_bp1_tid
	#print "Formed sr:", temp

	if sr_bp1_tid != sr_bp2_tid:	                                                      

		temp.tid_2 = sr_bp2_tid

	if temp in varHash and (varHash[temp].bp2_1 < sr_bp2 < varHash[temp].bp2_2 or varHash[temp].bp3_1 < sr_bp2 < varHash[temp].bp3_2):

		varNum = varHash[temp].num
		varType = varHash[temp].typeO
		#print "Found", varNum, varType, varHash[temp].bp2_1, varHash[temp].bp2_2, sr_bp2, varHash[temp].bp3_1, varHash[temp].bp3_2

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

		if temp in varHash and (varHash[temp].bp2_1 < sr_bp1 < varHash[temp].bp2_2 or varHash[temp].bp3_1 < sr_bp1 < varHash[temp].bp3_2):
			varNum = varHash[temp].num
			varType = varHash[temp].typeO
		        #print "Found", varNum, varType, varHash[temp].bp2_1, varHash[temp].bp2_2, sr_bp2, varHash[temp].bp3_1, varHash[temp].bp3_2
	
			if varNum not in matchingFrags:
				newbp = SRVar()
				newbp.bp1 = sr_bp2
				newbp.bp2 = sr_bp1

	#check DEL
	if varType == 0 and swap==0 and minsr.is_reverse == maxsr.is_reverse:

		match = 1
	#elif varType == 0:
		#print sr_bp1_tid, sr_bp1, sr_bp2_tid, sr_bp2, swap, minsr.is_reverse, maxsr.is_reverse

	#check TD
	elif varType == 1 and swap==1 and minsr.is_reverse == maxsr.is_reverse:

		match = 1

	#check INV
	elif varType == 2 and swap==0 and minsr.is_reverse != maxsr.is_reverse:
	
		match = 1

	#check INS
	elif (varType == 3 and minsr.is_reverse == maxsr.is_reverse) or (varType == 4 and minsr.is_reverse != maxsr.is_reverse):

		match = 1
		print "Match with INS", varNum
		
		# Set new bp3 of insertion if unset
		if varNum in matchingFrags and matchingFrags[varNum][0].bp3 == -1:

			temp2 = matchingFrags[varNum][0]
			if abs(sr_bp1 - temp2.bp1) > 2*SR_SLOP and abs(sr_bp1 - temp2.bp2) > 2*SR_SLOP and varHash[temp].bp3_1 < sr_bp1 < varHash[temp].bp3_2:
			
				matchingFrags[varNum][0].bp3 = sr_bp1

			elif abs(sr_bp2 - temp2.bp1) > 2*SR_SLOP and abs(sr_bp2 - temp2.bp2) > 2*SR_SLOP and varHash[temp].bp3_1 < sr_bp2 < varHash[temp].bp3_2:
		
				matchingFrags[varNum][0].bp3 = sr_bp2

			# insertion bp2 should be < bp3
                        if matchingFrags[varNum][0].bp3 != -1 and matchingFrags[varNum][0].bp3 < matchingFrags[varNum][0].bp2:
			
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

			#Do not look for inverted insertions due to ambiguity of swap parameter
			#Also, INS_I would be off by ~.5 RDL or less due to ambiguity of swapped alignment: whether ref start or end alignment point is breakpoint
                        if SRVarHash[temp].typeO == "DEL_INS" or SRVarHash[temp].typeO == "DEL":
				
				if l_orient == r_orient and swap==1 and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2) and not ( (l_orient != SRVarHash[temp].l_orient and SRVarHash[temp].r_orient == r_orient) or (l_orient == SRVarHash[temp].l_orient and SRVarHash[temp].r_orient != r_orient) ):
					SRVarHash[temp].typeO = "INS" 
                       	
					SRVarHash[temp].bp3 = other_bp
					SRVarHash[temp].bp3tid = other_bp_tid

					SRVarHash[temp].count+=1
	                                SRVarHash[temp].support.append(SRFrag)

					SRVarHash[temp].n_changes+=1
					claimed = 1
					print "DEL_INS"

				# cannot update type based on this info	
				elif swap==0 and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient == r_orient == SRVarHash[temp].l_orient:

					SRVarHash[temp].count+=1
	                                SRVarHash[temp].support.append(SRFrag)

				elif swap==0 and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient == r_orient and l_orient != SRVarHash[temp].l_orient:
					#If 2 bps supported by 2 mutually reverse mappings, count as DEL. Following gave many FPs -- they may be partially formed INSs
					#if SRVarHash[temp].typeO == "DEL_INS":
						#SRVarHash[temp].typeO = "DEL"
						#SRVarHash[temp].n_changes+=1
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
				
				#INS_C
				#elif not swap and l_orient == r_orient and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2):
        		#TD_I
                        elif SRVarHash[temp].typeO == "TD_I":
				#cannot update type here
				if swap==1 and l_orient == r_orient and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2:

	                                SRVarHash[temp].count+=1
        	                        SRVarHash[temp].support.append(SRFrag)
					#print "TD 1", SRVarHash[temp].count

				elif l_orient == r_orient and swap==0 and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or SRVarHash[temp].bp2 < temp.bp < other_bp or other_bp < temp.bp < SRVarHash[temp].bp2) and not ( (l_orient != SRVarHash[temp].l_orient and SRVarHash[temp].r_orient == r_orient) or (l_orient == SRVarHash[temp].l_orient and SRVarHash[temp].r_orient != r_orient) ):

					print "TD 2"
					SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
					SRVarHash[temp].typeO = "INS"
					SRVarHash[temp].bp3 = other_bp
					SRVarHash[temp].bp3tid = other_bp_tid

					SRVarHash[temp].n_changes+=1

                        elif SRVarHash[temp].typeO[:3] == "INV" and l_orient != r_orient:
				#cannot update type if INV_POSS currently, as only one half of inversion reported thus far
				if SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient == SRVarHash[temp].l_orient:
					SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)

				elif SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient != SRVarHash[temp].l_orient:	
					#$banks on more inversions existing than inverted insertions -- statistically no use of SWAP here seems to work
					if SRVarHash[temp].typeO == "INV_POSS":
						SRVarHash[temp].typeO = "INV"
						SRVarHash[temp].n_changes+=1
                                	SRVarHash[temp].count+=1
                                	SRVarHash[temp].support.append(SRFrag)

				elif not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or SRVarHash[temp].bp2 < temp.bp < other_bp or other_bp < temp.bp < SRVarHash[temp].bp2) and not ((l_orient != SRVarHash[temp].l_orient and SRVarHash[temp].r_orient == r_orient) or (l_orient == SRVarHash[temp].l_orient and SRVarHash[temp].r_orient != r_orient)):
					SRVarHash[temp].typeO = "INS_I"
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
					SRVarHash[temp].bp3 = other_bp
					SRVarHash[temp].bp3tid = other_bp_tid
					SRVarHash[temp].n_changes+=1
					print "INV INS"

			elif SRVarHash[temp].typeO == "INS":
				print "INS"
				if l_orient == r_orient and (SRVarHash[temp].bp3 == -1 or SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or SRVarHash[temp].bp3-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp3+SR_SLOP/2):
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)
					if SRVarHash[temp].bp3 == -1 and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or SRVarHash[temp].bp2 < temp.bp < other_bp or other_bp < temp.bp < SRVarHash[temp].bp2):
                                                SRVarHash[temp].bp3 = other_bp
						SRVarHash[temp].bp3tid = other_bp_tid
						SRVarHash[temp].n_changes+=1
					print SRVarHash[temp]

			#$ if 3rd bp absent write as inversion (swap can be regular INV due to swap error)
			elif SRVarHash[temp].typeO == "INS_I":
		
				print "INS_I"
                                if l_orient != r_orient and (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or SRVarHash[temp].bp3 == -1 or SRVarHash[temp].bp3-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp3+SR_SLOP/2):
                                        SRVarHash[temp].count+=1
                                        SRVarHash[temp].support.append(SRFrag)

					if SRVarHash[temp].bp3 == -1 and not (SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 or SRVarHash[temp].bp2 < temp.bp < other_bp or other_bp < temp.bp < SRVarHash[temp].bp2):
						SRVarHash[temp].bp3 = other_bp
						SRVarHash[temp].bp3tid = other_bp_tid
						SRVarHash[temp].n_changes+=1
					elif temp.tid == temp.tid_2 and SRVarHash[temp].bp2-SR_SLOP/2 <= other_bp <= SRVarHash[temp].bp2+SR_SLOP/2 and l_orient != SRVarHash[temp].l_orient:
						#if 3rd bp unset till end, then call it an inversion if this condition is fulfilled
						SRVarHash[temp].instoinv=1
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
			if temp.tid == temp.tid_2 and swap==0:

				newVariant.swapped = 0

				if l_orient == r_orient:
					newVariant.typeO = "DEL_INS"
				else:
					newVariant.typeO = "INV_POSS"
					
				#print "DELINS"
			elif temp.tid == temp.tid_2 and swap==1:

				#print "TD formed", temp, sr_bp2
				newVariant.typeO = "TD_I"
				newVariant.swapped = 1
				#$handle this case in INS_I matches
				if l_orient != r_orient:
					newVariant.typeO = "INS_I"
				#print "TD_I"

			elif temp.tid != temp.tid_2:
				newVariant.typeO = "INS"
                                newVariant.swapped = 0
				if l_orient != r_orient:
					newVariant.typeO = "INS_I"

			newVariant.isoriginal = 1
			if temp not in SRVarHash:
				SRVarHash[temp] = newVariant
				#print "Set original", SRVarHash[temp]

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

                        for x in list2:
				newVariant2 = newSRVar()
                                newVariant2.bp2 = sr_bp1
                                newVariant2.support.append(SRFrag)
                                newVariant2.l_orient = newVariant.r_orient
                                newVariant2.r_orient = newVariant.l_orient
                                newVariant2.typeO = newVariant.typeO
                                newVariant2.swapped = newVariant.swapped
                                #temp2 = copy.deepcopy(temp)
				temp2 = Variant()
                                temp2.tid = temp.tid_2
                                temp2.tid_2 = temp.tid
                                temp2.bp = x
				newSRList.append(temp2)
				if newSRList[-1] not in SRVarHash:
	                                SRVarHash[newSRList[-1]] = newVariant2
	

    # print new variant map
    fp1.seek(0)


    for line3 in fp2:
        
        line3_split = line3.split()
        varNum = int(line3_split[0])
        #print "In loop", varNum
 
        for line in fp1:
            
            line_split = line.split()
	    #if varNum in matchingFrags and line_split[1] == "INS" or line_split[1] == "INS_C" or line_split[1] == "INS_I" or line_split[1] == "INS_C_I":
                #print "Before:", line_split[3], line_split[4], line_split[6], line_split[7], line_split[9], line_split[10]

            if varNum in matchingFrags:
                line_split[11] = line_split[11] + "_SR"
		line_split[3] = str(matchingFrags[varNum][0].bp1)
		line_split[4] = str(matchingFrags[varNum][0].bp1 + 1)
		if len(matchingFrags[varNum]) > SR_support_thresh and matchingFrags[varNum][0].bp3 == -1:
			if int(line_split[6]) < matchingFrags[varNum][0].bp2 < int(line_split[7]):
				line_split[6] = str(matchingFrags[varNum][0].bp2)
				line_split[7] = str(matchingFrags[varNum][0].bp2 + 1)
			elif int(line_split[9]) < matchingFrags[varNum][0].bp2 < int(line_split[10]):
				line_split[9] = str(matchingFrags[varNum][0].bp2)
                                line_split[10] = str(matchingFrags[varNum][0].bp2 + 1)

			if (line_split[1] == "DEL" or line_split[1] == "TD" or line_split[1] == "INV") and int(line_split[3]) > int(line_split[7]):
				temp = line_split[3]
				temp2 = line_split[4]
				line_split[3] = line_split[6]
				line_split[4] = line_split[7] 
				line_split[6] = temp
				line_split[7] = temp2
		# full insertion matches
		elif len(matchingFrags[varNum]) > SR_support_thresh:
			line_split[6] = str(matchingFrags[varNum][0].bp2)
                        line_split[7] = str(matchingFrags[varNum][0].bp2 + 1)
			line_split[9] = str(matchingFrags[varNum][0].bp3)
			line_split[10] = str(matchingFrags[varNum][0].bp3 + 1)
           
	    #if varNum in matchingFrags and line_split[1] == "INS" or line_split[1] == "INS_C" or line_split[1] == "INS_I" or line_split[1] == "INS_C_I":
		#print "After:", line_split[3], line_split[4], line_split[6], line_split[7], line_split[9], line_split[10] 

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
			
		#This leads to low precision
		#if SRVarHash[item].typeO == "DEL_INS":
			#SRVarHash[item].typeO = "DEL"
		#print item, SRVarHash[item], SRVarHash[item].support
		
		bp_temp = item.bp
		bp_temp_mate = SRVarHash[item].bp2
		chosen_var = Variant()
		neighborSupport = []
		orig_sv = item
		
		#if SRVarHash[item].isoriginal == 1:
			#print "Original:", item, SRVarHash[item], SRVarHash[item].typeO

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
					#print "Write 1", chosen_var.bp
					chosen_var = temp
					SRVarHash[temp].write = 1
				else:
					#print "Write 0:", temp, SRVarHash[temp] 
					SRVarHash[temp].write = 0
					for ns in SRVarHash[temp].support[1:]:
						neighborSupport.append(ns)
				if SRVarHash[temp].isoriginal == 1:
                                        orig_sv = temp
					#print "Marked original", SRVarHash[temp]
		   #if none, check mate hash pairs and transfer to one with most changes in type (bona fide), else pick 1 randomly
			
		   for x in range(int(bp_temp_mate-SR_SLOP -1), int(bp_temp_mate + SR_SLOP +1)):
			#print "Second loop:", x
			temp = Variant()
                        temp.tid = item.tid_2
                        temp.tid_2 = item.tid
                        temp.bp = x

                        if temp in SRVarHash and SRVarHash[temp].support[0] == SRVarHash[item].support[0]: 
				if SRVarHash[temp].n_changes > 0 and chosen_var.bp == -1:
					#print "Write 1 loop 2:", chosen_var.bp
					chosen_var = temp
					SRVarHash[temp].write = 1
				else:
					#print "Write 0 loop 2:", temp, SRVarHash[temp]
					SRVarHash[temp].write = 0
					for ns in SRVarHash[temp].support[1:]:
						neighborSupport.append(ns)        
				if SRVarHash[temp].isoriginal == 1:
					orig_sv = temp

		   #if all neighbors unchanged, pick the original one among all neighbors
		   if chosen_var.bp == -1:
			#print "Choosing", orig_sv
			chosen_var = orig_sv
			SRVarHash[chosen_var].write = 1	
		   #print "Chosen:", chosen_var, SRVarHash[chosen_var]
		   SRVarHash[chosen_var].support = SRVarHash[chosen_var].support + list(set(neighborSupport) - set(SRVarHash[chosen_var].support))
		   SRVarHash[chosen_var].count = len(SRVarHash[chosen_var].support)
	
		if SRVarHash[item].write == 1 and SRVarHash[item].count >= MIN_CS_SR:

			if SRVarHash[item].typeO == "INS_I" and SRVarHash[item].instoinv == 1 and SRVarHash[item].bp3 == -1:

				#$can make this "INV" from "INV_POSS" if wish to be liberal, check FP increase
				SRVarHash[item].typeO == "INV"

			elif SRVarHash[item].typeO == "INS_I" and SRVarHash[item].bp3 == -1:

                                #$can make this "INV" from "INV_POSS" if wish to be liberal, check FP increase
                                SRVarHash[item].typeO == "INV_POSS"
			

			if (SRVarHash[item].typeO == "INS_I" or SRVarHash[item].typeO == "INS") and item.tid_2 == SRVarHash[item].bp3tid and abs(SRVarHash[item].bp2 - SRVarHash[item].bp3) < SR_INS_THRESH:
				SRVarHash[item].typeO = "INS_POSS"

			elif (SRVarHash[item].typeO == "INS_I" or SRVarHash[item].typeO == "INS") and item.tid_2 == SRVarHash[item].bp3tid and SRVarHash[item].bp2 > SRVarHash[item].bp3:

				swapper = SRVarHash[item].bp2
				SRVarHash[item].bp2 = SRVarHash[item].bp3
				SRVarHash[item].bp3 = swapper

			if (SRVarHash[item].typeO == "DEL_INS" or SRVarHash[item].typeO == "TD_I" or SRVarHash[item].typeO == "INV" or SRVarHash[item].typeO == "INV_POSS") and item.bp > SRVarHash[item].bp2:
				fp3.write("%s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(k+varNum+1,SRVarHash[item].typeO, item.tid_2, SRVarHash[item].bp2, SRVarHash[item].bp2+1, item.tid, item.bp, item.bp + 1, SRVarHash[item].bp3tid, SRVarHash[item].bp3, SRVarHash[item].bp3 + 1, "SR", SRVarHash[item].count))

			else:
				fp3.write("%s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(k+varNum+1,SRVarHash[item].typeO, item.tid, item.bp, item.bp + 1, item.tid_2, SRVarHash[item].bp2, SRVarHash[item].bp2 + 1, SRVarHash[item].bp3tid, SRVarHash[item].bp3, SRVarHash[item].bp3 + 1, "SR", SRVarHash[item].count))

			fp4.write("%s" %(k+varNum+1))
          		for elem in SRVarHash[item].support:
				fp4.write(" %s" %elem)
			fp4.write("\n") 
                                                                                     
    bamfile.close()
    
