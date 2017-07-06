# checks for genuine conflict
# put deletions first in merged bedpe variant file
import sys

BUFFER_L = int(sys.argv[2])

class ConflictCluster(object):

	def __init__(self):
		self.SVtype1 = -1
		self.SVtype2 = -1
		self.SVtype3 = -1
		self.confl_code = -1
		self.limit_l = -1
		self.limit_r = -1
		self.clusterNums = [-1,-1,-1]
	def __str__(self):
		self.clusterNums = sorted(self.clusterNums)
		return "%s %s %s %s" %(self.confl_code, self.clusterNums[0], self.clusterNums[1], self.clusterNums[2])	

def isConflicting(bpA, bpB, typeSV, pc):

	conflict = 0
	
	if pc.SVtype1 == "DEL" and pc.SVtype2 == "DEL":
		
		if typeSV == "DEL" or typeSV == "TD" or typeSV == "INV":

			if bpA > pc.limit_l + BUFFER_L and bpB < pc.limit_r - BUFFER_L:
				
				conflict = 1

		# INS condn here will depend on copy BPs

	#elif pc.SVtype1 == "INV" and pc.SVtype2 == "INV":
		
		#if typeSV == "TD" or typeSV == "INV":
			
			#conflict = 1	
	
	if conflict:
		pc.clusterNums[2] = k
      		pc.SVtype3 = typeV2
        	pc.confl_code = map_code(pc.SVtype3) + map_code(pc.SVtype1) + map_code(pc.SVtype2)	

def overlap(locA_1, locA_2, locB_1, locB_2, SV1, SV2, SVn1, SVn2):

	temp = ConflictCluster()

	if (SV1 == "DEL" and SV2 == "DEL"): #or (SV1 =="INV" and SV2 == "INV"):
	
		if locA_1 > locB_1 and locA_2 < locB_2:
			temp.limit_l = locB_1
			temp.limit_r = locB_2
			temp.SVtype1 = SV1
			temp.SVtype2 = SV2	
			temp.clusterNums[0] = SVn1
			temp.clusterNums[1] = SVn2	
			return 1, temp

		if locA_1 < locB_1 and locA_2 > locB_2:
			temp.limit_l = locB_1
			temp.limit_r = locB_2
			temp.SVtype1 = SV1
			temp.SVtype2 = SV2
			temp.clusterNums[0] = SVn1
			temp.clusterNums[1] = SVn2
			return 1, temp
	
	return 0, 0	

def map_code(sv_type):

        if sv_type == "DEL":
                return 1
        elif sv_type == "TD":
                return 10
        elif sv_type == "INV":
                return 100
        elif sv_type[0:3] == "INS":
                return 1000
        else:
                return 0

if __name__=="__main__":

	varList = []
	f = open(sys.argv[1], "r") # all bedpe's merged into 1 file
	n_conflict_v=0

	for line in f:
		varList.append(line)

	for g,elem in enumerate(varList):
		typeV1 = elem.split()[-1]
		bpA_1 = elem.split()[1]
		bpA_2 = elem.split()[5]
		TID_1 = elem.split()[0]
		pot_conflicts = []
		conflict = 0
		print g
		for k,item in enumerate(varList):
		   
			typeV2 = item.split()[-1]
			TID_2 = item.split()[0]
			bpB_1 = item.split()[1]
			bpB_2 = item.split()[5]

			#$include INS with bp2 and bp3 as conflict candidate, and check if all cluster orderings accounted for
			if (typeV1 == "DEL" or typeV1 == "TD" or typeV1 == "INV") and (typeV2 == "DEL" or typeV2 == "TD"or typeV2 == "INV"):

				if TID_1 == TID_2:
				   [overlaps, confl_cl] = overlap(int(bpA_1),int(bpA_2), int(bpB_1), int(bpB_2), typeV1, typeV2, g, k)
				   if overlaps:
				      #print confl_cl
				      for h, pc in enumerate(pot_conflicts):
					if isConflicting(int(bpB_1), int(bpB_2), typeV2, pc):
						
						print "Conflict:", pc						
						print elem
						print item, "\n"
						conflict = 1
						n_conflict_v+=1
						#break # just count how many variants are involved in conflict
				      pot_conflicts.append(confl_cl)
			#if conflict:
				#break
	print n_conflict_v, "total conflicts."	
