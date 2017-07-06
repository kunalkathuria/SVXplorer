import sys

def overlap(locA_1, locA_2, locB_1, locB_2):

	if locA_1 > locB_1 and locA_1 < locB_2:
		if locA_2 < locB_2:
			region = [locA_1, locA_2]
		else:
			region = [locA_1, locB_2]
		return 1, region

	if locA_2 > locB_1 and locA_2 < locB_2:
		if locA_1 < locB_1:
			region = [locB_1, locA_2]
		else:
			region = [locA_1, locA_2]
		return 1, region

	if locA_1 < locB_1 and locA_2 > locB_2:
		region = [locB_1, locB_2]
		return 1, region
	
	return 0, [0,0]

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
		#print "Variant", 1+g
		count=0
		Regions = []
		varTagged = []
		conflict = 0
		for k,item in enumerate(varList):
		  if k > g: 
			typeV2 = item.split()[-1]
			TID_2 = item.split()[0]
			bpB_1 = item.split()[1]
			bpB_2 = item.split()[5]
		# set type etc.bp's are extrema points.convert to int.bp1_1 etc
			if (typeV1 == "DEL" or typeV1 == "TD" or typeV1 == "INV") and (typeV2 == "DEL" or typeV2 == "TD"or typeV2 == "INV"):

				if TID_1 == TID_2:
				   [overlaps, reg_overlap] = overlap(int(bpA_1),int(bpA_2), int(bpB_1), int(bpB_2))
				   if overlaps:
				      for h,region_c in enumerate(Regions):
					if overlap(reg_overlap[0], reg_overlap[1], region_c[0], region_c[1]):
						
						#print "Conflict:", 1+g, 1+k, 1+varTagged[h]
						conflict_code = map_code(typeV1) + map_code(typeV2) + map_code(varList[varTagged[h]].split()[-1])
						print conflict_code
						print elem
						print varList[varTagged[h]]
						print item, "\n"
						conflict = 1
						n_conflict_v+=1
						break # just count how many variants are involved in conflict
				      Regions.append(reg_overlap)
				      varTagged.append(k)
			if conflict:
				break
	print n_conflict_v, "total conflicts."	
