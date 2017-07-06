# 12/16/16
# Set Cover Algorithms to Minimize number of SVs that account for all PEMs pointing to them

import math
import os
import sys
from collections import Counter

FILE = "../results/text/ClassifiedVariantMap.txt"
CALC_LIMIT = 500000

def DisjointAlg(fragmentList, nSetsR):

    varNum = []
    Sets = []
    #results = [0]*(max(fragmentList) + 1)

    # fragment numbering starts with 1

    #claimed = [0]*(1+max(fragmentList)) # change to n unique frags if saves time

    fp=open(FILE,"r")
    fp3=open("../results/text/VariantProfile_match_count_test.txt","w")

    # Obtain disjointness count
    # Assume any given fragment only appears once in set list
    for k,line in enumerate(fp):

	print "Set", k, "written."
        currentSet = map(int, line.split())
        varNum.append(currentSet[0])
	currentSet = currentSet[1:] 	
	Sets.append(set(currentSet))

    matchCount0 = []
    matchCount1 = []
    matchCount2 = []

    for g,elem in enumerate(Sets):
	print "Set", g, "analyzed."
	
	if g == CALC_LIMIT:
		break

	for h,item in enumerate(Sets):
		if h == 0:
                        matchCount0.append(0)
			matchCount1.append(0)
			matchCount2.append(0)
		if g != h:
			# if identical
			if elem == item:
				matchCount0[g]+= 1

			# if proper subset	
			elif not elem.isdisjoint(item) and len(elem - item) == 0 and len(item-elem) > 0:
				matchCount1[g]+=1

			# if superset
			elif not elem.isdisjoint(item) and len(elem-item) > 0 and len(item-elem)==0:
				matchCount2[g]+=1
		
				 	

    for g,item in enumerate(Sets):

	if g < CALC_LIMIT:
		fp3.write("%s %s %s %s\n" %(varNum[g], matchCount0[g], matchCount1[g], matchCount2[g]))

    fp3.close()
    return -1

    #potSVs = [item for item in disjointness if item > disjThresh]        
    
def ReadFile(filename):

    allFrags = []
    counter=0
    f=open(filename, 'r')
    for line in f:

        parsed = map(int, line.split())
        counter = counter + 1
        for x in parsed[1:]:
            allFrags.append(x)

    f.close()
    return allFrags, counter
      
    
if __name__ == "__main__":  

   [allFrags, nSetsR] = ReadFile(FILE)
   
   nSV =  DisjointAlg(allFrags, nSetsR)
   
   f= open("../results/text/NSVs_2.txt","w")
   
   f.write("%s\n" %nSV)
   f.close()

   #WriteMappings(DisjFM, "FragMapD.txt")
 
   
