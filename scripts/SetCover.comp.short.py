# 12/16/16
# Set Cover Algorithms to Minimize number of SVs that account for all PEMs pointing to them

import math
import os
import sys
from collections import Counter

FILE = "../results/text/VariantMap.txt"

def DisjointAlg(fragmentList, nSetsR):

    varNum = []
    Sets = []
    #results = [0]*(max(fragmentList) + 1)

    # fragment numbering starts with 1

    #claimed = [0]*(1+max(fragmentList)) # change to n unique frags if saves time

    fp=open(FILE,"r")
    fp3=open("../results/text/VariantProfile.txt","w")

    # Obtain disjointness count
    # Assume any given fragment only appears once in set list
    for k,line in enumerate(fp):

	print "Set", k, "written."
        currentSet = map(int, line.split())
        varNum.append(currentSet[0])
	currentSet = currentSet[1:] 	
	Sets.append(set(currentSet))

    setType = []

    for g,elem in enumerate(Sets):
	print "Set", g, "analyzed."
	for h,item in enumerate(Sets):
		if h == 0:
			setType.append(3)
		if g != h:
			# if identical
			if item == elem:
				print 0
				#setType[-1] = 0 
				break
			# if proper subset	
			elif len(elem - item) == 0 and len(item-elem) > 0:
				print 1
				#setType[-1] = 1
				break

			# if superset
			elif len(elem-item) > 0 and len(item-elem) == 0:
				print 2
				#setType[-1] = 2
				break
		
				 	

    #for g,item in enumerate(Sets):

	#fp3.write("%s %s\n" %(varNum[g], setType[g]))

    #fp3.close()
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
   
   f= open("../results/text/NSVs.txt","w")
   
   f.write("%s\n" %nSV)
   f.close()

   #WriteMappings(DisjFM, "FragMapD.txt")
 
   
