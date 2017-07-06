# 12/16/16
# Set Cover Algorithms to Minimize number of SVs that account for all PEMs pointing to them

import math
import os
import sys
from collections import Counter

FILE = "../results/text/VariantMap.txt"

nSets = 5
nFragments = 5
maxSetEntries = 4
ignoreRD = 1
RD_VAR_NUM = 3000000


def DisjointAlg(fragmentList, nSetsR):

    N_potSVs = 0
    counter = 0
    taggedList=[] # keeps track of shared elements
    setCover = []
    listSize = []
    disjointness = [0]*nSetsR
    varNum = []
    varSupport = []
	
    # fragment numbering starts with 1

    #claimed = [0]*(1+max(fragmentList)) # change to n unique frags if saves time

    fp=open(FILE,"r")
    fp2=open("../results/text/DisjSetCover_Unique.txt","w")
    fp3=open("../results/text/UniqueSupport.txt","w")

    x = Counter(fragmentList)
	
    # Obtain disjointness count
    # Assume any given fragment only appears once in set list
    for line in fp:

        currentSet = map(int, line.split())
        varNum.append(currentSet[0])
	currentSet = currentSet[1:] 	

        for elem in currentSet:
		
		if x[elem] == 1:

	            disjointness[counter]+=1

	varSupport.append(len(currentSet))
	counter+=1

    for g,item in enumerate(disjointness):

	fp3.write("%s %s %s\n" %(varNum[g], varSupport[g], disjointness[g]))

	# if all fragments supporting variant are uniquely supporting it
        if varSupport[g] == disjointness[g] and varNum[g] < RD_VAR_NUM:
           fp2.write("%s\n" %varNum[g])
	   N_potSVs+=1
	if varNum[g] >= RD_VAR_NUM and not ignoreRD:
           fp2.write("%s\n" %varNum[g]) 
	   N_potSVs+=1

    return N_potSVs

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
 
   
