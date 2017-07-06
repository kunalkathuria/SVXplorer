# 12/16/16
# Set Cover Algorithms to Minimize number of SVs that account for all PEMs pointing to them

import math
import os
import sys
from collections import Counter

FILE = "../results/text/VariantMap.txt"

disjThresh = int(sys.argv[1])
MAP_THRESH = int(sys.argv[2])
ignoreRD = 1
RD_VAR_NUM = int(sys.argv[3]) #3000000
RD_Frag = int(sys.argv[4])
MQ_HASH = {}

def formHash():

	f=open("../results/text/All_Discords_P_S.txt","r")
	for line in f:
		temp = int(line.split()[0])
		temp2 = int(line.split()[-1])
		#include only primary alignments; assume if fragment unique for a variant, then it is primary
		#we are ensured that primary alignment of the read clears MAP_THRESH, and relative AS of secondary is very close by default
		if temp not in MQ_HASH and temp2 >= MAP_THRESH:
			MQ_HASH[temp] = 1
	
def DisjointAlg(fragmentList, nSetsR):

    N_potSVs = 0
    counter = 0
    taggedList=[] # keeps track of shared elements
    setCover = []
    listSize = []
    disjointness = [0]*nSetsR
    varNum = []
    #results = [0]*(max(fragmentList) + 1)

    # fragment numbering starts with 1

    #claimed = [0]*(1+max(fragmentList)) # change to n unique frags if saves time

    formHash()
    fp=open(FILE,"r")
    fp2=open("../results/text/DisjSetCover_S.txt","w")
    fp3=open("../results/text/UniqueSupport.txt","w")
    fp4=open("../results/text/mq0_variants.txt","w")

    x = Counter(fragmentList)
    # Obtain disjointness count
    # Assume any given fragment only appears once in set list
    for line in fp:

        currentSet = map(int, line.split())
        varNum.append(currentSet[0])
	currentSet = currentSet[1:] 	

        for elem in currentSet:
		
		#pick those supported by RD automatically		
		if elem >= RD_Frag:
		
			disjointness[counter] = disjThresh

		#if SR, no secondaries so is above MQ_THRESH automatically
		elif elem < 0 or (x[elem] == 1 and elem in MQ_HASH):

	            disjointness[counter]+=1

		if disjointness[counter] == disjThresh:
			break
	counter+=1

    for g,item in enumerate(disjointness):


	# RD CNVs do not use same support threshold
        if item == disjThresh and varNum[g] < RD_VAR_NUM:
           fp2.write("%s\n" %varNum[g]) 
	   N_potSVs+=1
	elif varNum[g] < RD_VAR_NUM:
	   fp4.write("%s\n" %varNum[g])

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
 
   
