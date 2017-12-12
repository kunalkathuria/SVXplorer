# 12/16/16
#Filter variants having minimum unique support for each support category (SR, PE, MIXED)

import sys
import math
from collections import Counter

FILE = "../results/text/VariantMap.txt"

def readStats():
	f=open("../results/text/bam_stats.txt","r")
	for i,line in enumerate(f):
		if i==3:	
			break
		elif i==2:
			sigIL = float(line.split()[0])
	return float(line.split()[0]), sigIL

#disjThresh = int(sys.argv[1])
PE_THRESH_MAX = 6
SR_THRESH_MAX = 6
PE_THRESH_MIN = 3
SR_THRESH_MIN = 3
MIX_THRESH = 4
PE_LOW = 3
PE_HIGH = 6
COVG_LOW = 5
COVG_HIGH = 50
SR_LOW = 4
SR_HIGH = 6
IL_BOOST = 0
IL_LOW1 = 20
IL_LOW2 = 30
IL_HIGH1 = 90
IL_HIGH2 = 105

[COVG,SIG_IL] = readStats()

#LINEAR MODEL TO CALCULATE SUPPORT THRESHOLDS
if IL_LOW1 < SIG_IL < IL_LOW2:
	IL_BOOST = -.5
elif IL_HIGH1 <= SIG_IL <= IL_HIGH2:
	IL_BOOST = 1
elif SIG_IL > IL_HIGH2:
	IL_BOOST = 2

SR_THRESH = 4#math.floor(SR_LOW + (COVG-COVG_LOW)*1.0*(SR_HIGH - SR_LOW)/(COVG_HIGH - COVG_LOW))
PE_THRESH = 4#round(IL_BOOST + PE_LOW + (COVG-COVG_LOW)*1.0*(PE_HIGH - PE_LOW)/(COVG_HIGH - COVG_LOW))
if PE_THRESH > PE_THRESH_MAX:
	PE_THRESH = PE_THRESH_MAX
if SR_THRESH > SR_THRESH_MAX:
        SR_THRESH = SR_THRESH_MAX
if PE_THRESH < PE_THRESH_MIN:
        PE_THRESH = PE_THRESH_MIN
if SR_THRESH < SR_THRESH_MIN:
        SR_THRESH = SR_THRESH_MIN
print "SR, PE threshes, covg are:", SR_THRESH, PE_THRESH, COVG

MAP_THRESH = int(sys.argv[2])
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
    pickV = [0]*nSetsR
    varNum = []
    #results = [0]*(max(fragmentList) + 1)

    # fragment numbering starts with 1

    #claimed = [0]*(1+max(fragmentList)) # change to n unique frags if saves time

    formHash()
    fp=open(FILE,"r")
    fp2=open("../results/text/DisjSetCover_S.txt","w")
    fp3=open("../results/text/All_Variants.txt","r")
    fp4=open("../results/text/mq0_variants.txt","w")
    fp5=open("../results/text/uniqueCount.txt","w")

    print "In SC fxn."
    x = Counter(fragmentList)
    # Obtain disjointness count
    # Assume any given fragment only appears once in set list
    for line in fp:

        currentSet = map(int, line.split())
        varNum.append(currentSet[0])
	currentSet = currentSet[1:] 
	
	#3 variables for variant type for clarity	
	SRvar = 0
	PEvar = 0
	mixedVar = 0
	disjThresh = -1

	for line in fp3:
		if line.split()[11].find("PE") == -1 and line.split()[11].find("SR") != -1:
			SRvar = 1
			disjThresh = SR_THRESH
		elif line.split()[11].find("PE") != -1 and line.split()[11].find("SR") == -1:
			PEvar = 1
			disjThresh = PE_THRESH
		elif line.split()[11].find("PE") != -1 and line.split()[11].find("SR") != -1:
			mixedVar = 1
			disjThresh = MIX_THRESH

		#print line, disjThresh
		break

	
        for elem in currentSet:
		
		#pick those supported by RD automatically		
		if elem >= RD_Frag:
		
			if x[elem] == 1:
				disjointness[counter]+=1
			pickV[counter] = 2

		#if SR, no secondaries so is above MQ_THRESH automatically
		elif x[elem] == 1:

	            disjointness[counter]+=1

		    #require at least 1 MQ > some_thresh support
		    if elem in MQ_HASH or elem < 0:
			pickV[counter] = 1

		if disjointness[counter] == disjThresh:
			break
	counter+=1

    fp3.seek(0)
    for g,item in enumerate(disjointness):

	for line in fp3:
                if line.split()[11].find("PE") == -1 and line.split()[11].find("SR") != -1:
                        SRvar = 1
			disjThresh = SR_THRESH
                elif line.split()[11].find("PE") != -1 and line.split()[11].find("SR") == -1:
                        PEvar = 1
			disjThresh = PE_THRESH
                elif line.split()[11].find("PE") != -1 and line.split()[11].find("SR") != -1:
                        mixedVar = 1
			disjThresh = MIX_THRESH

                break

	# RD CNVs do not use same support threshold
        if (item >= disjThresh and pickV[g] == 1) or (item >= 1 and pickV[g] == 2):
           fp2.write("%s\n" %varNum[g]) 
	   N_potSVs+=1
	#make list of variants with unique support till thresh
	if varNum[g] < RD_VAR_NUM:
	   fp5.write("%s %s\n" %(varNum[g], disjointness[g]))

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
 
   
