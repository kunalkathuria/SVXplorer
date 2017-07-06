# Can include check to see if map type = 1 if other documented end of random insertion was claimed or not. If not, belongs to another variant.
# Use formatBP.py to make bedpe documentation amenable to this reading format

import sys
import pysam
from collections import Counter
BAM = sys.argv[2]
DEL_THRESH = float(sys.argv[3])
DUP_THRESH = float(sys.argv[4])

def file_len(f):
    
    for i, l in enumerate(f):
        pass
    return i + 1

if __name__ == "__main__":

	f1 = open("../results/text/All_Variants.txt","r")
        f8 = open(sys.argv[1],"r")
	f11 = open("../results/text/allPositives.txt","w")
	f12 = open("../results/text/allPositives.bedpe","w")
	f13 = open("../results/text/deletions.bedpe","w")
	f13b = open("../results/text/deletions_01.bedpe","w")
        f14 = open("../results/text/tandemDuplications.bedpe","w")
        f15 = open("../results/text/inversions.bedpe","w")
        f16 = open("../results/text/insertions.bedpe","w")
        f17 = open("../results/text/unknowns.bedpe","w")

	fo = open("../results/text/bam_stats.txt","r")		
	for i,line in enumerate(fo):
		if i==3:
			break
	COVERAGE = float(line[:-1])
	print "Coverage is:", COVERAGE

	DisjSC = []

        for line in f8:
				
                DisjSC.append(int(line))
		
        x = Counter(DisjSC)
	#print x
	#print Claimed
	counter = -1
	samfile = pysam.AlignmentFile(BAM, "rb" )

	for line2 in f1:
		counter+=1
		if counter % 100 == 0:
			print "Writing Variant", counter
		line2_split = line2.split()
		num = int(line2_split[0])

		if x[num] > 0:
                    
			#print num
			#f11.write("%s\n" %line2)
			write = 0
			if (line2_split[1] == "DEL_INS" or (line2_split[1]== "TD_I" and line2_split[11] == "SR")) and int(line2_split[4]) < int(line2_split[6]):
				chr_n = line2_split[2]
				start = int(line2_split[4])
				stop = int(line2_split[6])
				covLoc = 0
				counter2 = 0
				for pileupcolumn in samfile.pileup(chr_n, start, stop):
					covLoc = covLoc + pileupcolumn.n
					counter2+=1

				if counter2:	
					covLoc = 1.0*(covLoc/counter2)
			
				if line2_split[1] == "TD_I" and covLoc != 0 and covLoc/COVERAGE > DUP_THRESH:
				
					write = 1		
					print "TD confirmed (pileup)"	

				elif covLoc != 0 and covLoc/COVERAGE < DEL_THRESH:
				
					write = 1
					print "DEL confirmed (pileup)"

			#elif len(line2_split[1]) > 2 and line2_split[1][0:3] == "INS" and line2_split[1] != "INS_C" and line2_split[1] != "INS_C_I"):

				#chr_n = line2_split[5]
                                #start = line2_split[7]
                                #stop = line2_split[9]
                                #samfile = pysam.AlignmentFile(BAM, "rb" )
                                #covLoc = 0
                                #counter = 0
                                #for pileupcolumn in samfile.pileup(chr_n, start, stop):
                                        #covLoc = covLoc + pileupcolumn.n
                                        #counter+=1
                                #covLoc = 1.0*(covLoc/counter)

				#if covLoc/COVERAGE > DUP_THRESH:
					#write = 1

			# If not insertion as bp3 is not set
			if line2_split[1] == "DEL" or (line2_split[1] == "DEL_INS" and (line2_split[11].find("RD") != -1 or write)):
                            
                            f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"DEL"))
                        #$Comment out second condition and next elif if leads to low precision due to SR TD_I's and INS_I's.    
                        elif line2_split[1] == "TD" or (line2_split[1] == "TD_I" and (line2_split[11][:2] == "PE" or write)):

			    [bp1_s, bp1_e] = min(int(line2_split[3]),int(line2_split[6])), min(int(line2_split[4]), int(line2_split[7]))
			    [bp2_s, bp2_e] = max(int(line2_split[3]),int(line2_split[6])), max(int(line2_split[4]), int(line2_split[7]))
 
			    f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], bp1_s, bp1_e, line2_split[5], bp2_s, bp2_e,"TD"))

			#See whether this works better for PE only vs PE and SR both-- unlikely to have inv TD in chromosome...
			elif line2_split[1] == "INS_I" and line2_split[2] == line2_split[5] and line2_split[8] == "-1" and line2_split[11][:2] == "PE":
			
			    [bp1_s, bp1_e] = min(int(line2_split[3]),int(line2_split[6])), min(int(line2_split[4]), int(line2_split[7]))
                            [bp2_s, bp2_e] = max(int(line2_split[3]),int(line2_split[6])), max(int(line2_split[4]), int(line2_split[7]))

                            f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], bp1_s, bp1_e, line2_split[5], bp2_s, bp2_e,"TD_INV"))

			elif line2_split[1] == "INV":
                            f15.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1]))

			# this is read from cluster file, so has 2 bps only; INS_C w/ only 2 clusters supporting; INV w/ only 1
                        elif line2_split[1] == "Unknown" or line2_split[1] == "INS_POSS" or line2_split[1] == "TD_I" or line2_split[1] == "INV_POSS" or ( (line2_split[1] == "INS" or line2_split[1] == "INS_I" or line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I") and (line2_split[9] == "-1" or line2_split[6] == "-1")) or (line2_split[1] == "INS_C" and line2_split[12] == "2"):

                            f17.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"BND", line2_split[1]))

                        elif len(line2_split[1]) > 2 and line2_split[1][0:3] == "INS":

			    # two lines for insertion as in bedpe format; bp 1 and bp3 were flanks of bp2 by convention in INS_C classification unless confirmed further
			    [bp2_s, bp2_e] = min(int(line2_split[6]),int(line2_split[9])), min(int(line2_split[7]), int(line2_split[10]))
                            [bp3_s, bp3_e] = max(int(line2_split[6]),int(line2_split[9])), max(int(line2_split[7]), int(line2_split[10]))
			   
			    # in case of an artefact in insertion classfication
			    if bp2_e > bp3_s:
				temp = bp2_e
				bp2_e = bp3_s
				bp3_s = temp 
						
                            f16.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line2_split[5], bp2_e, bp3_s, line2_split[2], line2_split[3], line2_split[4], line2_split[1]))

                        else:
				f17.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"BND", line2_split[1]))
 
	    	if line2_split[1] == "DEL_uc":

                            f13b.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1]))
	
	f1.close()
	f8.close()
        f11.close()
	f12.close()
	f13.close()
	f14.close()
	f15.close()
	f16.close()
	f17.close()
	samfile.close()



                                        
                        
