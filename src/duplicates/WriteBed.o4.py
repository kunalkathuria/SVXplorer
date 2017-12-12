# Can include check to see if map type = 1 if other documented end of random insertion was claimed or not. If not, belongs to another variant.
# Use formatBP.py to make bedpe documentation amenable to this reading format

import sys
import pysam
from collections import Counter
BAM = sys.argv[2]
DEL_THRESH = float(sys.argv[3])
DUP_THRESH = float(sys.argv[4])
DEL_THRESH2 = float(sys.argv[5])
DUP_THRESH2 = float(sys.argv[6])
PILEUP_THRESH = float(sys.argv[7])
MIN_PILEUP_THRESH = float(sys.argv[8])
RPT_REGIONS_FILE =  sys.argv[9]
GOOD_REG_THRESH=.5
chrHash = {}

def formHash():

	print "Forming pile-up hash table..."
        prevTID = "*" # invalid value
        fo=open(RPT_REGIONS_FILE,"r")

        for line in fo:
                line_s = line.split()
                currentTID = line_s[0]
                if currentTID != prevTID:
			#global var but okay here since adding to it, not redefining
                        chrHash[currentTID] = {}

                for x in range(int(line_s[1]), int(line_s[2])+1):
                        chrHash[currentTID][x] = 1

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
		
        y = Counter(DisjSC)
	counter = -1
	samfile = pysam.AlignmentFile(BAM, "rb" )
	if RPT_REGIONS_FILE != "none":
		formHash()

	for line2 in f1:
		counter+=1
		if counter % 100 == 0:
			print "Writing Variant", counter
		line2_split = line2.split()
		num = int(line2_split[0])

		if y[num] > 0:
                   
			if line2_split[11].find("RD") == -1: 
				#print num
				#f11.write("%s\n" %line2)
				swap = 0
				GT=""
				if (line2_split[1] == "DEL_INS" or line2_split[1] == "DEL" or line2_split[1][0:2]== "TD") and int(line2_split[4]) < int(line2_split[6]):
					chr_n = line2_split[2]
					start = int(line2_split[4])
					stop = int(line2_split[6])
					covLoc = 0
					counter2 = 0
					x = start
					for pileupcolumn in samfile.pileup(chr_n, start, stop):

						if chr_n in chrHash and x in chrHash[chr_n]:
							covLoc = covLoc + pileupcolumn.n
							counter2+=1
							if counter2 > PILEUP_THRESH:
								break
						x+=1
						if x-start > 5*PILEUP_THRESH:
							break

					if counter2 > MIN_PILEUP_THRESH and (counter2 > PILEUP_THRESH or counter2 > GOOD_REG_THRESH*(stop-start)):	
						covLoc = 1.0*(covLoc/counter2)
				
						if line2_split[1][0:2] == "TD" and covLoc/COVERAGE > DUP_THRESH:
						
							print "TD confirmed (pileup)"	
							line2_split[1] = "TD"
							#if covLoc/COVERAGE > DUP_THRESH2:
								#GT="GT:1/1"
						
						elif line2_split[1][0:2] == "TD" and covLoc/COVERAGE < DEL_THRESH:
					
							line2_split[1] = "BND"
				
						elif line2_split[1][0:3] == "DEL" and covLoc/COVERAGE < DEL_THRESH:
						
							print "DEL confirmed (pileup)"
							line2_split[1] = "DEL"
							if covLoc/COVERAGE < DEL_THRESH2:
								GT="GT:1/1"
							elif covLoc/COVERAGE > 3*DEL_THRESH2:	
								GT="GT:0/1"

						elif line2_split[1][0:3] == "DEL" and covLoc/COVERAGE > DUP_THRESH:

							# since bp3 = -1, this will be written as a BND event
							line2_split[1] = "INS"
					
				# can add this and regular TDs also
				elif len(line2_split[1]) > 2 and (line2_split[1][0:3] == "INS" or line2_split[1] == "INS_I") and int(line2_split[7]) < int(line2_split[9]):

					chr_n = line2_split[5]
					start = int(line2_split[7])
					stop = int(line2_split[9])
					covLoc = 0
                                        counter2 = 0
					x = start

                                        for pileupcolumn in samfile.pileup(chr_n, start, stop):

                                                if chr_n in chrHash and x in chrHash[chr_n]:
                                                        covLoc = covLoc + pileupcolumn.n
                                                        counter2+=1
                                                        if counter2 > PILEUP_THRESH:
                                                                break
                                                x+=1
						if x-start > 5*PILEUP_THRESH:
                                                        break

                                        if counter2 > MIN_PILEUP_THRESH and (counter2 > PILEUP_THRESH or counter2 > GOOD_REG_THRESH*(stop-start)):

						covLoc = 1.0*(covLoc/counter2)

						if covLoc/COVERAGE < DEL_THRESH:
							line2_split[1] = "INS_C_P"
						elif covLoc/COVERAGE < 1.1:
							line2_split[1] = "INS_C"

				elif len(line2_split[1]) > 4 and line2_split[1][0:4] == "INS_C":

					chr_n = line2_split[5]
					start1 = int(line2_split[4])
					stop1 = int(line2_split[6])
					start2 = int(line2_split[7])
					stop2 = int(line2_split[9])
					if start1 > stop2:

						start1 = int(line2_split[10])
						stop2 = int(line2_split[3])

					if start1 < start2 < stop2:

						covLoc = 0
                                                counter2 = 0
						x = start1
                                        	for pileupcolumn in samfile.pileup(chr_n, start1, stop1):

                                                	if chr_n in chrHash and x in chrHash[chr_n]:
                                                        	covLoc = covLoc + pileupcolumn.n
                                                        	counter2+=1
                                                        	if counter2 > PILEUP_THRESH:
                                                                	break
                                                	x+=1
							if x-start > 5*PILEUP_THRESH:
                                                        	break
						covLoc = 1.0*(covLoc/counter2)

						covLoc2 = 0
						counter3 = 0
						x = start2
                                                for pileupcolumn in samfile.pileup(chr_n, start2, stop2):

                                                        if chr_n in chrHash and x in chrHash[chr_n]:
                                                                covLoc2 = covLoc2 + pileupcolumn.n
                                                                counter3+=1
                                                                if counter3 > PILEUP_THRESH:
                                                                        break
                                                        x+=1
							if x-start > 5*PILEUP_THRESH:
                                                        	break
						covLoc2 = 1.0*(covLoc2/counter3)

						# don't use pile up to sort v small variants
						if counter2 > MIN_PILEUP_THRESH and (counter2 > PILEUP_THRESH or counter2 > GOOD_REG_THRESH*(stop-start)):

						   if counter3 > MIN_PILEUP_THRESH and (counter3 > PILEUP_THRESH or counter3 > GOOD_REG_THRESH*(stop-start)):

							if covLoc/COVERAGE < DUP_THRESH and covLoc2/COVERAGE > 1.2*DUP_THRESH:
								line2_split[1] = "INS"
								#$ call another diploid deletion here

							elif covLoc2/COVERAGE < DEL_THRESH or covLoc/COVERAGE < DEL_THRESH:
								line2_split[1] = "BND_INS_C"

							elif (line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I") and covLoc/COVERAGE > DUP_THRESH and covLoc2/COVERAGE < DUP_THRESH:
								swap = 1
								if covLoc/COVERAGE > 1.2*DUP_THRESH:
									line2_split[1] = "INS"
								#$ call another diploid deletion here
					else:
						line2_split[1] = "BND_INS_C"	

			# If not insertion as bp3 is not set
			if line2_split[1] == "DEL":
                            
                            f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"DEL",GT))
                        #$Comment out second condition and next elif if leads to low precision due to SR TD_I's and INS_I's.    
                        elif line2_split[1] == "TD":

			    [bp1_s, bp1_e] = min(int(line2_split[3]),int(line2_split[6])), min(int(line2_split[4]), int(line2_split[7]))
			    [bp2_s, bp2_e] = max(int(line2_split[3]),int(line2_split[6])), max(int(line2_split[4]), int(line2_split[7]))
 
			    f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], bp1_s, bp1_e, line2_split[5], bp2_s, bp2_e,"TD",GT))

			#See whether this works better for PE only vs PE and SR both-- unlikely to have inv TD...
			elif line2_split[1] == "INS_I" and line2_split[2] == line2_split[5] and line2_split[8] == "-1" and line2_split[11][:2] == "PE":
			
			    [bp1_s, bp1_e] = min(int(line2_split[3]),int(line2_split[6])), min(int(line2_split[4]), int(line2_split[7]))
                            [bp2_s, bp2_e] = max(int(line2_split[3]),int(line2_split[6])), max(int(line2_split[4]), int(line2_split[7]))

                            f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], bp1_s, bp1_e, line2_split[5], bp2_s, bp2_e,"TD_INV",GT))

			elif line2_split[1] == "INV":
                            f15.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1],GT))

			# this is read from cluster file, so has 2 bps only; INS_C w/ only 2 clusters supporting; INV w/ only 1
                        elif line2_split[1] == "Unknown" or line2_split[1] == "INS_POSS" or line2_split[1] == "TD_I" or line2_split[1] == "INV_POSS" or ( (line2_split[1] == "INS" or line2_split[1] == "INS_I" or line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I") and (line2_split[9] == "-1" or line2_split[6] == "-1")) or (line2_split[1] == "INS_C" and line2_split[12] == "2"):

                            f17.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"BND", line2_split[1],GT))

                        elif len(line2_split[1]) > 2 and line2_split[1][0:3] == "INS":

			    # two lines for insertion as in bedpe format; bp 1 and bp3 were flanks of bp2 by convention in INS_C classification unless confirmed further
			    [bp1_s, bp1_e]= int(line2_split[3]), int(line2_split[4])
			    [bp2_s, bp2_e] = int(line2_split[6]),int(line2_split[7])
                            [bp3_s, bp3_e] = int(line2_split[9]),int(line2_split[10])
			    
			    if swap:
				temp = bp1_e
				bp1_e = bp3_e
				bp3_e = temp
				bp1_s = bp3_s
				
			    if not (line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I") and bp2_s > bp3_e and not swap:
				temp = bp2_s
				bp2_s = bp3_e
				bp3_e = temp 

			    # it was INS_C
			    elif bp2_s > bp3_e:

				bp2_s = bp3_s
				bp3_e = bp2_e
						
                            f16.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line2_split[5], bp2_s, bp3_e, line2_split[2], bp1_s, bp1_e, line2_split[1],GT))

                        else:
				f17.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"BND", line2_split[1],GT))
 
	    	if line2_split[1] == "DEL_uc":

                            f13b.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1],GT))
	
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



                                        
                        
