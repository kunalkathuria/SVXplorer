# Can include check to see if map type = 1 if other documented end of random insertion was claimed or not. If not, belongs to another variant.
# Use formatBP.py to make bedpe documentation amenable to this reading format

import sys
import pysam
from collections import Counter
from bitarray import bitarray
BAM = sys.argv[2]
DEL_THRESH = float(sys.argv[3])
DUP_THRESH = float(sys.argv[4])
DEL_THRESH2 = float(sys.argv[5])
DUP_THRESH2 = float(sys.argv[6])
PILEUP_THRESH = float(sys.argv[7])
MIN_PILEUP_THRESH = float(sys.argv[8])
RPT_REGIONS_FILE =  sys.argv[9]
GOOD_REG_THRESH=.8 # to trust pile-up depth in given region, this percentage should return data
SECOND_SWEEP_THRESH=.75
MINPU_MULT=1
VER_BUFFER=100
PE_THRESH_H=3

class Variant(object):

    #$ add inv, ins also (bp 3)
    def __init__(self):
        self.bp = -1

        self.tid = -1

    def __str__(self):

        return "%s %s" %(self.bp, self.tid)

def formHash(RDL, chrHash):

	print "Forming pile-up hash table..."
        fo=open(RPT_REGIONS_FILE,"r")
	prev_start = -1
	prev_stop = -1

        for k,line in enumerate(fo):
		
                line_s = line.split()
                currentTID = line_s[0]
		start = int(line_s[1])
		stop = int(line_s[2])+1
		#print currentTID, start

		if currentTID not in chrHash:
				prev_stop = -1
                                chrHash[currentTID] = bitarray()
                                #bed file is 1-based
                                for y in range(1,start):
                                        chrHash[currentTID].append(0)

		for x in range(start, stop):

			chrHash[currentTID].append(1)

		# make hash table of unreliable regions if greater than RDL (almt would be doubtful there)
		if prev_stop != -1 and currentTID == prevTID:

			addBit = 1
			if start - prev_stop > RDL:
				addBit = 0
	                for x in range(prev_stop, start):
				if currentTID not in chrHash:
					chrHash[currentTID] = bitarray()
				
				chrHash[currentTID].append(addBit)


		prev_start = start
		prev_stop = stop
		prevTID = currentTID

	print "Done"

def file_len(f):
    
    for i, l in enumerate(f):
        pass
    return i + 1

if __name__ == "__main__":

	f1 = open("../results/text/All_Variants.txt","r")
	f2=open("../results/text/VariantMap.txt","r")
        f8 = open(sys.argv[1],"r")
	f11 = open("../results/text/allPositives.txt","w")
	f12 = open("../results/text/allPositives.bedpe","w")
	f13 = open("../results/text/deletions.bedpe","w")
	f13b = open("../results/text/deletions_01.bedpe","w")
        f14 = open("../results/text/tandemDuplications.bedpe","w")
        f15 = open("../results/text/inversions.bedpe","w")
        f16 = open("../results/text/insertions.bedpe","w")
        f17 = open("../results/text/unknowns.bedpe","w")

	print "Writing final bedpe files using coverage information..."
	fo = open("../results/text/bam_stats.txt","r")		
	RDL = -1
	SD = -1
	for i,line in enumerate(fo):
		if i == 0:
			RDL=float(line[:-1])
		elif i==2:
			SD=float(line[:-1])
		elif i==3:
			break
	COVERAGE = float(line[:-1])
	print "Coverage is:", COVERAGE
	UNIV_VAR_THRESH=100
	INS_VAR_THRESH=50
	SR_DEL_THRESH=250
	PE_DEL_THRESH_S=250
	PE_DEL_THRESH_L=200
	SD_S = 15
	SD_L = 20
	# calculate min PE size based on insert length standard deviation under simple empirical linear model, and nature of small calls that fits generally well. Aligner tends to mark concordants as discordants when SD is small unless distr v good, seemingly sharp effects around sig_IL of 20 and lower for Poisson and other related distributions.
	PE_DEL_THRESH=PE_DEL_THRESH_S + int((SD-SD_S)*(PE_DEL_THRESH_L-PE_DEL_THRESH_S)/(SD_L-SD_S))
	if PE_DEL_THRESH > PE_DEL_THRESH_S:
		PE_DEL_THRESH = PE_DEL_THRESH_S
	elif PE_DEL_THRESH < PE_DEL_THRESH_L:
		PE_DEL_THRESH = UNIV_VAR_THRESH

	DisjSC = []
	chrHash = {}

        for line in f8:
				
                DisjSC.append(int(line))
		
        y = Counter(DisjSC)
	counter = -1
	samfile = pysam.AlignmentFile(BAM, "rb" )
	if RPT_REGIONS_FILE != "none":
		formHash(RDL, chrHash)

	nv_PEposs = 0
	nv_PEconf = 0
	nv_SRposs = 0
	nv_SRconf = 0

	support = 0
	for entry in f2:
		support = len(entry.split()) - 1
		break

	for line2 in f1:
		counter+=1
		if counter % 100 == 0:
			print "Writing Variant", counter
		line2_split = line2.split()
		num = int(line2_split[0])

		if y[num] > 0:
               
			if line2_split[1][0:3] == "TD" or (len(line2_split) > 2 and (line2_split[1][0:3] == "DEL" or line2_split[1][0:3] == "INV")):
				if int(line2_split[6])-int(line2_split[3]) < UNIV_VAR_THRESH:
					continue
			elif len(line2_split[1]) > 2 and line2_split[1][0:3] == "INS" and not (line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I"):

				if (0 < int(line2_split[10])-int(line2_split[6]) < INS_VAR_THRESH) or (line2_split[2] == line2_split[5] and ((0 < int(line2_split[7])-int(line2_split[3]) < UNIV_VAR_THRESH) or (0 < int(line2_split[4])-int(line2_split[9]) < UNIV_VAR_THRESH))):
                                        continue

  			elif line2_split == "INS_C" or line2_split == "INS_C_I":
	
				if (0 < int(line2_split[10])-int(line2_split[6]) < INS_VAR_THRESH):
                                        continue
 
			if line2_split[11].find("RD") == -1: 
				#f11.write("%s\n" %line2)
				swap = 0
				GT=""
				if (line2_split[1] == "DEL_INS" or line2_split[1] == "DEL" or line2_split[1][0:2]== "TD") and int(line2_split[4]) + MINPU_MULT*MIN_PILEUP_THRESH < int(line2_split[6]):
					chr_n = line2_split[2]
					gap = int(line2_split[6])-int(line2_split[4])
					start = .25*gap + int(line2_split[4])
					stop = min(start+.5*gap,start +3*PILEUP_THRESH)
					covLoc = 0
					counter2,counter_v1,counter_v2,confM, confL, confR = 0,0,0,0,0,0
					#print "In DEL condition"
					if stop > start:					
						for pileupcolumn in samfile.pileup(chr_n, start, stop):

                                                        if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
                                                                #print "DEL hash", pileupcolumn.pos, pileupcolumn.n
                                                                covLoc = covLoc + pileupcolumn.n
                                                                counter2+=1
                                                                if counter2 > PILEUP_THRESH:
                                                                        break

					if (counter2 > MIN_PILEUP_THRESH and (counter2 > GOOD_REG_THRESH*(stop-start) or counter2 > PILEUP_THRESH)):
                                                        confM = 1
							covLoc = (1.0*covLoc)/(1.0*counter2)

					if line2_split[1][0:2] != "TD":	

						ver_stop1 = int(line2_split[3])	
						ver_start1 = ver_stop1 - VER_BUFFER 
						ver_start2 = int(line2_split[7])
						ver_stop2 = ver_start2 + VER_BUFFER
						covLoc_v1, covLoc_v2 = 0, 0

						if stop > start:
						 for pileupcolumn in samfile.pileup(chr_n, ver_start1, ver_stop1):

							if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
								#print "DEL hash", pileupcolumn.pos, pileupcolumn.n
								covLoc_v1 = covLoc_v1 + pileupcolumn.n
								counter_v1+=1
								if counter_v1 > PILEUP_THRESH:
									break
						if stop > start:
						 for pileupcolumn in samfile.pileup(chr_n, ver_start2, ver_stop2):

							if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
								#print "DEL hash", pileupcolumn.pos, pileupcolumn.n
								covLoc_v2 = covLoc_v2 + pileupcolumn.n
								counter_v2+=1
								if counter_v2 > PILEUP_THRESH:
									break
						#print 1.0*covLoc/(.001+counter2), "is local coverage", counter2, counter3, covLoc

						if (counter_v1 > MIN_PILEUP_THRESH and (counter_v1 > GOOD_REG_THRESH*(ver_stop1-ver_start1) or counter_v1 > PILEUP_THRESH)):
							confL = 1
							covLoc_v1 = (1.0*covLoc_v1)/(1.0*counter_v1) 

						if (counter_v2 > MIN_PILEUP_THRESH and (counter_v2 > GOOD_REG_THRESH*(ver_stop2-ver_start2) or counter_v2 > PILEUP_THRESH)):
							confR = 1
							covLoc_v2 = (1.0*covLoc_v2)/(1.0*counter_v2)

					if confM == 1 or (confL == 1 and confR == 1):

						#print confM, confL, confR	
						if line2_split[1][0:2] == "TD" and confM == 1 and covLoc/COVERAGE > DUP_THRESH:
						
							print "TD confirmed (pileup)"	
							line2_split[1] = "TD"
							#if covLoc/COVERAGE > DUP_THRESH2:
								#GT="GT:1/1"
							if line2_split[11].find("PE") == -1 and line2_split[11].find("SR") == 1:
								nv_SRconf+=1
							elif line2_split[11].find("PE") == 1:
								nv_PEconf+=1

						elif line2_split[1][0:2] == "TD" and confM == 1:
							
							if line2_split[11].find("PE") == -1 and line2_split[11].find("SR") == 1:
								nv_SRposs+=1
							elif line2_split[11].find("PE") == 1:
								nv_PEposs+=1

							if covLoc/COVERAGE < 1.0:
								line2_split[1] = "BND"

				
						elif line2_split[1][0:3] == "DEL" and ((confM == 1 and covLoc/COVERAGE < DEL_THRESH) or (confL == 1 and covLoc_v1/COVERAGE < DEL_THRESH and confR == 1 and covLoc_v2/COVERAGE < DEL_THRESH)):
						
							print "DEL confirmed (pileup)"
							#print confM, covLoc/COVERAGE, confL, covLoc_v1/COVERAGE, confR, covLoc_v2/COVERAGE
							line2_split[1] = "DEL"
							if covLoc/COVERAGE < DEL_THRESH2:
								GT="GT:1/1"
							elif covLoc/COVERAGE > 3*DEL_THRESH2:	
								GT="GT:0/1"
					
							#if confM and covLoc/COVERAGE > 1.0:
                                                                # since bp3 = -1, this will be written as a BND event
                                                                #line2_split[1] = "INS"
								#print "Oops, not confirmed. Barely missed..."
	
							#if line2_split[11].find("PE") == -1 and line2_split[11].find("SR") != -1:
								#nv_SRconf+=1
								#nv_SRposs+=1
							#elif line2_split[11].find("PE") != -1:
								#nv_PEconf+=1
								#nv_PEposs+=1

						elif line2_split[1][0:3] == "DEL":

							if line2_split[11].find("PE") != -1 and line2_split[11].find("SR") == -1:
								if support < PE_THRESH_H:
									line2_split[1] = "DEL_POSS"
								#nv_SRposs+=1
							#elif line2_split[11].find("PE") != -1:
								#nv_PEposs+=1

							if confM and covLoc/COVERAGE > 1.0:
								# since bp3 = -1, this will be written as a BND event
								line2_split[1] = "INS"
					
				# can add this and regular TDs also
				elif False and len(line2_split[1]) > 2 and line2_split[11].find("PE") != -1 and (line2_split[1] == "INS" or line2_split[1] == "INS_I") and int(line2_split[7]) + MINPU_MULT*MIN_PILEUP_THRESH < int(line2_split[9]) and line2_split[8] != "-1":

					#bp2-3
					chr_n = line2_split[5]
					gap = int(line2_split[9])-int(line2_split[7])
                                        start = .25*gap + int(line2_split[7])
                                        stop = min(start+.5*gap,start +3*PILEUP_THRESH)
                                        covLoc = 0
                                        counter2 = 0
					if stop > start:
                                         for pileupcolumn in samfile.pileup(chr_n, start, stop):
						if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
                                                        covLoc = covLoc + pileupcolumn.n
                                                        counter2+=1
                                                        if counter2 > PILEUP_THRESH:
                                                                break
					 if counter2 > 0:
					 	covLoc = 1.0*covLoc/counter2
					
					if counter2 > MIN_PILEUP_THRESH and (counter2 > GOOD_REG_THRESH*(stop-start) or counter2 > PILEUP_THRESH) and covLoc/COVERAGE < 1.15:

						line2_split[1] = "INS_C_P"

					#if on same chromosome
					elif line2_split[2] == line2_split[5] and not (int(line2_split[6]) < int(line2_split[3]) < int(line2_split[10])):
						#$put all this in a function
						#bp 1-2 and 1-3 or 3-1 and 2-1 (dep on upstream vs downstream INS)
						if int(line2_split[6])-int(line2_split[4]) > 0:
							gap = int(line2_split[6])-int(line2_split[4])
							start = .25*gap + int(line2_split[4])
							gap2 = int(line2_split[9])-int(line2_split[4])
							start2 = .25*gap2 + int(line2_split[4])
						else:
							gap = int(line2_split[3])-int(line2_split[10])
							start = .25*gap + int(line2_split[10])
							gap2 = int(line2_split[3])-int(line2_split[7])
                                                        start2 = .25*gap2 + int(line2_split[7])

						stop = min(start+.5*gap,start +3*PILEUP_THRESH)
						covLoc = 0
						counter2 = 0
						stop2 = min(start2+.5*gap2,start2 +3*PILEUP_THRESH)
                                                covLoc_l = 0
                                                counter_l = 0
						#print start, stop, gap
						#print line2
						if stop > start and stop2 > start2:
						 for pileupcolumn in samfile.pileup(chr_n, start, stop):
							if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
								covLoc = covLoc + pileupcolumn.n
								counter2+=1
								if counter2 > PILEUP_THRESH:
									break


						if stop > start and stop2 > start2:
                                                 for pileupcolumn in samfile.pileup(chr_n, start2, stop2):
							if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
                                                                covLoc_l = covLoc_l + pileupcolumn.n
                                                                counter_l+=1
                                                                if counter_l > PILEUP_THRESH:
                                                                        break
						
						if counter2 > MIN_PILEUP_THRESH and (counter2 > GOOD_REG_THRESH*(stop-start) or counter2 > PILEUP_THRESH):

							covLoc = (1.0*covLoc)/(1.0*counter2)
							if covLoc/COVERAGE < DEL_THRESH:
								line2_split[1] = "DEL"
								#write deletion
								f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], int(start -.25*gap) - (int(line2_split[4]) - int(line2_split[3])), int(start-.25*gap), line2_split[2], int(start+.75*gap), int(start+.75*gap) + int(line2_split[4]) - int(line2_split[3]),"DEL",GT))
							
								if counter_l > MIN_PILEUP_THRESH and (counter_l > GOOD_REG_THRESH*(stop-start) or counter_l > PILEUP_THRESH):

									covLoc_l = (1.0*covLoc_l)/(1.0*counter_l)

									if covLoc_l/COVERAGE >= 1.0:
										line2_split[1] = "TD"
										#write TD
										f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], int(start2 -.25*gap2) - (int(line2_split[4]) - int(line2_split[3])), int(start2-.25*gap2), line2_split[2], int(start2+.75*gap2), int(start2+.75*gap2) + int(line2_split[4]) - int(line2_split[3]),"TD"))
								continue

						elif counter_l > MIN_PILEUP_THRESH and (counter_l > GOOD_REG_THRESH*(stop-start) or counter_l > PILEUP_THRESH):
								covLoc_l = (1.0*covLoc_l)/(1.0*counter_l)

								if covLoc_l/COVERAGE >= DUP_THRESH:
									line2_split[1] = "TD"
									#write TD
									f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], int(start2 -.25*gap2) - (int(line2_split[4]) - int(line2_split[3])), int(start2-.25*gap2), line2_split[2], int(start2+.75*gap2), int(start2+.75*gap2) + int(line2_split[4]) - int(line2_split[3]),"TD"))
									
									continue	


				elif False and len(line2_split[1]) > 4 and line2_split[11].find("PE") != -1 and line2_split[1][0:4] == "INS_C" and line_split[8] != "-1":

				      	del_23 = 0
                                      	del_21 = 0
                                      	dup_23 = 0
                                      	dup_21 = 0

					chr_n = line2_split[5]
                                        gap = int(line2_split[9])-int(line2_split[7])
                                        start = .25*gap + int(line2_split[7])
                                        stop = min(start+.5*gap,start +3*PILEUP_THRESH)
                                        covLoc = 0
                                        counter2 = 0
                                        if stop > start:
						start, stop = stop, start
                                        for pileupcolumn in samfile.pileup(chr_n, start, stop):
                                                if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
                                                        covLoc = covLoc + pileupcolumn.n
                                                        counter2+=1
                                                        if counter2 > PILEUP_THRESH:
                                                                break
                                        if counter2 > 0:
                                               covLoc = 1.0*covLoc/counter2

					#bp 1-2
					gap = int(line2_split[6])-int(line2_split[4])
					start = .25*gap + int(line2_split[4])
					stop = min(start+.5*gap,start +3*PILEUP_THRESH)
					covLoc3 = 0
					counter3 = 0

					if stop > start:
						start, stop = stop, start
					if line2_split[2] == line2_split[5]:
					 for pileupcolumn in samfile.pileup(chr_n, start, stop):
						if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
							covLoc3 = covLoc3 + pileupcolumn.n
							counter3+=1
							if counter3 > PILEUP_THRESH:
								break
					if counter3 > 0:
                                                covLoc3 = 1.0*covLoc3/counter3

					if counter2 > MIN_PILEUP_THRESH and (counter2 > GOOD_REG_THRESH*(stop-start) or counter2 > PILEUP_THRESH):

						if covLoc/COVERAGE > DUP_THRESH:
							dup_23 = 1
						elif covLoc/COVERAGE < DEL_THRESH:
							del_23 =1

					if counter3 > MIN_PILEUP_THRESH and (counter3 > GOOD_REG_THRESH*(stop-start) or counter3 > PILEUP_THRESH):

						if covLoc3/COVERAGE > DUP_THRESH:
                                                        dup_21 = 1
                                                elif covLoc3/COVERAGE < DEL_THRESH:
                                                        del_21 =1


					confIC = 0
					if (line2_split[1] == "INS_C" or line2_split[1] == "INS_C_I"): 

						if dup_23:
							line2_split[1]+="_P"
							confIC = 1

						elif dup_21:
							line2_split[1]+="_P"
							confIC = 1	
							swap = 1
				
					#if both bp intervals contain deletion read depth	
					if not confIC and not (int(line2_split[6]) < int(line2_split[3]) < int(line2_split[10])):
						#calc1-3 pu
#						gap = int(line2_split[3])-int(line2_split[10])
#                                                start = .25*gap + int(line2_split[10])
#						if stop > start:
#	                                                start, stop = stop, start
#        	                                 for pileupcolumn in samfile.pileup(chr_n, start, stop):
#                	                                if chr_n in chrHash and pileupcolumn.pos < len(chrHash[chr_n]) and chrHash[chr_n][pileupcolumn.pos]:
#                        	                                covLoc = covLoc + pileupcolumn.n
#                                	                        counter2+=1
#                                        	                if counter2 > PILEUP_THRESH:
#                                                	                break
#                                        	 if counter2 > 0:
#                                                	covLoc = 1.0*covLoc/counter2
#
#						dup_13 = 0
#						if counter2 > MIN_PILEUP_THRESH and (counter2 > GOOD_REG_THRESH*(stop-start) or counter2 > PILEUP_THRESH) and covLoc/COVERAGE > DUP_THRESH:
#							dup_13 =1						
#						
#						if dup_13:
#							#write TD
#							f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"DEL",GT))

						if line2_split[2] == line2_split[5] and del_23 and del_21:

							if int(line2_split[3]) < int(line2_split[6]):
						
								if int(line2_split[9]) - int(line2_split[6]) > PE_DEL_THRESH:
									f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[6], line2_split[7], line2_split[2], line2_split[8], line2_split[9],"DEL",GT))

								if int(line2_split[7]) - int(line2_split[3]) > PE_DEL_THRESH:
									f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[2], line2_split[6], line2_split[7],"DEL",GT))

							else:

								if int(line2_split[9]) - int(line2_split[6]) > PE_DEL_THRESH:
								f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[8], line2_split[9], line2_split[2], line2_split[6], line2_split[7],"DEL",GT))
                                                                f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[6], line2_split[7], line2_split[5], line2_split[3], line2_split[4],"DEL",GT))

                                                  	continue

						elif del_23:
							#write 1 DEL
							if int(line2_split[6]) < int(line2_split[9]):

                                                                f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[5], line2_split[6], line2_split[7], line2_split[5], line2_split[8], line2_split[9],"DEL",GT))

                                                        else:
                                                                f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[5], line2_split[8], line2_split[9], line2_split[5], line2_split[6], line2_split[7],"DEL",GT))
							continue
						elif line2_split[2] == line2_split[5] and del_21:
							#write 1 DEL
							if int(line2_split[3]) < int(line2_split[6]):

                                                                f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[2], line2_split[6], line2_split[7],"DEL",GT))

                                                        else:
                                                                f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[6], line2_split[7], line2_split[5], line2_split[3], line2_split[4],"DEL",GT))
							continue

			if line2_split[1] == "DEL":
                           
			    if line2_split[11].find("PE") == -1 and line2_split[11].find("SR") != -1 and int(line2_split[7]) - int(line2_split[3]) < SR_DEL_THRESH:
                                continue
                            elif line2_split[11].find("SR") == -1 and line2_split[11].find("PE") != -1 and int(line2_split[7]) - int(line2_split[3]) < PE_DEL_THRESH:
                                continue
 
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

                        elif len(line2_split[1]) > 2 and line2_split[1][0:3] == "INS" and not (line2_split[1]== "INS_C" or line2_split[1] == "INS_C_I"):

			    # two lines for insertion as in bedpe format; bp 1 and bp3 were flanks of bp2 by convention in INS_C classification unless confirmed further
			    [bp1_s, bp1_e]= int(line2_split[3]), int(line2_split[4])
			    [bp2_s, bp2_e] = int(line2_split[6]),int(line2_split[7])
                            [bp3_s, bp3_e] = int(line2_split[9]),int(line2_split[10])
			    
			    if swap:
				temp2 = bp1_s
				bp1_e = bp3_e
				bp1_s = bp3_s
				bp3_s = temp2

				if bp3_s < bp2_s:
					bp2_s = bp3_s
					bp3_e = bp2_e
						
                            f16.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line2_split[5], bp2_s, bp3_e, line2_split[2], bp1_s, bp1_e, line2_split[1],GT))

                        else:
				f17.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],"BND", line2_split[1],GT))
 
	    	if line2_split[1] == "DEL_uc":

                            f13b.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1],GT))
	
	#print "PE good call ratio:", 1.0*nv_PEconf/nv_PEposs
        #print "SR good call ratio:", 1.0*nv_SRconf/nv_SRposs	
	
	#chrHash.clear()
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



                                        
                        
