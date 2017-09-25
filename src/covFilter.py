import sys
import pysam

if __name__== "__main__":

	CONC_BUFFER = 280 #ERRLIB
	DISC_BUFFER = 120
	RDL=101
	R_THRESH = 1.20
	BAMFILE= sys.argv[1]
	BEDFILE = sys.argv[2]
	f=open(BEDFILE,"r")
	samfile = pysam.AlignmentFile(BAMFILE, "rb" )	
	fw1=open("../results/text/poss_bad_dels.txt","w")
	fw2=open("../results/text/poss_good_dels.txt","w")
	bad_count = 0

	for k,line2 in enumerate(f):
		print k
		line2_split = line2.split()
#		chr_n = line2_split[0]
#		start = int(line2_split[1]) - DISC_BUFFER
#		stop = int(line2_split[1])
#		covLoc = 0
#		counter2 = 0
#		for pileupcolumn in samfile.pileup(chr_n, start, stop):
#				#print "DEL hash", pileupcolumn.pos, pileupcolumn.n
#				covLoc = covLoc + pileupcolumn.n
#				counter2+=1
#		covLoc1 = covLoc/(1.0*counter2)

		lmin = int(line2_split[-5])
		lmax = int(line2_split[-4])
		n_reads = int(line2_split[1])
		covLoc1 = (RDL*n_reads)/(1.0*(lmax-lmin))

		chr_n = line2_split[3]
		start = int(line2_split[4]) - CONC_BUFFER - 200
		stop = int(line2_split[4]) - CONC_BUFFER
		covLoc = 0
		covLoc2 = 0
		counter2 = 0

		if start > 0:
			for pileupcolumn in samfile.pileup(chr_n, start, stop):
				#print "DEL hash", pileupcolumn.pos, pileupcolumn.n
				covLoc = covLoc + pileupcolumn.n
				counter2+=1

		if counter2 > 0:
			covLoc2 = covLoc/(1.0*counter2)
			ratio = abs(covLoc1 - covLoc2)/(1.0*min(covLoc1, covLoc2))
		else:
			print "Zero counter:", line2, chr_n, start, stop
			ratio = 0


		if covLoc1 > 0 and covLoc2 > 0 and ratio > R_THRESH:

			bad_count+=1
			print "bad count is", bad_count
			#print ratio
			fw1.write("%s %s %s %s\n" %(line2[:-1], ratio, covLoc1, covLoc2))

		elif covLoc1 > 0 and covLoc2 > 0:

			#print "small ratio", line2, covLoc1, covLoc2, ratio
                        fw2.write("%s %s %s %s\n" %(line2[:-1], ratio, covLoc1, covLoc2))

	print k, bad_count, bad_count/(1.0*(k+1))
