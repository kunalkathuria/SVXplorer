import sys
import pysam

if __name__== "__main__":

	CONC_BUFFER = 280 #ERRLIB
	RDL=101
	R_THRESH = .10
	BAMFILE= sys.argv[1]
	BEDFILE = sys.argv[2]
	f=open(BEDFILE,"r")
	samfile = pysam.AlignmentFile(BAMFILE, "rb" )	
	fw1=open("../results/text/poss_bad_dels.txt","w")
	fw2=open("../results/text/poss_good_dels.txt","w")
	bad_count = 0
	clist = []
	dlist = []

	for k,line2 in enumerate(f):
		print k
		line2_split = line2.split()

		lmin = int(line2_split[-5])
		lmax = int(line2_split[-4])
		n_frags_d = int(line2_split[1])

		chr_n = line2_split[3]
		start = lmin - RDL - (lmax - lmin + RDL)
		stop = lmin - RDL # there should be no discordant reads before lmin - RDL
		n_frags_c = 0

		if start > 0 and line2_split[2] == "01":
			for read in samfile.fetch(chr_n, start, stop):
				#print "DEL hash", pileupcolumn.pos, pileupcolumn.n
				if read.is_proper_pair and read.reference_start >= start and read.reference_end <= stop and read.is_read1 > 0 and not read.is_secondary and not read.is_supplementary:
					n_frags_c+=1

		clist.append(n_frags_c)
		dlist.append(n_frags_d)
		#n_frags_c = n_frags_c*.5
		if n_frags_c > 0:
			print n_frags_c, n_frags_d
			ratio = abs(n_frags_d-n_frags_c)/(1.0*min(n_frags_c, n_frags_d))
		else:
			#print "Zero counter:", line2, chr_n, start, stop
			ratio = 0


		if ratio > R_THRESH:

			bad_count+=1
			#print "bad count is", bad_count
			#print ratio
			fw1.write("%s\t%s\t%s\n" %(line2_split[3], line2_split[4], line2_split[9]))

		elif ratio != 0:

			#print "small ratio", ratio
			fw1.write("%s\t%s\t%s\n" %(line2_split[3], line2_split[4], line2_split[9]))

	print k, bad_count, bad_count/(1.0*(k+1))
