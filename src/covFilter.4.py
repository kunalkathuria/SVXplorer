import sys
import pysam
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

if __name__== "__main__":

	MIN_CS = 2
	SBUFFER = 0
	RLENGTH = 100
	RDL=101
	R_THRESH = .35
	SPLITM = .65*RDL
	BAMFILE= sys.argv[1]
	BEDFILE = sys.argv[2]
	f=open(BEDFILE,"r")
	samfile = pysam.AlignmentFile(BAMFILE, "rb" )	
	fw1=open("../results/text/poss_bad_dels.txt","w")
	fw2=open("../results/text/poss_good_dels.txt","w")
	fwc=open("../results/text/nconc_frags.txt","w")
	fwd=open("../results/text/ndisc_frags.txt","w")
	bad_count = 0
	clist = []
	dlist = []

	for k,line2 in enumerate(f):
		print k
		line2_split = line2.split()

		lmin = int(line2_split[-5])
		lmax = int(line2_split[-4])
		#n_frags_d = int(line2_split[1])

		chr_n = line2_split[3]
		stop = lmin - RDL - SBUFFER# there should be no discordant reads before lmin - RDL
		stopd = lmax
		if lmax - lmin + RDL >= RLENGTH:
			start = stop - RLENGTH
			startd = stopd - RLENGTH
		else:
			start = -1

		n_frags_c = 0
		n_frags_d = 0 
		count_c = 0
		if start > 0 and line2_split[2] == "01":
			for read in samfile.fetch(chr_n, start, stop):
				#print "DEL hash", pileupcolumn.pos, pileupcolumn.n
				if read.is_proper_pair and (start <= read.reference_start <= stop - SPLITM) and read.template_length > 0 and not read.is_secondary and not read.is_supplementary:
					n_frags_c+=1

			for read in samfile.fetch(chr_n, startd, stopd):
                                #print "DEL hash", pileupcolumn.pos, pileupcolumn.n
                                if not read.is_proper_pair and (startd <= read.reference_start <= stopd - SPLITM) and read.template_length > 0 and not read.is_secondary and not read.is_supplementary:
                                        n_frags_d+=1

				# don't use if conc reads exist in interval
				if read.is_proper_pair and startd <= read.reference_start:
					count_c+= 1	


		#n_frags_c = n_frags_c*.5
		if start > 0 and n_frags_c > MIN_CS and n_frags_d > MIN_CS and count_c <= 5:
			print n_frags_c, n_frags_d
			ratio = abs(n_frags_d-n_frags_c)/(1.0*min(n_frags_c, n_frags_d))
			clist.append(n_frags_c)
                        dlist.append(n_frags_d)
			fwc.write("%s\n" %n_frags_c)
			fwd.write("%s\n" %n_frags_d)
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
			fw2.write("%s\t%s\t%s\n" %(line2_split[3], line2_split[4], line2_split[9]))

	n, bins, patches = plt.hist(clist,10,normed=1,alpha = .5, label = 'c')
	n2, bins2, patches2 = plt.hist(dlist,10,normed=1,alpha=.5, label = 'd')

	plt.xlabel('Number of Fragments in Window')
        plt.ylabel('Frequency')
        plt.axis([0, 100, 0, .1])
        plt.show()
	print k, bad_count, bad_count/(1.0*(k+1))
