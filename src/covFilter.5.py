import sys
import pysam
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

if __name__== "__main__":

	SBUFFER = 215
	RDL=102
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
		print "Cluster", 1+k
		line2_split = line2.split()

		bpstart = int(line2_split[4])

		chr_n = line2_split[3]
		stop = bpstart
		start = stop - SBUFFER

		n_frags_c = 0
		n_frags_d = 0 
		count_c = 0

		if start > 0 and line2_split[2] == "01":
			for read in samfile.fetch(chr_n, start, stop):
				#print "DEL hash", pileupcolumn.pos, pileupcolumn.n
				if read.is_proper_pair and (start <= read.reference_start <= stop) and not read.is_reverse and not read.is_secondary and not read.is_supplementary and stop < read.next_reference_start:
					n_frags_c+=1

                                #print "DEL hash", pileupcolumn.pos, pileupcolumn.n
                                if not read.is_proper_pair and (start <= read.reference_start <= stop) and not read.is_reverse and not read.is_secondary and not read.is_supplementary and stop < read.next_reference_start:
                                        n_frags_d+=1

				# don't use if conc reads exist in interval
				if read.is_proper_pair and read.is_reverse and (stop <= read.reference_start <= stop + 950):
					count_c+= 1	


		#n_frags_c = n_frags_c*.5
		if start > 0:
			print n_frags_c, n_frags_d, 1.0*n_frags_c/n_frags_d, count_c


