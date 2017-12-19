import pysam

pileThresh = 20
BAM = sys.argv[1]

def addDelTag(file1):

   bamfile = pysam.AlignmentFile(file1,"rb")

   f1 = open("../results/text/All_Variants.txt","r")
   f2 = open("../results/text/All_Variants_DC.txt", "w")
   f3 = open("../results/text/bam_stats.txt", "r")

   temp = f3.readline()
   MEAN_IL = float(f3.readline())

   for line in f1:

	ls = line.split()
	chr_n = ls[2]
	tag = 0
	counter = 0

	if ls[1] == "DEL":

		start = int(ls[4])
		stop = int(ls[6])

		for read in bamfile.fetch(chr_n, start, start + MEAN_IL):
		
			if read.is_proper_pair and read.is_reverse and read.reference_start - start < MEAN_IL:
			
				counter+= 1	
			
			if counter > pileThresh:
				tag = 1
				break

		for read in bamfile.fetch(chr_n, stop - MEAN_IL, stop):

                        if read.is_proper_pair and not read.is_reverse and stop - read.reference_end < MEAN_IL:

                                counter+= 1

                        if counter > pileThresh:
                                tag = 1
                                break

	if tag:
	
		ls.append("TBD")
	
	varLine = " ".join(ls)
	f2.write("%s\n" %varLine)
		
		
   bamfile.close()

if __name__== "__main__":					

	addDelTag(BAM)

