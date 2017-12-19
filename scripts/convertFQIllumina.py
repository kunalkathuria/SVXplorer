import sys
FILE1 = sys.argv[1]
FILE2 = sys.argv[2]
#include gzip
import gzip

if __name__=="__main__":

	f1=gzip.open(FILE1, "r")
	f2=gzip.open(FILE2, "r")
	f3=open(str(sys.argv[3]+"/1_ed.fq"),"w")
	f4=open(str(sys.argv[3]+"/2_ed.fq"),"w")

	for line in f1:
		line_split = line.split()

		if line_split[0][0:4] == "@SRR" or line_split[0][0:4] == "+SRR":
			f3.write("%s\n" %(line_split[0]+"/1"))
		else:
			f3.write("%s" %line)

	for line in f2:
                line_split = line.split()

                if line_split[0][0:4] == "@SRR" or line_split[0][0:4] == "+SRR":
                        f4.write("%s\n" %(line_split[0]+"/2"))
                else:
                        f4.write("%s" %line)
