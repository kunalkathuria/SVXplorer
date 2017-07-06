#kk7t
import sys

if __name__=="__main__":

	CHROM="20"
	NINS=2500
	search=0
	fw=open(str("ins2.vcf"), "w")

	for counter in range(1,NINS+1):
		start = counter*9000
		end  = start + 1000 

		#half copy-paste, half TD
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(CHROM, str(start+5002), str(counter), ".", "<DUP:ISP>", ".", "PASS","SVTYPE=DUP;SVLEN="+str(end-start) + ";CHR2="+ CHROM+ ";POS2=" + str(start+4000)+";END2=" + str(start+5000)))

	fw.close()
