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
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(CHROM, str(start), str(counter*3-2), ".", "<DEL:TRA>", ".", "PASS","SVTYPE=DEL;SVLEN=-"+str(end-start) + ";TRAID=" + str(counter)))
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(CHROM, str(start+2000), str(counter*3-1), ".", "<DUP:TRA>", ".", "PASS","SVTYPE=DUP;SVLEN="+str(end-start) + ";CHR2="+ CHROM+ ";POS2=" + str(start)+";END2=" + str(end) +";TRAID=" + str(counter)))

		#half copy-paste, half TD
		if counter % 2 == 0:
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(CHROM, str(start+6000), str(counter*3), ".", "<DUP:ISP>", ".", "PASS","SVTYPE=DUP;SVLEN="+str(end-start) + ";CHR2="+ CHROM+ ";POS2=" + str(start+4000)+";END2=" + str(start+5000)))
		else:
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(CHROM, str(start+6000), str(counter*3), ".", "<DUP:TANDEM>", ".", "PASS","SVTYPE=DUP;SVLEN="+str(end-start)))
	fw.close()
