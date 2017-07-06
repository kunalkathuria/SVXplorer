#kk7t
import sys

if __name__=="__main__":

	shift = 1000
	search=0
	f=open(sys.argv[1], "r")
	fw=open(str(sys.argv[1]+".ed.vcf"), "w")

	for line in f:
		line_split = line.split()
		if line_split[0] == "#CHROM":
			search = 1
			continue
		if search:
			eflag = 0
			end = -1

			if len(line_split) > 7 and line_split[4]. find("<DUP:ISP>") != -1::
				eflag =line_split[7].find("CHR2")
				if eflag != -1:
					end = int(line_split[7].find("END2")])
				line_split[7] = line_split[7][:eflag] + line_split[7][end:]
				
				line_split[1] = str(int(line_split[1]) - shift)
				line_split[4] = "<DUP:TANDEM>"
				fw.write("%s\n" % ("\t".join(line_split)))

	fw.close()	
