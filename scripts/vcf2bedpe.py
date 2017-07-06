import sys

if __name__=="__main__":

	search=0
	f=open(sys.argv[1], "r")
	fw=open(str(sys.argv[1]+".bedpe"), "w")

	for line in f:
		line_split = line.split()
		if line_split[0] == "#CHROM":
			search = 1
			continue
		if search:
			chrom= line_split[0]
			if len(chrom) > 2 and chrom[0:3] == "chr":
				chrom = chrom[3:]
			#svlen= line_split[7](line_split[7].find("SVLEN") + 6)
			short_str = line_split[7][line_split[7].find(";CIPOS"):]
			cipos = short_str[7:short_str.find(",")]
			ciend=0
			if line_split[4] == "<DEL>":
				short_str = line_split[7][line_split[7].find(";CIEND"):]
	                        ciend = short_str[7:short_str.find(",")]
				short_str = line_split[7][line_split[7].find(";END"):]
				start2 = int(short_str[5: 1 + short_str[1:].find(";")]) - abs(int(ciend))
				stop2 = int(short_str[5: 1 + short_str[1:].find(";")]) + abs(int(ciend))
			else:
				start2 = -1
				stop2 = -1

			start1 = int(line_split[1]) - abs(int(cipos))
			stop1 = int(line_split[1]) + abs(int(cipos))

			#print line, cipos, ciend, start1, stop1, start2, stop2
			if line_split[4] == "<DEL>":

				fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chrom, start1, stop1, chrom, start2, stop2, line_split[4], line_split[6]))
			
			elif line_split[4] == "<INS>":

				fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chrom, start2, stop2, chrom, start1, stop1, line_split[4], line_split[6]))
			
