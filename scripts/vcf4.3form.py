#kk7t
import sys

if __name__=="__main__":

	search=0
	f=open(sys.argv[1], "r")
	fw=open(str(sys.argv[1]+".ed"), "w")

	for line in f:
		line_split = line.split()
		if line_split[0] == "#CHROM":
			search = 1
			continue
		if search:
			eflag = 0
			end = -1
			start = -1

			if len(line_split) > 7:
				eflag =line_split[7].find(";END")
				if eflag != -1:
					short_str = line_split[7][eflag:]
					end = int(short_str[5: 1 + short_str[1:].find(";")])
				start = int(line_split[1])

			svlen = end - start

			if len(line_split) > 7:

				if line_split[4]. find("<DEL>") != -1:
					svlen = -1*svlen
				line_split[7] = "SVLEN=" + str(svlen) + ";" + line_split[7]
				fw.write("%s\n" % ("\t".join(line_split)))

	fw.close()	
