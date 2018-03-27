#kk7t
import sys

if __name__=="__main__":

	const = int(sys.argv[2])
	search=0
	f=open(sys.argv[1], "r")
	fw=open(str(sys.argv[1]+".bedpe"), "w")
	first = 0

	for line in f:
		line_split = line.split()
		if line_split[0][0] != "#":
			search = 1
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


	
		if (search and line_split[0][0] != "#"):
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], start - const/2, start + const/2, line_split[0], end - const/2, end+const/2))

	fw.close()	
