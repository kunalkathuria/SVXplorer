import sys

if __name__=="__main__":

	f=open(sys.argv[1],"r")
	for i in range(int(sys.argv[2]),int(sys.argv[2])+50):	
		supp_thresh = i
		f.seek(0)
		fw=open(sys.argv[1]+".minsupp."+str(i),"w")
		for line in f:
			if line[0] != "#":
				suppT=0
				try:
					prindex, srindex = -1,-1
					line_split = line.split("\t")[8].split(":")
					#print line_split
					for index, item in enumerate(line_split):
						if line_split[index] == "PR":
							prindex = index
						elif line_split[index] == "SR":
							srindex = index

					if prindex != -1 and srindex != -1:
						suppT+= int(line.split()[9].split(":")[prindex].split(",")[1]) + int(line.split()[9].split(":")[srindex].split(",")[1])
					elif prindex != -1:
						suppT+= int(line.split()[9].split(":")[prindex].split(",")[1])

					elif srindex != -1:
						suppT+= int(line.split()[9].split(":")[srindex].split(",")[1])

					if srindex != -1:
						print line, srindex, prindex, suppT
				except:
					pass
				
	
				if suppT >= supp_thresh:
                                	fw.write("%s" %line)
			else:
				fw.write("%s" %line)

		fw.close()
