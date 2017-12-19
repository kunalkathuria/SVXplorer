#!/usr/bin/env python

import sys

if __name__=="__main__":

	f=open(sys.argv[1],"r")
	for i in range(int(sys.argv[2]),int(sys.argv[2])+50):	
		supp_thresh = i
		f.seek(0)
		fw=open(sys.argv[1]+"minsupp."+str(i),"w")
		for line in f:
			if line[0] != "#":
				suppT=0
				try:
					pe=int(line.split(";PE=")[1].split(";")[0])
					suppT+= pe
				except:
					pass
				
				try:
					sr=int(line.split(";SR=")[1].split(";")[0])
					suppT+= sr

				except:
					pass

				if suppT >= supp_thresh:
                                                fw.write("%s" %line)
			else:
				fw.write("%s" %line)

		fw.close()
