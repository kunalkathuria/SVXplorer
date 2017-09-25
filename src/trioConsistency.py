import sys

if __name__=="__main__":

	#mother, father, child
	f1=open(sys.argv[1],"r")
	f2=open(sys.argv[2],"r")
	f3=open(sys.argv[3],"r")
	
	for line in f1:
		chrm = line.split()[0]
		startm = line.split()[1]

		
