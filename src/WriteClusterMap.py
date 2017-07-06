# Identify all discordant reads

# Notes: can include different methods (SR, PEM) and use this algorithm to resolve ambiguities

import sys

FILE = "../results/text/VariantMapInp.txt"
MIN_CLUSTER_SIZE = int(sys.argv[1])

if __name__ == "__main__":

    f1 = open(FILE,"r")
    f2 = open("../results/text/VariantMap.txt","w")
    
    prevVars = []	
    temp = f1.readline().split()
    ptemp1 = int(temp[1])
    prevVars.append(int(temp[0]))
    #new = 1
    #start = 0
    sim_counter = 1	

    #Ensure cluster supported by > 1 fragment
    for line in f1:
        temp = line.split()

        if(int(temp[1]) != ptemp1):
            if sim_counter > MIN_CLUSTER_SIZE:
	
		f2.write("%s " %ptemp1)
		for t in range(len(prevVars)):
                   f2.write("%s " %(prevVars[t]))	
                f2.write("\n")
	    sim_counter = 1
	    prevVars = []
            prevVars.append(temp[0])
        else:
	    sim_counter+=1		
            prevVars.append(temp[0])	   
        
	ptemp1 = int(temp[1])

    if sim_counter > MIN_CLUSTER_SIZE:
	f2.write("%s " %temp[1])

	for t in range(len(prevVars)):
           f2.write("%s " %(prevVars[t]))	
        f2.write("\n")
    f1.close()
    f2.close()
    
    
