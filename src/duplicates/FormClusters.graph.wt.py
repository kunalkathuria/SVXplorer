import math
import time
import sys
import pysam
import numpy as np
import matplotlib as mp
import copy
import networkx as nx

counter =1

def READ_BAM_STATS(file1):
	
	stats = []	
	fp = open(file1,"r")
	for line in fp:
		stats.append(line)
	fp.close()
	return float(stats[0]), float(stats[1]), float(stats[2]), float(stats[5]), float(stats[6]), float(stats[7])

STAT_FILE = sys.argv[1]
MIN_CLUSTER_SIZE = int(sys.argv[2])
[RDL, MEAN_D, SIG_D, DISC_D, DIST_PEN, DIST_END] = READ_BAM_STATS(STAT_FILE)
SIG_MULT = int(sys.argv[3])
DISC_ENHANCER = float(sys.argv[9])# 1.67 default, wil multiply 3 sigma disc distance by a "safety" factor to not miss larger sampling IL clusters
DISC_ENHANCER2 = 1.0 #DISC_ENHANCER*(5.0/5.0) # safety factor for cluster margin (e.g. if normal distr, 3 sig for disc marking, 4 sig for cluster margin, 5 sig for cluster inclusion)
CLUSTER_D = MEAN_D + DISC_ENHANCER*DISC_D - 2*RDL
BP_MARGIN = int(sys.argv[4])
SIG_BOUND = int(sys.argv[5])
CLD_THRESH = float(sys.argv[6])
RDL_FACTOR = float(sys.argv[7])
CL_BREAK = int(sys.argv[8])
LIB_MULT = float(sys.argv[10]) #liberal distance multiplier to include reads in given cluster even if far from beginning of cluster; used if reads piling up close; cluster will be broken later based on density
BINDIST_HASH = {}
print RDL, MEAN_D, SIG_D, DISC_D, CLUSTER_D   
Weights = []
MINBPDIST = 100

def readDistHash():
	TOTAL_ENTRIES = 0
	f=open("../results/text/bindist.txt","r")
	for line in f:
		ls1 = int(line.split()[1])
		BINDIST_HASH[int(line.split()[0])]=ls1
		TOTAL_ENTRIES+= ls1
	print TOTAL_ENTRIES
	return TOTAL_ENTRIES

def calcEdgeWeight(f1lpos, f1rpos, f2lpos, f2rpos, totalent, c_type):

	BIN_SIZE = 10
	pen_atmargin = .75

	if c_type == "01" or c_type == "10":
		ildist = abs(BIN_SIZE*int((f2lpos - f1lpos + f1rpos - f2rpos)/(1.0*BIN_SIZE)))
	elif c_type == "00" or c_type == "11":
		ildist = abs(BIN_SIZE*int((f2lpos - f1lpos + f2rpos - f1rpos)/(1.0*BIN_SIZE)))
	
	distpen = 1 + (pen_atmargin - 1)*(abs(f2lpos-f1lpos)+2*RDL-DIST_PEN)*1.0/(DIST_END - DIST_PEN)

	if distpen > 1:
		distpen = 1
	if ildist in BINDIST_HASH and distpen > 0: #and f1lpos < f2rpos - MINBPDIST and f2lpos < f1rpos - MINBPDIST:
		weight = distpen*BINDIST_HASH[abs(ildist)]/(1.0*totalent)
	else:
		weight = 0

	#print "Weight:", weight
	return weight

def swap(a, b):
	
	temp = a
	a= b
	b = temp
	
def calculateMargin(item):

		l_orient = int(item.C_type[0])
		r_orient = int(item.C_type[1])
		
		cl_margin_l = MEAN_D + 2*DISC_ENHANCER*DISC_D - 2*RDL - (item.lmax - item.lmin)
		cl_margin_r = MEAN_D + 2*DISC_ENHANCER*DISC_D - 2*RDL - (item.rmax - item.rmin)

		cl_margin = min(cl_margin_l, cl_margin_r)
		cl_margin = int(math.ceil(cl_margin))

		if cl_margin <= BP_MARGIN:
			# use small margin
			cl_margin = BP_MARGIN

		ERROR_MARGIN = int(BP_MARGIN/2)
		# for SVSim like test where tolerance is only added to start of SV

		#$test-- remove later. Results v good w/SVSim tally on this.
		#cl_margin = ERROR_MARGIN

		if l_orient == 0:
                        item.l_end = item.l_bound + cl_margin
                        item.l_start = item.l_bound - ERROR_MARGIN

                else:
                        item.l_start = item.l_bound - cl_margin
                        item.l_end = item.l_bound + ERROR_MARGIN

		if r_orient == 0:
			item.r_end = item.r_bound + cl_margin
			item.r_start = item.r_bound - ERROR_MARGIN
		else:
			item.r_start = item.r_bound - cl_margin
			item.r_end = item.r_bound + ERROR_MARGIN

		if item.l_start < 0:
			item.l_start = 0
			item.l_end = max(item.l_start + 1, item.l_end)

		#if diff chr, right can be 0
		if item.r_start < 0:
                        item.r_start = 0
                        item.r_end = max(item.r_start + 1, item.r_end)

		if item.ltid == item.rtid: 
			if not l_orient and r_orient:

				item.l_end = min(item.l_end, item.r_end-1)
				item.r_start = max(item.r_start, item.l_start+1)

			elif not l_orient and not r_orient:

                                item.r_start = max(item.lmax+1, item.r_start)

			elif l_orient and r_orient:

				item.l_end = min(item.rmin-1, item.l_end)

# Class definitions may be local for small classes even though used in different modules
# See ReadDiscordant for variable details
class Cluster(object):

    # loosely used for individual fragments to keep record as well
    def __init__(self):
            
            self.l_bound=None
            self.r_bound=None
            
            self.count=1
            self.C_type=None 

            self.ltid = -1
            self.rtid = -1
	    self.lmin = -1
	    self.lmax = -1
	    self.rmin = -1
            self.rmax = -1

	    # final cluster locations	
	    self.l_start = -1 
            self.l_end = -1
            self.r_start = -1
            self.r_end = -1
	    self.clsmall = -1
	    self.fragNum = -1

	    self.used = 0
	    self.blockID = -1	
	
    def __str__(self):
        return "%s %s %s %s %s %s %s %s %s" % (self.count, self.C_type, self.ltid, self.l_start, self.l_end, self.rtid, self.r_start, self.r_end, self.clsmall)

def writeClusters(G, fragHash, fc, ft, f3, f4, preserveList):

    currentCliques = []

    for connected_comp in list(nx.connected_component_subgraphs(G)):
	preserveComp = 0
	if len(connected_comp) >= MIN_CLUSTER_SIZE:

	   currentCliquesC = []
  	   for node in connected_comp:
	      #keep those nodes/components intact that are connected to recently added nodes on reference line
	      if node in preserveList:
			preserveComp = 1
			break
	      if G.degree(node) >= MIN_CLUSTER_SIZE - 1:
		edgeC_T = 0
		maxavgW = 0
		for k,cl in enumerate(currentCliquesC):
			sumW = 0
			edgeC = 0
			for item in cl:
				if G.has_edge(node,item):
					sumW = sumW + G[node][item]['weight']
					edgeC+=1
					edgeC_T+=1
				if edgeC_T == G.degree(node):
					break
			if edgeC > 0:
				avgW = sumW/(1.0*edgeC)
				if avgW > maxavgW:
					maxavgW = avgW
					maxindex = k
			#elif avgW == maxavgW and len(cl) > < len(currentCliques[maxindex]):
				#maxindex = k
					
		if edgeC_T == 0:
			currentCliquesC.append([node])
		else:
			currentCliquesC[maxindex].append(node)
	   
	   currentCliques.extend(currentCliquesC)

	if len(preserveList) > 0 and len(connected_comp) < MIN_CLUSTER_SIZE:
		for nodeC in connected_comp:
			if nodeC in preserveList:
				preserveComp = 1
				break
	if not preserveComp:
                G.remove_nodes_from(connected_comp)


    #print "Max cliques:", max_clique_list2
    global counter
    for elem in currentCliques:

        elem0 = elem[0]
        ltid = fragHash[elem0].ltid
        rtid = fragHash[elem0].rtid
        lbp = fragHash[elem0].l_bound
        rbp = fragHash[elem0].r_bound
	lmin = fragHash[elem0].l_bound
	lmax = fragHash[elem0].l_bound + 1
	rmin = fragHash[elem0].r_bound
	rmax = fragHash[elem0].r_bound + 1
        fclist = []
        count = 0
        print "Clique size is:", len(elem)
        if len(elem) >= MIN_CLUSTER_SIZE:
        	for item in elem:
        	   if not fragHash[item].used:
                        count+=1
                        fclist.append(item)
                        fragHash[item].used = 1

                        if fragHash[item].C_type[0] == "0" and fragHash[item].l_bound > lbp:
	                        lbp = fragHash[item].l_bound
			elif fragHash[item].C_type[0] == "1" and fragHash[item].l_bound < lbp:
				lbp = fragHash[item].l_bound
                        if fragHash[item].C_type[1] == "0" and fragHash[item].r_bound > rbp:
                                rbp = fragHash[item].r_bound
			elif fragHash[item].C_type[1] == "1" and fragHash[item].r_bound < rbp:
                                rbp = fragHash[item].r_bound

			if fragHash[item].l_bound < lmin:
				lmin = fragHash[item].l_bound
			elif fragHash[item].l_bound > lmax:
				lmax = fragHash[item].l_bound

			if fragHash[item].r_bound < rmin:
                                rmin = fragHash[item].r_bound
                        elif fragHash[item].r_bound > rmax:
                                rmax = fragHash[item].r_bound

                        #Compute and write sum of total edge weights for this node to all other nodes in clique
                        #fc.write

        	if count >= MIN_CLUSTER_SIZE:
                        newCl = Cluster()
                        newCl.count = count
                        newCl.ltid = ltid
                        newCl.rtid = rtid
                        newCl.l_bound = lbp
                        newCl.r_bound = rbp
                        newCl.C_type = fragHash[elem0].C_type
                        newCl.clsmall = fragHash[elem0].clsmall
			newCl.lmin = lmin
			newCl.rmin = rmin
			newCl.lmax = lmax
			newCl.rmax = rmax
			#print "Cluster is", newCl.lmin, newCl.lmax, newCl.rmin, newCl.rmax
                        calculateMargin(newCl)
			#print "Cluster now is:", newCl
                        #print "Cluster clique formed:", counter
                        fc.write("%s %s\n" %("@Cluster"+str(counter),newCl))
                        f3.write("%s %s\n" %(counter, newCl))
                        f4.write("%s" %counter)
                        for item in fclist:
                                #fc.write("%s %s %s\n" %(item, fragHash[item].l_bound, fragHash[item].r_bound))
                                f4.write(" %s" %item)
			f4.write("\n")
                        counter+=1

class Read(object):

	def __init__(self):
		name = -1
		l_read_bound = -1
		r_read_bound = -1
	def __str__(self):
		return "%s %s" %(self.l_read_bound, self.r_read_bound)

if __name__== "__main__":
    
    f1=open("../results/text/All_Discords_P_S.txt","r")
    #f2=open("../results/text/All_Discords_I_S.txt","r")
    f3=open("../results/text/All_Clusters.txt","w")
    f4=open("../results/text/VariantMap.txt","w")
    f5=open("../results/text/Time_FC_loop.txt","w")
    f6=open("../results/text/Time_FC_line.txt","w")
    fc = open("../results/text/cluster_cliques.txt","w")
    ft=open("../results/text/cc_size.txt","w") 

    print "Function start"
    
    FList = []
    newBlockS = 0

    totalent = readDistHash()
    fragHash = {}
    FG = nx.Graph()
    print "Cluster D:", CLUSTER_D
	
    for line_num,line in enumerate(f1):

        #start_fl = time.clock()
	if line_num % 100 == 0:
		print line_num	
        parsed = line.split()

        temp = Cluster()
        temp.l_bound = int(parsed[2])
        temp.r_bound = int(parsed[4])
	temp.lmin = temp.l_bound
	temp.lmax = temp.lmin + 1
	temp.rmin = temp.r_bound
	temp.rmax = temp.rmin + 1
        temp.C_type = parsed[5]
        temp.ltid = parsed[1]
        temp.rtid = parsed[3]
	temp.clsmall = parsed[-2]
        currentFrag = parsed[0]
	temp.fragNum = currentFrag
        Claimed = 0

        l_orient = int(temp.C_type[0])
        r_orient = int(temp.C_type[1])

	#ignore artefact read seen often
	if temp.l_bound == temp.r_bound and (temp.C_type == "00" or temp.C_type == "11"):
		continue
	
	#criss crossing read from small TD may appear like this
	if temp.l_bound > temp.r_bound and temp.ltid == temp.rtid and temp.C_type == "01":
		temp2 = l_orient
		l_orient = r_orient
		r_orient = temp2

		temp2 = temp.l_bound
                temp.l_bound = temp.r_bound
                temp.r_bound = temp2

		temp2 = temp.lmin
                temp.lmin =  temp.rmin
                temp.rmin = temp2

		temp2 = temp.lmax
                temp.lmax =  temp.rmax
                temp.rmax = temp2

		temp2 = l_orient
                l_orient = r_orient
                r_orient = temp2
		temp.C_type = "10"

        ld = len(FList)
        #start_for = time.clock()
	#print line_num
	EW_THRESH = 0

	FG.add_node(temp.fragNum)
	fragHash[temp.fragNum] = temp
	FList.append(temp)
		
        for x in range(ld):
	    #print "Member:", FList[x]
 	    if FList[x].ltid == temp.ltid and FList[x].rtid == temp.rtid and FList[x].C_type == temp.C_type and FList[x].clsmall == temp.clsmall:

		    #compare 2 frags
		    f1lpos = FList[x].l_bound
		    f1rpos = FList[x].r_bound	    
		    f2lpos = temp.l_bound
		    f2rpos = temp.r_bound
		    #print x, f1lpos, f1rpos, f2lpos, f2rpos
		    EW = calcEdgeWeight(f1lpos, f1rpos, f2lpos, f2rpos, totalent, temp.C_type)
		    #print "Weight is:", EW
		    #Weights.append(EW)


		    if EW > EW_THRESH:
			#print "Adding connection"
			if FList[x].fragNum != temp.fragNum:
				FG.add_edge(FList[x].fragNum, temp.fragNum, weight = EW)
		
	# refresh list	
        if len(FList) > 1:
		
	 	if temp.ltid != FList[-2].ltid:
			
			#print "New TID: Processing ltid", temp.ltid, line, "Size of list:", len(DList)
			newBlockS = 0
                        FList = [FList[-1]]
			writeClusters(FG, fragHash, fc, ft, f3, f4,[])
			FG.clear()
			FG.add_node(FList[0].fragNum)
			fragHash = {}
			fragHash[FList[0].fragNum] = FList[0]

		elif temp.ltid == FList[newBlockS].ltid and (temp.l_bound - FList[newBlockS].l_bound) > 2*CLUSTER_D:

                    #FList = FList[newBlockS:]
		    nodeList = {}
                    for frag in FList[newBlockS:]:
                        nodeList[frag.fragNum] = 1
	  	    writeClusters(FG, fragHash, fc, ft, f3, f4, nodeList)
		    del FList[0:newBlockS]
                    newBlockS = len(FList)-1
	
    writeClusters(FG, fragHash, fc, ft, f3, f4,[])

    FG.clear()
    fc.close()
    ft.close() 
    f3.close()
    f4.close()            
   
       
