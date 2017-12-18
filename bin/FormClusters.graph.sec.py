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
LIB_MULT = float(sys.argv[10]) #obsolete -- liberal distance multiplier to include reads in given cluster even if far from beginning of cluster; used if reads piling up close; cluster will be broken later based on density
WORK_DIR = sys.argv[12]
BINDIST_HASH = {}
print RDL, MEAN_D, SIG_D, DISC_D, CLUSTER_D   
MINBPDIST = 100

def readDistHash():
	TOTAL_ENTRIES = 0
	f=open(WORK_DIR+"/bindist.txt","r")
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
	
	distpen = 1 - (abs(f2lpos-f1lpos)+2*RDL-DIST_END)*1.0/(DIST_PEN - DIST_END)

	if distpen > 1:
		distpen = 1
	if ildist in BINDIST_HASH and distpen > 0 and (c_type == "10" or (f1lpos < f2rpos and f2lpos < f1rpos)):
		#strictly, the denominator should be ~.5*totalent^2 but it is an overall constant and thus won't matter
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

		MIN_PE_SIZE = 2
		SMALL_CL_MARGIN = 50
		
		l_orient = int(item.C_type[0])
		r_orient = int(item.C_type[1])
		
		cl_margin_l = MEAN_D + DISC_ENHANCER*DISC_D - 2*RDL - (item.lmax - item.lmin)
		cl_margin_r = MEAN_D + DISC_ENHANCER*DISC_D - 2*RDL - (item.rmax - item.rmin)

		cl_margin = min(cl_margin_l, cl_margin_r)
		cl_margin = int(math.ceil(cl_margin))

		#SR support for low support PE clusters may be significant primarily if lying close
		if item.count <= MIN_PE_SIZE:
			cl_margin = min(SMALL_CL_MARGIN,cl_margin)

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
	
    def __str__(self):
        return "%s %s %s %s %s %s %s %s %s" % (self.count, self.C_type, self.ltid, self.l_start, self.l_end, self.rtid, self.r_start, self.r_end, self.clsmall)

def writeClusters(G, fragHash, fc, ft, f3, f4, preserveList):

    max_clique_list = []

    for connected_comp in list(nx.connected_component_subgraphs(G)):
	preserveComp = 0
	#connected_comp_list = list(connected_comp)
	#eligible_edges = [(from_node,to_node,edge_attributes) for from_node,to_node,edge_attributes in my_network.edges(data=True) if edge_attributes['Weight'] > weightT]

        if len(connected_comp) >= MIN_CLUSTER_SIZE:
		for node in connected_comp:
			if node in preserveList:
                                preserveComp = 1
                                break
		#ft.write("%s %s %s\n" %(len(connected_comp),connected_comp.number_of_nodes(), connected_comp.number_of_edges()))
		if not preserveComp:
			max_clique_list.extend(list(nx.find_cliques(connected_comp)))

    	#print "Max cliques:", max_clique_list2
    	elif len(preserveList) > 0 and len(connected_comp) < MIN_CLUSTER_SIZE:
                for nodeC in connected_comp:
                        if nodeC in preserveList:
                                preserveComp = 1
                                break
        if not preserveComp:
                G.remove_nodes_from(connected_comp)

    global counter
    max_clique_list2 = sorted(max_clique_list, key=len, reverse=True)
    for k,elem in enumerate(max_clique_list2):
	if len(elem) < MIN_CLUSTER_SIZE:
		break
    if len(max_clique_list2) > 0 and k != len(max_clique_list2) - 1:
	    del max_clique_list2[k:]

    for elem in max_clique_list2:

	pickedFrags = {}
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
        #print "Clique size is:", len(elem)
        if len(elem) >= MIN_CLUSTER_SIZE:
        	for item in elem:

		   #do not put diff almts of same fragment in same cluster
		   uspos = item.find("_")
		   if uspos != -1:
			item_s = item[:item.find("_")]

		   else:
			item_s = item


		   if item_s not in pickedFrags and not fragHash[item].used:

			pickedFrags[item_s] = 1
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
				item_s = item
				uscore = item.find("_")
				if len(item) > 2 and uscore != -1:
					#print "Found duplicated", item
					item_s = item[:uscore]
				
                                f4.write(" %s" %item_s)
			f4.write("\n")
                        counter+=1

class Read(object):

	def __init__(self):
		name = -1
		l_read_bound = -1
		r_read_bound = -1
	def __str__(self):
		return "%s %s" %(self.l_read_bound, self.r_read_bound)

def doSubsample(fp):
	subsample = 0
	#subsample routine
	

	return subsample

if __name__== "__main__":
    
    f1=open(WORK_DIR+"/All_Discords_P_S.txt","r")
    #f2=open(WORK_DIR+"/All_Discords_I_S.txt","r")
    f3=open(WORK_DIR+"/All_Clusters.txt","w")
    f4=open(WORK_DIR+"/VariantMap.txt","w")
    f5=open(WORK_DIR+"/Time_FC_loop.txt","w")
    f6=open(WORK_DIR+"/Time_FC_line.txt","w")
    fc = open(WORK_DIR+"/cluster_cliques.txt","w")
    ft=open(WORK_DIR+"/cc_size.txt","w") 

    print "Function start"
    
    FList = []
    newBlockS = 0

    totalent = readDistHash()
    fragHash = {}
    edgeCounter = {}
    FG = nx.Graph()
    print "Cluster D:", CLUSTER_D
    Weights = []
    nocalcW = 1
    wtCalcT = 30000
    WP_LOW = 0
    WP_HIGH = .05
    SD_IL_HIGH = 70
    SD_IL_LOW = 15
    #if sigma_IL is high, chances of alignments starying into neighboring cluster/clique is higher
    WEIGHT_THRESHP = WP_LOW + (SIG_D-SD_IL_LOW)*1.0*(WP_HIGH-WP_LOW)/(SD_IL_HIGH-SD_IL_LOW) #.05
    if WEIGHT_THRESHP > WP_HIGH:
	WEIGHT_THRESHP = WP_HIGH
    elif WEIGHT_THRESHP < WP_LOW:
	WEIGHT_THRESHP = WP_LOW

    edgeStore = []
    EW_THRESH = -1

    BLOCK_GAP_THRESH = 25
    BLOCK_THRESH = 3*MIN_CLUSTER_SIZE
    block_lbound =-1000
    block_rbound = -1000
    N_READS_BLOCK = 0
    BLOCK_HASH = {}
    subsample = doSubsample(f1)
	
    for line_num,line in enumerate(f1):

        #start_fl = time.clock()
	if line_num % 100 == 0:
		print "Done with fragment", line_num	
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

	#subsample if alignments v close, for speed	
	if subsample:
		process = 0
		found = 0
		for z in range(temp.l_bound - BLOCK_GAP_THRESH, temp.l_bound + BLOCK_GAP_THRESH):
			for q in range(temp.r_bound - BLOCK_GAP_THRESH, temp.r_bound + BLOCK_GAP_THRESH):
				if (z,q) in BLOCK_HASH and BLOCK_HASH[(z,q)] > BLOCK_THRESH:
					found = 1
					BLOCK_HASH[(z,q)]+=1
					break
				elif (z,q) in BLOCK_HASH:
					found = 1
					process = 1
					BLOCK_HASH[(z,q)]+=1
					break

			if found:
				break

		if found and not process:
			#print "skipping", temp.fragNum, temp.l_bound, temp.r_bound
			continue
		elif not found:
			BLOCK_HASH[(temp.l_bound, temp.r_bound)] = 1
			

	FG.add_node(temp.fragNum)

	if temp.fragNum not in fragHash:
		fragHash[temp.fragNum] = temp
		edgeCounter[temp.fragNum] = 1
	else: 
		edgeCounter[temp.fragNum]+=1
		temp.fragNum=str(temp.fragNum) + "_" + str(edgeCounter[temp.fragNum])
		fragHash[temp.fragNum] = temp

	FList.append(temp)
	#print temp.ltid, temp.l_bound, temp.rtid, temp.r_bound, len(FList)
	#print "Appending" temp
	#for y in FList:
		#print y.fragNum, y.l_bound, y.r_bound
	if line_num == 0:
		t1 = time.time()

        for x in range(ld):
	    #print "Member:", FList[x]
 	    if FList[x].ltid == temp.ltid and FList[x].rtid == temp.rtid and FList[x].C_type == temp.C_type and FList[x].clsmall == temp.clsmall:
		    #print temp.l_bound, temp.r_bound, temp.fragNum, FList[x].l_bound, FList[x].r_bound, FList[x].fragNum
		    #print "Comparing 2 frags", temp, FList[x], line
		    #compare 2 frags
		    f1lpos = FList[x].l_bound
		    f1rpos = FList[x].r_bound	    
		    f2lpos = temp.l_bound
		    f2rpos = temp.r_bound
		    #print x, f1lpos, f1rpos, f2lpos, f2rpos
		    EW = calcEdgeWeight(f1lpos, f1rpos, f2lpos, f2rpos, totalent, temp.C_type)
		    #print "Weight is:", EW
		    if len(Weights) < wtCalcT and EW > 0:
		    	Weights.append(EW)
			edgeStore.append([FList[x].fragNum, temp.fragNum, EW])
			#print EW, len(Weights)
		    elif len(Weights) >= wtCalcT and nocalcW:
			Weights_S = sorted(Weights)
			indexW = int(WEIGHT_THRESHP*len(Weights))
			EW_THRESH = Weights_S[indexW]
			nocalcW = 0
			print "Edge weight threshold calculated to be", EW_THRESH, "at percentile:", WEIGHT_THRESHP
			for edgeS in edgeStore:
				#print edgeS
				EW = edgeS[2]
				frag1 = edgeS[0]
				frag2 = edgeS[1]
				if EW > EW_THRESH and edgeS[0] != edgeS[1]:
	                                FG.add_edge(frag1, frag2)
			edgeStore = []
		    #only add if calculation threshold to get EW_THRESH has been reached, else determine later
		    if not nocalcW and EW > EW_THRESH:
			#print "Adding connection", EW, EW_THRESH, FList[x].fragNum, temp.fragNum
			if FList[x].fragNum != temp.fragNum:
				FG.add_edge(FList[x].fragNum, temp.fragNum)

	if line_num == 1000:
		print "Loop end time (10000 runs):", -t1 +time.time(), len(FList)	
	# refresh list
        if len(FList) > 1:
		
	 	if temp.ltid != FList[-2].ltid or temp.rtid != FList[-2].rtid:
			
			print "New TID. Processing ltid: rtid: Size of list:", temp.ltid, temp.rtid, len(FList)
			newBlockS = 0
                        FList = [FList[-1]]
			if not nocalcW:
				writeClusters(FG, fragHash, fc, ft, f3, f4,[])
				FG.clear()
				FG.add_node(FList[0].fragNum)
				fragHash = {}
				fragHash[FList[0].fragNum] = FList[0]

		elif temp.ltid == FList[newBlockS].ltid and (temp.l_bound - FList[newBlockS].l_bound) > 2*CLUSTER_D:

                    #FList = FList[newBlockS:]
		    #print "Refresh 2"	
		    if not nocalcW:
		    	nodeList = {}
                    	for frag in FList[newBlockS:]:
                        	nodeList[frag.fragNum] = 1
                    	writeClusters(FG, fragHash, fc, ft, f3, f4, nodeList)
		    del FList[0:newBlockS]
                    newBlockS = len(FList)-1

    #write remaining graph entries
    for edgeS in edgeStore:
		EW = edgeS[2]
		frag1 = edgeS[0]
		frag2 = edgeS[1]
		if EW > 0 and edgeS[0] != edgeS[1]:
			FG.add_edge(frag1, frag2)

    writeClusters(FG, fragHash, fc, ft, f3, f4, [])

    fc.close()
    ft.close() 
    f3.close()
    f4.close()            
   
       
