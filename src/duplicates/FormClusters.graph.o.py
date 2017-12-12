import math
import time
import sys
import pysam
import numpy as np
import matplotlib as mp
import copy
import networkx as nx

def READ_BAM_STATS(file1):
	
	stats = []	
	fp = open(file1,"r")
	for line in fp:
		stats.append(line)
	fp.close()
	return float(stats[0]), float(stats[1]), float(stats[2]), float(stats[5]), float(stats[6])

STAT_FILE = sys.argv[1]
MIN_CLUSTER_SIZE = int(sys.argv[2])
[RDL, MEAN_D, SIG_D, DISC_D, DIST_PEN] = READ_BAM_STATS(STAT_FILE)
SIG_MULT = int(sys.argv[3])
DISC_ENHANCER = float(sys.argv[9])# 1.67 default, wil multiply 3 sigma disc distance by a "safety" factor to not miss larger sampling IL clusters
DISC_ENHANCER2 = 1.0 #DISC_ENHANCER*(5.0/5.0) # safety factor for cluster margin (e.g. if normal distr, 3 sig for disc marking, 4 sig for cluster margin, 5 sig for cluster inclusion)
DISC_ENHANCER3 = float(sys.argv[11])
#DISC_D_O = DISC_D
DISC_D = DISC_D*DISC_ENHANCER3
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

def readDistHash():
	TOTAL_ENTRIES = 0
	f=open("../results/text/bindist.txt","r")
	for line in f:
		ls1 = int(line.split()[1])
		BINDIST_HASH[int(line.split()[0])]=ls1
		TOTAL_ENTRIES+= ls1
	print TOTAL_ENTRIES
	return TOTAL_ENTRIES

def calcEdgeWeight(f1lpos, f1rpos, f2lpos, f2rpos, totalent):

	BIN_SIZE = 10
	pen_atmargin = .75
	#FR
	ildist = abs(BIN_SIZE*int((f2lpos - f1lpos + f1rpos - f2rpos)/(1.0*BIN_SIZE)))
	distpen = 1 + (pen_atmargin - 1)*(abs(f2lpos-f1lpos)-DIST_PEN)*1.0/(DISC_D - DIST_PEN)
	if ildist in BINDIST_HASH and distpen > 0:
		weight = distpen*BINDIST_HASH[abs(ildist)]/(1.0*totalent)
	else:
		weight = 0
	#print "Weight:", weight
	return weight

def swap(a, b):
	
	temp = a
	a= b
	b = temp
	
def CompareCluster(map_type, DiscordCluster, Cl2):

   iPos = Cl2.l_bound
   lPos = Cl2.l_bound
   rPos = Cl2.r_bound
	
   val = 100000 # some large value
  
   #print "Fxn 1"
 
   if Cl2.C_type == DiscordCluster.C_type and Cl2.ltid == DiscordCluster.ltid and Cl2.rtid == DiscordCluster.rtid and (not DiscordCluster.C_type == "01" or DiscordCluster.clsmall == Cl2.clsmall):# by itself last one this did not work well due to small TD clusters involving partial/split PE alignments
 
     if map_type == 0 and (lPos - DiscordCluster.l_bound_O < CLUSTER_D or (abs(lPos-DiscordCluster.l_bound) < RDL_FACTOR*RDL and lPos - DiscordCluster.l_bound_O < MEAN_D - 2*RDL + LIB_MULT*DISC_D)):

        #print "Fxn 2"
        if DiscordCluster.C_type=="01" or DiscordCluster.C_type=="10":
	    val = abs(lPos - DiscordCluster.l_bound_O - (rPos - DiscordCluster.r_bound_O))	
          
        elif DiscordCluster.C_type =="00" or DiscordCluster.C_type == "11":
            val = abs(lPos - DiscordCluster.l_bound_O - (DiscordCluster.r_bound_O - rPos))
       
	# if latter condition not satisfied, we will see a new cluster. Equal sign to be safe, e.g for small insertion location cluster.
	#print rPos, lPos, DiscordCluster.l_bound, DiscordCluster.r_bound, rPos >= DiscordCluster.l_bound and lPos <= DiscordCluster.r_bound
	if val < 2*DISC_D and (Cl2.ltid != Cl2.rtid or (rPos >= DiscordCluster.l_bound and lPos <= DiscordCluster.r_bound_O)):
	    #print "Fxn 3"
            return 1
            
     if map_type == 1 and (iPos - DiscordCluster.l_bound_O) < CLUSTER_D:
        return 1

   return 0

def ChangeBound(Cl1, Cl2, typeC):

     	[l_orient, r_orient, indivOrient] = [int(Cl2.C_type[0]), int(Cl2.C_type[1]), int(Cl2.C_type[0])]
	
	if r_orient != 2:
		map_type = 0
	else:
		map_type = 1

        # Can add different conditions for different cluster categories later. Right now, all fall into type 1.
	if typeC == 1:
			
        	if map_type == 0:
		
			# The first condition is redundant currently, but safe	
			if Cl2.l_bound < Cl1.lmin:
				Cl1.lmin = Cl2.l_bound
			elif Cl2.l_bound > Cl1.lmax:
				Cl1.lmax = Cl2.l_bound
			
			if Cl2.r_bound < Cl1.rmin:
				Cl1.rmin = Cl2.r_bound
			elif Cl2.r_bound > Cl1.rmax:
				Cl1.rmax = Cl2.r_bound

			# $ clean up and combine with above? can't, need track of min and max in all cases.
            		if (not l_orient) and Cl2.l_bound > Cl1.l_bound and Cl2.l_bound!=-1:
                		
                		Cl1.l_bound = Cl2.l_bound
                		                                   
            		elif l_orient and Cl2.l_bound < Cl1.l_bound and Cl2.l_bound!=-1:
                		
                		Cl1.l_bound = Cl2.l_bound
                		
            		if (not r_orient) and Cl2.r_bound > Cl1.r_bound and Cl2.r_bound!=-1:
                	
                		Cl1.r_bound = Cl2.r_bound
                		
            		elif r_orient and Cl2.r_bound < Cl1.r_bound and Cl2.r_bound!=-1:
                		
                		Cl1.r_bound = Cl2.r_bound
                		

        	if map_type == 1:
           
			# The first condition is redundant currently due to sorted comparison, but safe
                        if Cl2.l_bound < Cl1.lmin:
                                Cl1.lmin = Cl2.l_bound
                        elif Cl2.l_bound > Cl1.lmax:
                                Cl1.lmax = Cl2.l_bound
 	
            		if not indivOrient and Cl2.l_bound > Cl1.l_bound and Cl2.l_bound!=-1:
                		Cl1.l_bound = Cl2.l_bound
            		elif indivOrient and Cl2.l_bound < Cl1.l_bound and Cl2.l_bound!=-1:
                		Cl1.l_bound = Cl2.l_bound

        
	return Cl1
            
def calculateMargin(item):

		l_orient = int(item.C_type[0])
		r_orient = int(item.C_type[1])
		
		cl_margin_l = MEAN_D + DISC_ENHANCER2*DISC_D - 2*RDL - (item.lmax - item.lmin)
		cl_margin_r = MEAN_D + DISC_ENHANCER2*DISC_D - 2*RDL - (item.rmax - item.rmin)

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

    def __init__(self):
            
            self.l_bound=None
            self.r_bound=None
            self.l_bound_O=None  
            self.r_bound_O=None
            
            self.count=1
            self.C_type=None 
            self.Cmap_type = None 

            self.ltid = -1
            self.rtid = -1

	    # min and max locations by read
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

	    # locations of reads supporting cluster, and list of fragments supporting it, respectively	
	    self.readList = []	
	    self.support = []
			
    def __str__(self):
        return "%s %s %s %s %s %s %s %s %s" % (self.count, self.C_type, self.ltid, self.l_start, self.l_end, self.rtid, self.r_start, self.r_end, self.clsmall)

class Read(object):

	def __init__(self):
		name = -1
		l_read_bound = -1
		r_read_bound = -1
	def __str__(self):
		return "%s %s" %(self.l_read_bound, self.r_read_bound)

def writeSupport(fp,num,sList):

	writtenList = {}
	fp.write("%s" %num)
	for supportRead in sList:
		if supportRead not in writtenList:
			fp.write(" %s" %supportRead)
			writtenList[supportRead] = 1
	fp.write("\n")

def calcReadStats(CList):

	cluster_d_list = []

	for k,cluster in enumerate(CList):

		cluster.readList = sorted(cluster.readList,key = lambda member: member.l_read_bound)
		temp = cluster.readList
		prev = 0
		for k,item in enumerate(temp):
		   if k > 0:
			r_distance = item.l_read_bound - prev
			cluster_d_list.append(r_distance)
		   prev = item.l_read_bound
		
	cluster_d_list = sorted(cluster_d_list)
	f=open("../results/text/read_dist.txt","w")
	for item in cluster_d_list:
		f.write("%s\n" %item)
	f.close()

	bound1 = int(CLD_THRESH*len(cluster_d_list))
	bound2 = int(.65*len(cluster_d_list))
	if bound1 > 0:
		bound1 = bound1 - 1
	if len(cluster_d_list) == 0:
		print "WARNING: No clusters formed. Check input BAM file to see if it contains discordant reads."
		exit()

	print bound1, len(cluster_d_list)
	READ_BOUND1 = cluster_d_list[bound1]
	READ_BOUND2 = cluster_d_list[bound2]
	print "Read bounds are", READ_BOUND1, READ_BOUND2	
	return READ_BOUND1, READ_BOUND2

def writeClusters(CList, fp, fp2, READ_BOUND1):

	shift = 0
	prev_l_read = 0
	n_breaks = 0

	for k,cluster in enumerate(CList):

	     if CL_BREAK:
	   	cluster2 = -1
		broken = 0
		cluster_copy = copy.deepcopy(cluster)
 		for g,item in enumerate(cluster.readList):
			if item.l_read_bound - prev_l_read > READ_BOUND1 and g > 0:
				#print "Breaking cluster ", cluster, cluster.l_bound, cluster.r_bound
				#for elem in cluster.readList:
					#print elem
				# Break cluster by defining new one, and change support and bounds for both
				broken = 1
				cluster2 = Cluster()
				cluster2.ltid = cluster.ltid
				cluster2.rtid = cluster.rtid
				cluster2.C_type = cluster.C_type

				cluster.lmax = prev_l_read
				cluster.lmin = cluster.readList[0].l_read_bound
				cluster2.lmin = item.l_read_bound
				cluster2.lmax = cluster.readList[-1].l_read_bound
				cluster1_sorted = sorted(cluster.readList[:g],key = lambda member: member.r_read_bound)
				cluster2_sorted = sorted(cluster.readList[g:],key = lambda member: member.r_read_bound)
				cluster.rmax = cluster1_sorted[-1].r_read_bound
				cluster.rmin = cluster1_sorted[0].r_read_bound
				cluster2.rmax = cluster2_sorted[-1].r_read_bound
				cluster2.rmin = cluster2_sorted[0].r_read_bound

				if int(cluster.C_type[0]) == 0:
					cluster.l_bound = cluster.lmax
					cluster2.l_bound = cluster2.lmax
				elif int(cluster.C_type[0]) == 1:
					cluster.l_bound = cluster.lmin
					cluster2.l_bound = cluster2.lmin
				
				if int(cluster.C_type[1])==0:
					cluster.r_bound = cluster.rmax
					cluster2.r_bound = cluster2.rmax
				elif int(cluster.C_type[1])==1:
					cluster.r_bound = cluster.rmin
					cluster2.r_bound = cluster2.rmin
	
				cluster2.count = len(cluster.readList[g:])
				cluster.count = cluster.count - cluster2.count
				#print "into ", cluster, cluster.l_bound, cluster.r_bound, "and", cluster2, cluster2.l_bound, cluster2.r_bound
				cluster2.support = cluster.support[g:]
				cluster.support = cluster.support[:g]
				break
			prev_l_read = item.l_read_bound
		
		if not broken:
		   cluster.readList = sorted(cluster.readList,key = lambda member: member.r_read_bound)
		   for g,item in enumerate(cluster.readList):
                        if g > 0 and item.r_read_bound - prev_r_read > READ_BOUND1:
                                #print "Breaking cluster ", cluster, cluster.l_bound, cluster.r_bound
                                #for elem in cluster.readList:
                                        #print elem
                                # Break cluster by defining new one, and change support and bounds for both
                                cluster2 = Cluster()
                                cluster2.ltid = cluster.ltid
                                cluster2.rtid = cluster.rtid
                                cluster2.C_type = cluster.C_type

				cluster.rmax = prev_r_read
                                cluster.rmin = cluster.readList[0].r_read_bound
                                cluster2.rmin = item.r_read_bound
                                cluster2.rmax = cluster.readList[-1].r_read_bound
                                cluster1_sorted = sorted(cluster.readList[:g],key = lambda member: member.l_read_bound)
                                cluster2_sorted = sorted(cluster.readList[g:],key = lambda member: member.l_read_bound)
                                cluster.lmax = cluster1_sorted[-1].l_read_bound
                                cluster.lmin = cluster1_sorted[0].l_read_bound
                                cluster2.lmax = cluster2_sorted[-1].l_read_bound
                                cluster2.lmin = cluster2_sorted[0].l_read_bound

                                if int(cluster.C_type[0]) == 0:
                                        cluster.l_bound = cluster.lmax
                                        cluster2.l_bound = cluster2.lmax
                                elif int(cluster.C_type[0]) == 1:
                                        cluster.l_bound = cluster.lmin
					cluster2.l_bound = cluster2.lmin

                                if int(cluster.C_type[1])==0:
                                        cluster.r_bound = cluster.rmax
                                        cluster2.r_bound = cluster2.rmax
                                elif int(cluster.C_type[1])==1:
                                        cluster.r_bound = cluster.rmin
                                        cluster2.r_bound = cluster2.rmin

                                cluster2.count = len(cluster.readList[g:])
                                cluster.count = cluster.count - cluster2.count
                                #print "into ", cluster, cluster.l_bound, cluster.r_bound, "and", cluster2, cluster2.l_bound, cluster2.r_bound
                                cluster2.support = cluster.support[g:]
                                cluster.support = cluster.support[:g]
                                break
                        prev_r_read = item.r_read_bound
		# if everything seems to have gone OK
		if cluster2 != -1 and cluster.l_bound < cluster.r_bound and cluster2.l_bound < cluster2.r_bound:

			n_breaks+=1
			print "Writing broken clusters.", cluster.count, cluster2.count
			if cluster.count >= MIN_CLUSTER_SIZE:
				calculateMargin(cluster)
				fp.write("%s %s\n" %(1+k+shift,cluster))
				writeSupport(fp2,1+k+shift,cluster.support)
			if cluster2.count >= MIN_CLUSTER_SIZE:
				calculateMargin(cluster2)
				fp.write("%s %s\n" %(2+k+shift, cluster2))
				writeSupport(fp2,2+k+shift,cluster2.support)
				shift+=1
		elif cluster_copy.count >= MIN_CLUSTER_SIZE:
			
			calculateMargin(cluster_copy)
			fp.write("%s %s\n" %(1+k+shift,cluster_copy))
			writeSupport(fp2,1+k+shift,cluster_copy.support)

	     elif cluster.count >= MIN_CLUSTER_SIZE:
		calculateMargin(cluster)
		fp.write("%s %s\n" %(1+k+shift,cluster))
		writeSupport(fp2,1+k+shift,cluster.support)

	print n_breaks, "clusters were broken down and written.", CL_BREAK

if __name__== "__main__":
    
    f1=open("../results/text/All_Discords_P_S.txt","r")
    #f2=open("../results/text/All_Discords_I_S.txt","r")
    f3=open("../results/text/All_Clusters.txt","w")
    f4=open("../results/text/VariantMap.txt","w")
    f5=open("../results/text/Time_FC_loop.txt","w")
    f6=open("../results/text/Time_FC_line.txt","w")
    
    print "Function start"
    
    clusterList = []
    FList = []
    counter = 0
    newBlockS = 0
    offset = 0
    start_ref = 0

    totalent = readDistHash()
    fragHash = {}
    FG = nx.Graph()
    print "Cluster D:", CLUSTER_D
	
    for line_num,line in enumerate(f1):

        #start_fl = time.clock()
	
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
		
        print temp
        for x in range(ld):
	    print "Member:", FList[x]
 	    if FList[x].ltid == temp.ltid and FList[x].rtid == temp.rtid and FList[x].C_type == temp.C_type:

		    #compare 2 frags
		    f1lpos = FList[x].l_bound
		    f1rpos = FList[x].r_bound	    
		    f2lpos = temp.l_bound
		    f2rpos = temp.r_bound
		    #print x, f1lpos, f1rpos, f2lpos, f2rpos
		    EW = calcEdgeWeight(f1lpos, f1rpos, f2lpos, f2rpos, totalent)
		    #print "Weight is:", EW
		    Weights.append(EW)


		    if EW > EW_THRESH:
			print "Adding connection"
			FG.add_edge(FList[x].fragNum, temp.fragNum)
		
	# refresh list	
        if len(FList) > 1:
		
	 	if temp.ltid != FList[-2].ltid:
			
			#print "New TID: Processing ltid", temp.ltid, line, "Size of list:", len(DList)
			newBlockS = 0
                        FList = [FList[-1]]

		elif temp.ltid == FList[newBlockS].ltid and (temp.l_bound - FList[newBlockS].l_bound) > 2*CLUSTER_D:

                    FList = FList[newBlockS:]
                    newBlockS = len(FList)-1
	
            
    fc = open("../results/text/cluster_cliques.txt","w")
    #print "Graph is:", FG.nodes(),"\n", FG.edges()
    #nx.draw_networkx(FG)
    max_clique_lists = list(nx.find_cliques(FG))
    #max_clique_lists = list(nx.k_clique_communities(FG, MIN_CLUSTER_SIZE))
    #print "Max cliques:", max_clique_lists
    for k,elem in enumerate(max_clique_lists):
	elem0 = list(elem)[0]
	fc.write("%s\n" %("@Cluster"+str(k)))
	ltid = fragHash[elem0].ltid
	rtid = fragHash[elem0].rtid
	lbp = fragHash[elem0].l_bound
	rbp = fragHash[elem0].r_bound

	for item in list(elem):
		fc.write("%s\n" %str(fragHash[item]))
		if fragHash[item].l_bound > lbp:
			lbp = fragHash[item].l_bound
		if fragHash[item].r_bound < rbp:
                        rbp = fragHash[item].r_bound
		#Compute and write sum of total edge weights for this node to all other nodes in cliquw
		#fc.write
	fc.write("%s\t%s\t%s\t%s\n" %(ltid,lbp, rtid, rbp))

    #f2.close()

    #[RB1, RB2] = calcReadStats(clusterList)
    #writeClusters(clusterList, f3, f4, RB1)
 
    f3.close()
    f4.close()            
   
       
