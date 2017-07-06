# .6 s on small data set w/ refRate of 5

from collections import Counter
import sys
import time

MARGIN = int(sys.argv[1]) #30
SLOP = int(sys.argv[2]) #80 max(maxBPwidth,100)
refRate = int(sys.argv[3]) # 5
glVarNum = 1

f=open("../results/text/bam_stats.txt","r")
temp = f.readline()
MEAN_IL = int(f.readline().split()[0])

class Vertex:
    def __init__(self, node):
        self.id = node
        self.adjacent = {}

    def __str__(self):
        return str(self.id) + ' adjacent: ' + str([x.id for x in self.adjacent])

    def add_neighbor(self, neighbor, weight=0):
        self.adjacent[neighbor] = weight

    def get_connections(self):
        return self.adjacent

    def get_id(self):
        return self.id

    def get_weight(self, neighbor):
        return self.adjacent[neighbor]

class Graph:
    def __init__(self):
        self.vert_dict = {}
        self.num_vertices = 0

    def __iter__(self):
        return iter(self.vert_dict.values())

    def add_vertex(self, node):
        self.num_vertices = self.num_vertices + 1
        new_vertex = Vertex(node)
        self.vert_dict[node] = new_vertex
        return new_vertex

    def get_vertex(self, n):
        if n in self.vert_dict:
            return self.vert_dict[n]
        else:
            return None

    def add_edge(self, frm, to, cost = 0):
        if frm not in self.vert_dict:
            self.add_vertex(frm)
        if to not in self.vert_dict:
            self.add_vertex(to)

        self.vert_dict[frm].add_neighbor(to, cost)
        self.vert_dict[to].add_neighbor(frm, cost)

    def get_vertices(self):
        return self.vert_dict

class Cluster(object):

    def __init__(self, data):
        data_split = data.split()
        self.l_start = int(data_split[4])
        self.l_end = int(data_split[5])
        self.r_start = int(data_split[8])
        self.r_end = int(data_split[9])
        self.l_orient = int(data_split[2][0])
        self.r_orient = int(data_split[2][1])

        self.typeC = data_split[2]
        self.ltid = data_split[3]
        self.rtid = data_split[7]
        self.mapNum = int(data_split[0])

    def __str__(self):
        return "%s %s %s %s %s %s %s %s %s %s" % (self.mapNum, self.l_start,self.l_end, self.r_start, self.r_end, self.typeC, self.ltid, self.rtid, self.l_orient, self.r_orient)

class OverlappingCluster(object):

    # at most 4 clusters will overlap
    # $ move to init
    bp1_start = -1
    bp2_start = -1
    bp3_start = -1
    bp1_end = -1
    bp2_end = -1
    bp3_end = -1
    count = 2
    typeO = None

    # it is an insertion cluster if overlapping point has had both orientation of reads
    bp1f = 0
    bp1r = 0
    bp1tid = -1
    bp2tid = -1
    bp3tid = -1

    def __init__(self):
        self.clusterNums = [] # original cluster numbers that support this variant
        self.support = "PE"
        self.complete = 0 # All clusters that can match with this variant have arrived and matched
        self.varNum = None

    def __str__(self):
        #return "%s %s %s %s %s %s %s %s %s" % (self.typeO, self.bp1, self.bp2,self.bp3, self.bp1tid, self.bp2tid, self.bp3tid, self.count, self.support)
        return "%s %s %s %s %s %s %s %s %s %s %s %s" % (self.typeO, self.bp1tid, self.bp1_start, self.bp1_end, self.bp2tid , self.bp2_start, self.bp2_end, self.bp3tid, self.bp3_start, self.bp3_end, self.support, self.count)

def isOverlapping(CorV, cl1, cl2, LR):

    # Add small SLOP to one bp to account for sequencing errors
    l1_start = cl1.l_start - 0
    l1_end = cl1.l_end + 0
    r1_start = cl1.r_start - 0
    r1_end = cl1.r_end + 0
    
    if CorV == "C":

        #$ Add condition to check if contained fully within
        if LR=="LL":

            if (l1_start >= cl2.l_start and l1_start <= cl2.l_end) or (l1_end >= cl2.l_start and l1_end <= cl2.l_end) or abs(cl1.l_start - cl2.l_end) < SLOP or abs(cl2.l_start - cl1.l_end) < SLOP:

                return 1
            
        if LR=="LR":

            # $ Can write a function called intervalOverlap and pass these 4 points in
            if (l1_start >= cl2.r_start and l1_start <= cl2.r_end) or (l1_end >= cl2.r_start and l1_end <= cl2.r_end) or abs(cl1.l_start - cl2.r_end) < SLOP or abs(cl2.r_start - cl1.l_end) < SLOP:

                return 1
            
        if LR=="RL":

            if (r1_start >= cl2.l_start and r1_start <= cl2.l_end) or (r1_end >= cl2.l_start and r1_end <= cl2.l_end) or abs(cl1.r_start - cl2.l_end) < SLOP or abs(cl2.l_start - cl1.r_end) < SLOP:

                return 1

        if LR=="RR":

            if (r1_start >= cl2.r_start and r1_start <= cl2.r_end) or (r1_end >= cl2.r_start and r1_end <= cl2.r_end) or abs(cl1.r_start - cl2.r_end) < SLOP or abs(cl2.r_start - cl1.r_end) < SLOP:

                return 1

        return 0

    if CorV == "V":

        if LR=="L1":

            if (l1_start >= cl2.bp1_start and l1_start <= cl2.bp1_end) or (l1_end >= cl2.bp1_start and l1_end <= cl2.bp1_end) or abs(cl1.l_start - cl2.bp1_end) < SLOP or abs(cl2.bp1_start - cl1.l_end) < SLOP:

                return 1
            
        if LR=="L2":

            if (l1_start >= cl2.bp2_start and l1_start <= cl2.bp2_end) or (l1_end >= cl2.bp2_start and l1_end <= cl2.bp2_end) or abs(cl1.l_start - cl2.bp2_end) < SLOP or abs(cl2.bp2_start - cl1.l_end) < SLOP:

                return 1
            
        if LR=="L3":

            if (l1_start >= cl2.bp3_start and l1_start <= cl2.bp3_end) or (l1_end >= cl2.bp3_start and l1_end <= cl2.bp3_end) or abs(cl1.l_start - cl2.bp3_end) < SLOP or abs(cl2.bp3_start - cl1.l_end) < SLOP:

                return 1

        if LR=="R1":

            if (r1_start >= cl2.bp1_start and r1_start <= cl2.bp1_end) or (r1_end >= cl2.bp1_start and r1_end <= cl2.bp1_end) or abs(cl1.r_start - cl2.bp1_end) < SLOP or abs(cl2.bp1_start - cl1.r_end) < SLOP:

                return 1

        if LR=="R2":

            if (r1_start >= cl2.bp2_start and r1_start <= cl2.bp2_end) or (r1_end >= cl2.bp2_start and r1_end <= cl2.bp2_end) or abs(cl1.r_start - cl2.bp2_end) < SLOP or abs(cl2.bp2_start - cl1.r_end) < SLOP:

                return 1

        if LR=="R3":

            if (r1_start >= cl2.bp3_start and r1_start <= cl2.bp3_end) or (r1_end >= cl2.bp3_start and r1_end <= cl2.bp3_end) or abs(cl1.r_start - cl2.bp3_end) < SLOP or abs(cl2.bp3_start - cl1.r_end) < SLOP:

                return 1
            
        return 0

def setBPs(cl1, cl2, LR):


    if LR == "LL":

        sorted_bps = sorted([cl1.l_start, cl1.l_end, cl2.l_start, cl2.l_end])
        return sorted_bps[1], sorted_bps[2]

    elif LR == "LR":

        sorted_bps = sorted([cl1.l_start, cl1.l_end, cl2.r_start, cl2.r_end])
        return sorted_bps[1], sorted_bps[2]

    elif LR == "RL":

        sorted_bps = sorted([cl1.r_start, cl1.r_end, cl2.l_start, cl2.l_end])
        return sorted_bps[1], sorted_bps[2]

    elif LR == "RR":

        sorted_bps = sorted([cl1.r_start, cl1.r_end, cl2.r_start, cl2.r_end])
        return sorted_bps[1], sorted_bps[2]

    return -1,-1

def swapBPs(bpList, tidList):

    temp_s = bpList[0]
    temp_e = bpList[1]
    bpList[0] = bpList[2]
    bpList[1] = bpList[3]
    bpList[2] = temp_s
    bpList[3] = temp_e
    temp = tidList[0]
    tidList[0] = tidList[1]
    tidList[1] = temp
        
            
def populateCluster(fp, listName, LR):

    first = 1

    for line in fp:
        
        if line[0:2]!= "-1":
            
            temp = Cluster(line)

            if LR == "left":
                tid = temp.ltid
            elif LR == "right":
                tid = temp.rtid
            else:
                print "incorrect parameter to populateCluster LR"
                return
            
            if first:
                curentTID = tid
                first = 0
            if tid == currentTID:
                listName.append(temp)
            else:
                fp.seek(last_pos)
                break
            
            last_pos = fp.tell()

def writeVariants(varList, fpAV, fpCVM, hashedVM, offset):

    #print "In WriteVariants. VarList length:", len(varList)
    for g, item in enumerate(varList):

        num = offset + 1 + g
        fpAV.write("%s %s\n" %(num,item))
        fpCVM.write("%s" %(num))

        if item.typeO == "INS" or item.typeO == "INS_I" or item.typeO == "INS_C_P" or item.typeO == "INS_C_I_P":

                # The first condition should be ensured anyway
                if item.bp2tid == item.bp3tid and item.bp3_start > item.bp2_start:

                    bpList = [item.bp2_start, item.bp2_end, item.bp3_start, item.bp3_end]
                    swapBPs(bpList, [0, 0])
                    [item.bp2_start, item.bp2_end, item.bp3_start, item.bp3_end] = bpList

        elif item.typeO == "TD":

            # The first condition should be ensured anyway
            if item.bp1tid == item.bp2tid and item.bp1_start > item.bp2_start:

                bpList = [item.bp1_start, item.bp1_end, item.bp2_start, item.bp2_end]
                swapBPs(bpList, [0,0])
                [item.bp1_start, item.bp1_end, item.bp2_start, item.bp2_end] = bpList
            
                
        for elem in item.clusterNums:
                
         try:
             
            for elem2 in hashedVM[elem]:
                    
                #print "Writing", elem2, "from cluster", elem
                fpCVM.write(" %s" %elem2)
         except:
             print "Exception writing in Variant Map from cluster", elem
                
        fpCVM.write("\n")
            

def file_len(f):
    
    for i, l in enumerate(f):
        pass
    return i + 1

def compareCluster(cluster1, Clusters, Claimed, graph_c, graph_m, offset, LR, OCArray):

    debug = 0
    global glVarNum
    AnyMatch = 0 
    cloc1_s = cluster1.l_start
    cloc1_e = cluster1.l_end
    cloc2_s = cluster1.r_start
    cloc2_e = cluster1.r_end
    cloc1_tid = cluster1.ltid
    cloc2_tid = cluster1.rtid
 
    lorient1 = int(cluster1.l_orient)
    rorient1 = int(cluster1.r_orient)
    cl2 = 'c' + str(cluster1.mapNum)
 
    for x, item in enumerate(Clusters):

        if cluster1.mapNum == Clusters[x].mapNum:
            continue

        cl1 = 'c' + str(Clusters[x].mapNum)

        # If 2 clusters have been compared do not compare them again. See exception 1 at end of file. 
        if graph_c.get_vertex(cl2) != None:
            if cl1 in graph_c.get_vertex(cl2).get_connections():
                continue

        # If 2 clusters share variant, then comparing them is redundant. Exception 2 noted at end of file:
        # Any subtleties can be addressed when comparing to existing variant, e.g. to branch into "MQ-0" variant.
        skip = 0
        if graph_m.get_vertex(cl2) != None and graph_m.get_vertex(cl1) != None:
            for variant in graph_m.get_vertex(cl2).get_connections():
                if variant in graph_m.get_vertex(cl1).get_connections():
                    skip = 1
                    break

        if skip:
            continue
        
        newCl = 0 # this is flag for new variant

        # All comparisons

        cloc3_s = Clusters[x].l_start
        cloc3_e = Clusters[x].l_end
        cloc4_s = Clusters[x].r_start
        cloc4_e = Clusters[x].r_end
        cloc3_tid = Clusters[x].ltid
        cloc4_tid = Clusters[x].rtid
        
        lorient2 = Clusters[x].l_orient
        rorient2 = Clusters[x].r_orient

        sign1=0
        sign2=0
        sign3=0
        sign4=0
        bp_1_start = -1
        bp_1_end = -1
        bp_2_start = -1
        bp_2_end = -1
        bp_3_start = -1
        bp_3_end = -1
        bp1_f = -1
        bp1_r = -1
        bp_1_tid = -1
        bp_2_tid = -1
        bp_3_tid = -1

        newCl = 0
        type_O = "none"

        if cluster1.ltid == Clusters[x].ltid and isOverlapping("C", cluster1, Clusters[x], "LL"): 

            sign1 = 1

        if cluster1.ltid == Clusters[x].rtid and isOverlapping("C", cluster1, Clusters[x], "LR"): 

            sign2 = 1
            
        if cluster1.rtid == Clusters[x].ltid and isOverlapping("C", cluster1, Clusters[x], "RL"): 

            sign3 = 1

        if cluster1.rtid == Clusters[x].rtid and isOverlapping("C", cluster1, Clusters[x], "RR"): 

            sign4 = 1


        if (LR == "LL" and not sign1) or (LR == "LR" and not sign2) or (LR =="RL" and not sign3) or (LR == "RR" and not sign4):
            continue
            
        if not sign1 and not sign2 and not sign3 and not sign4:
            continue

        graph_c.add_edge(cl1, cl2, 1)
	
	Cl1_small = 0
	Cl2_small = 0
        
	# flanking cluster for all small SVs 
	if cluster1.l_orient == 0 and cluster1.r_orient == 1 and cluster1.r_end - cluster1.l_start < MEAN_IL:
		Cl1_small = 1
	if Clusters[x].l_orient == 0 and Clusters[x].r_orient == 1 and Clusters[x].r_end - Clusters[x].l_start < MEAN_IL:
		Cl2_small = 1

        #$ Add checks to see if TID same for INV, TD and if at least 2 of 3 same for INS		

        if sign1 and sign4 and cluster1.l_orient != Clusters[x].l_orient and cluster1.r_orient != Clusters[x].r_orient and cluster1.r_orient == cluster1.l_orient and cluster1.ltid == cluster1.rtid:
                
            type_O = "INV"
            newCl=1

            # bp1 is the left alignments
            [bp_1_start, bp_1_end] = setBPs(cluster1, Clusters[x], "LL")
            [bp_2_start, bp_2_end] = setBPs(cluster1, Clusters[x], "RR")
            [bp_1_tid, bp_2_tid]  = cloc1_tid, cloc2_tid

        #$ can add condition to check for opposite orientation of clusters[i] mates even though v unlikely to not be case
        # next 4 ifs are for small insertions;  
        elif sign2 and sign4 and not sign1 and not sign3 and Cl1_small and not Cl2_small:

            #$ strictly, this should be TD_I_I but for now clubbing all inverted duplications/insertions into INS_I

            if Clusters[x].l_orient == Clusters[x].r_orient:
                type_O = "INS_I"
            elif not (cluster1.ltid == cluster1.rtid == Clusters[x].ltid == Clusters[x].rtid):
                type_O = "INS"              
	    elif cluster1.l_orient == 0 and cluster1.r_orient == 1 and Clusters[x].l_orient == 0 and Clusters[x].r_orient == 1:
		type_O = "INS" 
            else:    
                type_O = "TD_I"
                # TD or INS
                
            newCl=1
            #bp_1 numbering is crucial only
            # $ can make more precise by evaluating which of 3 bp pairs yields least bp width
            [bp_1_start, bp_1_end] = setBPs(cluster1, Clusters[x], "LR")
            [bp_2_start, bp_2_end] = cloc3_s, cloc3_e
            bp1_f = not rorient2
            bp1_r = rorient2
            [bp_1_tid, bp_2_tid] = cloc1_tid, cloc3_tid
            
            
        elif sign3 and sign4 and not sign1 and not sign2 and Cl2_small and not Cl1_small:

            if cluster1.l_orient == cluster1.r_orient:
                type_O = "INS_I"
            elif not (cluster1.ltid == cluster1.rtid == Clusters[x].ltid == Clusters[x].rtid):
                type_O = "INS"
	    elif cluster1.l_orient == 0 and cluster1.r_orient == 1 and Clusters[x].l_orient == 0 and Clusters[x].r_orient == 1:
                type_O = "INS"
            else:    
                type_O = "TD_I"
                
            newCl=1

            [bp_1_start, bp_1_end] = setBPs(cluster1, Clusters[x], "RL")
            [bp_2_start, bp_2_end] = cloc1_s, cloc1_e
            bp1_f = not rorient1
            bp1_r = rorient1
            bp_1_tid = cloc2_tid
            bp_2_tid = cloc1_tid
            
        elif sign1 and sign3 and not sign2 and not sign4 and Cl1_small and not Cl2_small:

            if Clusters[x].l_orient == Clusters[x].r_orient:
                type_O = "INS_I"
            elif not (cluster1.ltid == cluster1.rtid == Clusters[x].ltid == Clusters[x].rtid):
                type_O = "INS"
	    elif cluster1.l_orient == 0 and cluster1.r_orient == 1 and Clusters[x].l_orient == 0 and Clusters[x].r_orient == 1:
                type_O = "INS"
            else:    
                type_O = "TD_I"
            
            newCl=1
            [bp_1_start, bp_1_end] = setBPs(cluster1, Clusters[x], "RL")
            [bp_2_start, bp_2_end] = cloc4_s, cloc4_e
            bp1_f = not lorient2
            bp1_r = lorient2
            bp_1_tid = cloc1_tid
            bp_2_tid = cloc4_tid
    
        elif sign1 and sign2 and not sign3 and not sign4 and Cl2_small and not Cl1_small:

            if cluster1.l_orient == cluster1.r_orient:
                type_O = "INS_I"
            if not (cluster1.ltid == cluster1.rtid == Clusters[x].ltid == Clusters[x].rtid):
                type_O = "INS"
	    elif cluster1.l_orient == 0 and cluster1.r_orient == 1 and Clusters[x].l_orient == 0 and Clusters[x].r_orient == 1:
                type_O = "INS"
            else:    
                type_O = "TD_I"
                
            newCl=1
            [bp_1_start, bp_1_end] = setBPs(cluster1, Clusters[x], "LR")
            [bp_2_start, bp_2_end] = cloc2_s, cloc2_e
            bp1_f = not lorient1
            bp1_r = lorient1
            bp_1_tid = cloc1_tid
            bp_2_tid = cloc2_tid
            
        # Large insertions
        elif sign1 and cluster1.l_orient != Clusters[x].l_orient and (cluster1.ltid == cluster1.rtid or cluster1.ltid == Clusters[x].rtid or cluster1.rtid == Clusters[x].rtid):
                    
            #$ Can add a check here to see if the breakpoints are not all on different chromosomes. Two must be on same chromosome. 3rd is bp1 if on diff chr, and this is, so far, copy-paste.
            if cluster1.l_orient == cluster1.r_orient and Clusters[x].l_orient == Clusters[x].r_orient:
                type_O = "INS_I"
            elif cluster1.l_orient != cluster1.r_orient and Clusters[x].l_orient != Clusters[x].r_orient:
                type_O = "INS"


            if type_O == "INS" or type_O == "INS_I":
                
                newCl=1
                [bp_1_start, bp_1_end] = setBPs(cluster1, Clusters[x], "LL")
                [bp_2_start, bp_2_end] = cloc2_s, cloc2_e
                [bp_3_start, bp_3_end] = cloc4_s, cloc4_e

                bp_1_tid = cloc1_tid
                bp_2_tid = cloc2_tid
                bp_3_tid = cloc4_tid

                # Call breakpoint that is one of the two cut bp's bp2. Later, bp1 and bp3 will be switched if necessary in RD and SR verification.
                # if middle breakpoint is overlapping, then has to be a cut insertion
                if bp_1_tid == bp_2_tid and bp_2_tid == bp_3_tid:

                    if bp_1_start > min(bp_2_end, bp_3_end) and bp_1_end < max(bp_2_start, bp_3_start):

                        if type_O == "INS":
                            type_O = "INS_C"
                        elif type_O == "INS_I":
                            type_O = "INS_C_I"

                        temp_s = bp_2_start
                        temp_e = bp_2_end
                        bp_2_start = bp_1_start
                        bp_2_end = bp_1_end

                        # this is conventional
                        bp_1_start = min(temp_s, bp_3_start)
                        bp_1_end = min(temp_e, bp_3_end)
                        bp_3_start = max(temp_s, bp_3_start)
                        bp_3_end = max(temp_e, bp_3_end)
                        
                        if cloc2_s < cloc4_s:
                            bp_1_tid = cloc2_tid
                            bp_3_tid = cloc3_tid
                        else:
                            bp_1_tid = cloc3_tid
                            bp_3_tid = cloc2_tid
                            
                        bp_2_tid = cloc1_tid
                
                # if overlapping breakpoint occurs on the same chromosome as the second bp but not the 3rd, this must be a cut
                elif bp_1_tid == bp_2_tid and bp_2_tid != bp_3_tid:

                    # breakpoint locations are confirmed
                    if type_O == "INS":
                        type_O = "INS_C_P"
                    elif type_O == "INS_I":
                        type_O = "INS_C_I_P"

                    # swap bp1 and bp3
                    bpList = [bp_1_start, bp_1_end, bp_3_start, bp_3_end]
                    tidList = [bp_1_tid, bp_3_tid]
                    swapBPs(bpList, tidList)
                    [bp_1_start, bp_1_end, bp_3_start, bp_3_end] = bpList
                    [bp_1_tid, bp_3_tid] = tidList


                elif bp_1_tid == bp_3_tid and bp_3_tid != bp_2_tid:
                    
                    if type_O == "INS":
                        type_O = "INS_C_P"
                    elif type_O == "INS_I":
                        type_O = "INS_C_I_P"

                    # swap bp1 and bp2
                    bpList = [bp_1_start, bp_1_end, bp_2_start, bp_2_end]
                    tidList = [bp_1_tid, bp_2_tid]
                    swapBPs(bpList, tidList)
                    [bp_1_start, bp_1_end, bp_2_start, bp_2_end] = bpList
                    [bp_1_tid, bp_2_tid] = tidList

                    

        elif sign4 and cluster1.r_orient != Clusters[x].r_orient and (cluster1.rtid == cluster1.ltid or cluster1.rtid == Clusters[x].ltid or cluster1.ltid == Clusters[x].ltid):

            if cluster1.l_orient == cluster1.r_orient and Clusters[x].l_orient == Clusters[x].r_orient:
                type_O = "INS_I"
            elif cluster1.l_orient != cluster1.r_orient and Clusters[x].l_orient != Clusters[x].r_orient:
                type_O = "INS"

            if debug:
                print type_O

            if type_O == "INS" or type_O == "INS_I":
            
                newCl=1
                [bp_1_start, bp_1_end] = setBPs(cluster1, Clusters[x], "RR")
                [bp_2_start, bp_2_end] = cloc1_s, cloc1_e
                [bp_3_start, bp_3_end] = cloc3_s, cloc3_e
                bp_1_tid = cloc2_tid
                bp_2_tid = cloc1_tid
                bp_3_tid = cloc3_tid
                

                if bp_1_tid == bp_2_tid and bp_2_tid == bp_3_tid:
                    
                    if bp_1_start > min(bp_2_end, bp_3_end) and bp_1_start < max(bp_2_start, bp_3_start):

                        if type_O == "INS":
                            type_O = "INS_C"
                        elif type_O == "INS_I":
                            type_O = "INS_C_I"

                        temp_s = bp_2_start
                        temp_e = bp_2_end
                        bp_2_start = bp_1_start
                        bp_2_end = bp_1_end
                        
                        #conventional, used later
                        bp_1_start = min(temp_s, bp_3_start)
                        bp_1_end = min(temp_e, bp_3_end)
                        bp_3_start = max(temp_s, bp_3_start)
                        bp_3_end = max(temp_e, bp_3_end)
                
                        if cloc1_s < cloc3_s:
                            bp_1_tid = cloc1_tid
                            bp_3_tid = cloc3_tid
                        else:
                            bp_1_tid = cloc3_tid
                            bp_3_tid = cloc1_tid
                            
                        bp_2_tid = cloc2_tid

                # if overlapping breakpoint occurs on the same chromosome as exactly one other breakpoint, this must be a cut
                elif bp_1_tid == bp_2_tid and bp_2_tid != bp_3_tid:

                    # breakpoint locations are confirmed
                    if type_O == "INS":
                        type_O = "INS_C_P"
                    elif type_O == "INS_I":
                        type_O = "INS_C_I_P"
                    
                    # swap bp1 and bp3
                    bpList = [bp_1_start, bp_1_end, bp_3_start, bp_3_end]
                    tidList = [bp_1_tid, bp_3_tid]
                    swapBPs(bpList, tidList)
                    [bp_1_start, bp_1_end, bp_3_start, bp_3_end] = bpList
                    [bp_1_tid, bp_3_tid] = tidList


                elif bp_1_tid == bp_3_tid and bp_3_tid != bp_2_tid:
                    
                    if type_O == "INS":
                        type_O = "INS_C_P"
                    elif type_O == "INS_I":
                        type_O = "INS_C_I_P"

                    # swap bp1 and bp2
                    bpList = [bp_1_start, bp_1_end, bp_2_start, bp_2_end]
                    tidList = [bp_1_tid, bp_2_tid]
                    swapBPs(bpList, tidList)
                    [bp_1_start, bp_1_end, bp_2_start, bp_2_end] = bpList
                    [bp_1_tid, bp_2_tid] = tidList
                    
##
##        elif sign2 and not sign3 and not sign4 and not sign1 and cluster1.l_orient != Clusters[x].r_orient and cluster1.l_orient != cluster1.r_orient:
##
##                    type_O = "INS_C"
##                    if Clusters[x].l_orient == Clusters[x].r_orient:
##                        type_O = "INS_C_I"
##                    newCl=1
##                    bp_2 = min(cloc1, cloc4)
##                    bp_3 = cloc2
##                    bp_2_tid = cloc2_tid
##                    bp_3_tid = cloc2_tid
##                    bp1 = cloc3
##                    bp_1_tid = cloc3_tid
##
##        elif sign3 and not sign2 and not sign4 and not sign1 and cluster1.r_orient != Clusters[x].l_orient and Clusters[x].l_orient != Clusters[x].r_orient:
##
##                    type_O = "INS_C"
##                    if cluster1.l_orient == cluster1.r_orient:
##                        type_O = "INS_C_I"
##                    newCl=1
##                    bp_2 = min(cloc2, cloc3)
##                    bp_3 = cloc4
##                    bp_2_tid = cloc2_tid
##                    bp_3_tid = cloc4_tid
##                    bp1 = cloc1
##                    bp_1_tid = cloc1_tid


        # What to do if matches

        if newCl:
            
            temp = OverlappingCluster()
            temp.typeO = type_O
            temp.bp1_start = bp_1_start
            temp.bp1_end = bp_1_end
            temp.bp2_start = bp_2_start
            temp.bp2_end = bp_2_end
            temp.bp1tid = bp_1_tid
            temp.bp2tid = bp_2_tid
            if bp_3_start != -1:
                temp.bp3_start = bp_3_start
                temp.bp3_end = bp_3_end
                temp.bp3tid = bp_3_tid
            if bp1_f!= -1:
                temp.bp1f = bp1_f
            if bp1_r!= -1:
                temp.bp1r = bp1_r
    
            temp.clusterNums.append(cluster1.mapNum)
            temp.clusterNums.append(Clusters[x].mapNum)
            temp.varNum = glVarNum
            glVarNum+=1
            OCArray.append(temp)
            
            Claimed[Clusters[x].mapNum] = 1

            AnyMatch = 1

            v1 = 'v' + str(temp.varNum)
            graph_m.add_edge(cl1,v1,1)
            graph_m.add_edge(cl2,v1,1)
            graph_c.add_edge(cl1,v1,1)
            graph_c.add_edge(cl2,v1,1)

            #print "Match:", temp, " from", cluster1, "and", Clusters[x]

            #print len(graph_m.get_vertex(v1).get_connections()), "after comparison"
                           
            # $ Can break here if do comparison for partial bp overlap (with 1 bp or more) within each variant and call it a new variant if partially overlaps properly.
            #break

    if AnyMatch:
        Claimed[cluster1.mapNum] = 1

def compareVariant(cluster1, varList, Claimed, graph_c, graph_m, offset, LR):

    AnyMatch = 0
    cloc1_s = cluster1.l_start
    cloc1_e = cluster1.l_end
    cloc2_s = cluster1.r_start
    cloc2_e = cluster1.r_end
    cloc1_tid = cluster1.ltid
    cloc2_tid = cluster1.rtid
 
    lorient1 = int(cluster1.l_orient)
    rorient1 = int(cluster1.r_orient)
    cl1 = 'c' + str(cluster1.mapNum)
    
    for g,elem in enumerate(varList):

        match = 0

        v1 = 'v' + str(elem.varNum)
        
        #print len(graph_m.get_vertex(v1).get_connections()), "before comparison"

        # If cluster compared to any cluster through variant, do not compare them.            
        # For later, when implement variant branching as noted at end of file. Reduced TPs by 2, same FPs, in small data set.
##        for connection in graph_m.get_vertex(v1).get_connections():
##            graph_c.add_edge(cl1,connection,1)

        # If cluster has been compared to variant, do not compare
        if graph_c.get_vertex(cl1) != None and graph_c.get_vertex(v1) != None:
            if cl1 in graph_c.get_vertex(v1).get_connections():
                continue
        
        # If cluster compared to all clusters currently in variant, no need to compare again. Implement variant branching later to resolve unlikely exceptions to this, if any.
        skip = 0

        if graph_c.get_vertex(cl1) !=None and graph_m.get_vertex(v1) != None:

            skip = 1
            
            for connection in graph_m.get_vertex(v1).get_connections():

                if connection not in graph_c.get_vertex(cl1).get_connections():

                    skip = 0
                    break

        if skip:
            continue

        # Don't compare if exceeds left-sorted comparison bounds. Will appear again for right bound comparison later.
        if LR == "L" and (elem.bp1tid != cluster1.ltid or abs(elem.bp1_start - cluster1.l_start) > MARGIN) and (elem.bp2tid != cluster1.ltid or abs(elem.bp2_start - cluster1.l_start) > MARGIN) and (elem.bp3tid != cluster1.ltid or abs(elem.bp3_start - cluster1.l_start) > MARGIN):
            continue
        
        # Analogous to above
        if LR == "R" and (elem.bp1tid != cluster1.rtid or abs(elem.bp1_start - cluster1.r_start) > MARGIN) and (elem.bp2tid != cluster1.rtid or abs(elem.bp2_start - cluster1.r_start) > MARGIN) and (elem.bp3tid != cluster1.rtid or abs(elem.bp3_start - cluster1.r_start) > MARGIN):
            continue

        graph_c.add_edge(cl1,v1,1)

       
        if elem.typeO == "TD_I":
            
            # small-medium insertions
            if isOverlapping("V", cluster1, elem, "L1")  and isOverlapping("V", cluster1, elem, "R1") and cloc1_tid == elem.bp1tid and cloc2_tid == elem.bp1tid:

                elem.count+= 1
                match = 1
                AnyMatch = 1

            elif isOverlapping("V", cluster1, elem, "L2") and isOverlapping("V", cluster1, elem, "R2") and cloc1_tid == elem.bp2tid and cloc2_tid == elem.bp2tid:

                elem.typeO = "TD"
                elem.count+=1
                match = 1
                AnyMatch = 1
                elem.complete = 1

                bpList = [elem.bp1_start, elem.bp1_end, elem.bp2_start, elem.bp2_end]
                tidList = [elem.bp1tid, elem.bp2tid]
                swapBPs(bpList, tidList)
                [elem.bp1_start, elem.bp1_end, elem.bp2_start, elem.bp2_end] = bpList
                [elem.bp1tid, elem.bp2tid] = tidList

            # all insertions
            elif isOverlapping("V", cluster1, elem, "L1") and cloc1_tid == elem.bp1tid:
                
                # this is a good check for diploid (non-)conflicts if insertion occurs in both chromosomes in same location
                if (not lorient1 and elem.bp1r) or (lorient1 and elem.bp1f):

                    elem.bp1f = 1
                    elem.bp1r = 1
                    elem.typeO = "INS"

                    if cluster1.l_orient == cluster1.r_orient:
                        elem.typeO ="INS_I"

                    #if elem.bp3 == -1:
                    [elem.bp3_start, elem.bp3_end] = cloc2_s, cloc2_e
                    elem.bp3tid = cloc2_tid
                    elem.count+=1
                    match = 1
                    AnyMatch = 1

                # may be redundant b/c one of them is already set    
                elif not lorient1:

                    elem.bp1f = 1

                else:

                    elem.bp1r = 1                   
                

            elif isOverlapping("V", cluster1, elem, "R1") and cloc2_tid == elem.bp1tid:

                if (not rorient1 and elem.bp1r) or (rorient1 and elem.bp1f):

                    elem.bp1f = 1
                    elem.bp1r = 1
                    elem.typeO = "INS"

                    if cluster1.l_orient == cluster1.r_orient:
                        elem.typeO ="INS_I"

                    #if elem.bp3 == -1:
                    [elem.bp3_start, elem.bp3_end] = cloc1_s, cloc1_e
                    elem.bp3tid = cloc1_tid
                    elem.count+=1
                    match = 1
                    AnyMatch = 1

                # may be redundant b/c one of them is already set    
                elif not lorient1:

                    elem.bp1f = 1

                else:

                    elem.bp1r = 1
                    
            if elem.bp3_start == -1:

                if isOverlapping("V", cluster1, elem, "L2") and cloc1_tid == elem.bp2tid and (not isOverlapping("V", cluster1, elem, "R2") or cloc2_tid != elem.bp2tid) and cluster1.l_orient != cluster1.r_orient:

                    match = 1
                    elem.count+=1
                    AnyMatch = 1
                    [elem.bp3_start, elem.bp3_end] = cloc2_s, cloc2_e
                    elem.bp3tid = cloc2_tid
                    if elem.typeO == "INS_I":
                        elem.typeO = "INS_C_I"
                    elif elem.typeO == "INS":
                        elem.typeO = "INS_C"

                # won't really occur but for future just in case

                elif not isOverlapping("V", cluster1, elem, "L2") and cloc1_tid == elem.bp2tid and not( isOverlapping("V", cluster1, elem, "R2") or cloc2_tid != elem.bp2tid) and cluster1.l_orient != cluster1.r_orient:

                    match = 1
                    elem.count+=1
                    AnyMatch = 1
                    [elem.bp3_start, elem.bp3_end] = cloc1_s, cloc1_e
                    elem.bp3tid = cloc1_tid
                    if elem.typeO == "INS_I":
                        elem.typeO = "INS_C_I"
                    elif elem.typeO == "INS":
                        elem.typeO = "INS_C"


        if elem.typeO == "INS_C" or elem.typeO == "INS_C_I":

            # small or medium insertions
            if isOverlapping("V", cluster1, elem, "L1") and cloc1_tid == elem.bp1tid and isOverlapping("V", cluster1, elem, "R1") and cloc2_tid == elem.bp1tid and cluster1.l_orient != cluster1.r_orient:

                match = 1
                elem.count+=1
                AnyMatch = 1
                
                if elem.typeO == "INS_C":
                    elem.typeO = "INS_C_P"
                else:
                    elem.typeO = "INS_C_I_P"

                elem.complete = 1
                    
            # the _P subscript denotes that the breakpoints are confirmed. bp1 is indeed the pasted location for this INS_C.
                
            elif isOverlapping("V", cluster1, elem, "L3") and cloc1_tid == elem.bp3tid and isOverlapping("V", cluster1, elem, "R3") and cloc2_tid == elem.bp3tid and cluster1.l_orient != cluster1.r_orient:

                match = 1
                elem.count+=1
                AnyMatch = 1
                
                if elem.typeO == "INS_C":
                    elem.typeO = "INS_C_P"
                else:
                    elem.typeO = "INS_C_I_P"

                # make bp1 the paste location since now it is known
                # swap bp1 and bp3
                bpList = [elem.bp1_start, elem.bp1_end, elem.bp3_start, elem.bp3_end]
                tidList = [elem.bp1tid, elem.bp3tid]
                swapBPs(bpList, tidList)
                [elem.bp1_start, elem.bp1_end, elem.bp3_start, elem.bp3_end] = bpList
                [elem.bp1tid, elem.bp3tid] = tidList

                elem.complete = 1
                

        if elem.typeO == "TD":

            if (isOverlapping("V", cluster1, elem, "L1") and cloc1_tid == elem.bp1tid) or (isOverlapping("V", cluster1, elem, "R1") and cloc2_tid == elem.bp1tid) or (isOverlapping("V", cluster1, elem, "L2") and cloc1_tid == elem.bp2tid) or (isOverlapping("V", cluster1, elem, "R2") and cloc2_tid == elem.bp2tid):

                match = 1
                elem.count+=1
                AnyMatch = 1
                elem.complete = 1
           

        if elem.typeO == "INS" or elem.typeO == "INS_I":

            if isOverlapping("V", cluster1, elem, "L1") and cloc1_tid == elem.bp1tid and isOverlapping("V", cluster1, elem, "R1")  and cloc2_tid == elem.bp1tid:

                match = 1
                elem.count+=1
                AnyMatch = 1

            # applies to INS_I only
            elif elem.bp3_start == -1:

                if isOverlapping("V", cluster1, elem, "L2") and cloc1_tid == elem.bp2tid and (not isOverlapping("V", cluster1, elem, "R2") or cloc2_tid != elem.bp2tid ) and cluster1.l_orient != cluster1.r_orient:

                    match = 1
                    elem.count+=1
                    AnyMatch = 1
                    [elem.bp3_start, elem.bp3_end] = cloc2_s, cloc2_e
                    if elem.typeO == "INS_I":
                        elem.typeO = "INS_C_I"
                    elif elem.typeO == "INS":
                        elem.typeO = "INS_C"
                        # won't really occur but for future just in case

                elif (not isOverlapping("V", cluster1, elem, "L2") or cloc1_tid != elem.bp2tid) and isOverlapping("V", cluster1, elem, "R2") and cloc2_tid == elem.bp2tid and cluster1.l_orient != cluster1.r_orient:

                    match = 1
                    elem.count+=1
                    AnyMatch = 1
                    [elem.bp3_start, elem.bp3_end] = cloc1_s, cloc1_e
                    if elem.typeO == "INS_I":
                        elem.typeO = "INS_C_I"
                    elif elem.typeO == "INS":
                        elem.typeO = "INS_C"

                elif isOverlapping("V", cluster1, elem, "L1") and cloc1_tid == elem.bp1tid and (not isOverlapping("V", cluster1, elem, "R1") or cloc2_tid != elem.bp1tid):

                    match = 1
                    elem.count+=1
                    AnyMatch = 1
                    [elem.bp3_start, elem.bp3_end] = cloc2_s, cloc2_e

                elif (not isOverlapping("V", cluster1, elem, "L1") or cloc1_tid != elem.bp1tid) and isOverlapping("V", cluster1, elem, "R1") and cloc2_tid == elem.bp1tid:

                    match = 1
                    elem.count+=1
                    AnyMatch = 1
                    [elem.bp3_start, elem.bp3_end] = cloc1_s, cloc1_e
                   

            else:

                # This comparison is conventional and ignores TID. Consistent with cloc1, cloc2 convention of "left" and "right."    
                if elem.bp2_start < elem.bp3_start:
                    
                    if isOverlapping("V", cluster1, elem, "L2") and cloc1_tid == elem.bp2tid and isOverlapping("V", cluster1, elem, "R3") and cloc2_tid == elem.bp3tid and cluster1.l_orient != cluster1.r_orient:
                        match = 1
                else:
                    if isOverlapping("V", cluster1, elem, "L3") and cloc1_tid == elem.bp3tid and isOverlapping("V", cluster1, elem, "R2") and cloc2_tid == elem.bp2tid and cluster1.l_orient != cluster1.r_orient:
                        match = 1
                
                if match:

                    elem.count+=1
                    AnyMatch = 1
                    if elem.typeO == "INS_I":
                        elem.typeO = "INS_C_I"
                    elif elem.typeO == "INS":
                        elem.typeO = "INS_C"

        if elem.typeO == "INV":

            if isOverlapping("V", cluster1, elem, "L1") and cloc1_tid == elem.bp1tid and isOverlapping("V", cluster1, elem, "R2") and cloc2_tid == elem.bp2tid:

                match = 1
                elem.count+=1
                AnyMatch = 1

        if match:

            #print "Cluster", cluster1, "and variant", offset+g, "matched~"
            #print elem.clusterNums, cluster1.mapNum
            elem.count+=1
            if cluster1.mapNum not in elem.clusterNums:
                elem.clusterNums.append(cluster1.mapNum)
                
            graph_m.add_edge(cl1,v1,1)

    if AnyMatch:
        #print "Changing CLAIMED list"
        Claimed[cluster1.mapNum] = 1

    
if __name__ == "__main__":

    # input form files sorted by left TID and position and right TID and position for faster comparison
    # There is double comparison in some cases but that is hard to avoid
    # Remove columns in All_Clusters.txt before sorting and change Cluster class formation code to speed up
    fpL = open("../results/text/All_Clusters_LS.txt","r")
    fpR = open("../results/text/All_Clusters_RS.txt","r")
    fp2 = open("../results/text/doubleOccurences.txt","w")
    fp3 = open("../results/text/All_Variants.txt","w")
    fp4 = open("../results/text/VariantMap.txt", "r")
    fp5 = open("../results/text/ClassifiedVariantMap.txt", "w")
    
    Clusters_L = []
    Clusters_R = []
    OCArray = []
    OCArray_C = []
    tidListL = {}
    tidListR = {}
    tidListL[-1] = 1
    tidListR[-1] = 1 # set keys for unset tid's in variant to facilitate comparisons below
  
    start = 0
    loopCount = 0
    offset = 0
    clusterNum = -1
    Claimed = {}
    fpL.seek(0)
    newBlockL = 0
    newBlockR = 0
    matchGraph = Graph()
    # stores matches
    compareGraph = Graph()
    # This graph stores all comparisons between variants and clusters as undirected edges
    loopNum = 0
    
    # remove -1 last line from All_Clusters in code
    while fpL:

        L_Match = 0
        R_Match = 0
        
        # Both files are same size so will reach end at same time
        lineL = fpL.readline()
        if lineL == "\n" or not (len(lineL) > 0):
            break
       
        if len(lineL) > 0 and lineL[0] == "-":
            continue
        
        clusterNum+= 1
        tempL = Cluster(lineL)
        tempL.number = clusterNum
        lineR = fpR.readline()
        tempR = Cluster(lineR)
        loopCount+=1
        
        if clusterNum % 100 == 0:
            print clusterNum
##            print len(OCArray), "curent variants."
##            print len(Clusters_L), "current left-sorted clusters."
##            print len(Clusters_R), "current right-sorted clusters."
        
        cl_L = 'c' + str(tempL.mapNum)
        cl_R = 'c' + str(tempR.mapNum)

        if tempL.ltid not in tidListL:
            tidListL[tempL.ltid] = 1
        
        if tempR.rtid not in tidListR:
            tidListR[tempR.rtid] = 1

        if len(OCArray) > 0:

            if len(OCArray) % refRate == 0:

                counter=0
                loopNum+=1

                for h,item in enumerate(OCArray):


                    varTIDs = set([item.bp1tid, item.bp2tid, item.bp3tid])

                    # if TIDs in variant have occurred in both cluster lists
                    if varTIDs.issubset(tidListL) and varTIDs.issubset(tidListR):

                        if tempL.ltid not in varTIDs and tempR.rtid not in varTIDs:

                            OCArray_C.append(item)
                            del OCArray[h]
                            #print "Appending", item, "Len:", len(OCArray), "Counter:", counter

                        elif (tempL.l_start - item.bp1_start > MARGIN or tempL.ltid != item.bp1tid) and (tempR.rtid != item.bp1tid or tempR.r_start - item.bp1_start > MARGIN):
                            if (tempL.l_start - item.bp2_start > MARGIN or tempL.ltid != item.bp2tid) and (tempR.rtid != item.bp2tid or tempR.r_start - item.bp2_start > MARGIN):
                                if (tempL.l_start - item.bp3_start > MARGIN or tempL.ltid != item.bp3tid) and (tempR.rtid != item.bp3tid or tempR.r_start - item.bp3_start > MARGIN):

                                    OCArray_C.append(item)
                                    #print "Appending", item, "Len:", len(OCArray), "Counter:", counter
                                    del OCArray[h]
                                

                #for x in OCArray[:counter]:
                    #print "Extracting", x, len(OCArray), "Counter:", counter
            
            compareVariant(tempL, OCArray, Claimed, compareGraph, matchGraph, offset, "L")
            compareVariant(tempR, OCArray, Claimed, compareGraph, matchGraph, offset, "R")
            # Need both so as to compare clusters against later relevant variants too

        if len(Clusters_L) > 0:

            counter=0
	    # $ can start search in middle of Clusters_L if list big to save time?	    
            for item in Clusters_L:

                if item.ltid in tidListL and item.ltid in tidListR:

                    if tempL.ltid != item.ltid and tempR.rtid != item.ltid:

                        counter+=1

                    # Maintaining 2 separate lists and doing separate comparisons may speed up code
                    # Doubling MARGIN below did not change results for small data set
                    elif (tempL.ltid != item.ltid or tempL.l_start - item.l_start > MARGIN) and (tempR.rtid != item.ltid or tempR.r_start - item.l_start > MARGIN):

                        counter+=1
                    else:
                        break

                Clusters_L = Clusters_L[counter:]

        if len(Clusters_R) > 0:
                           
            counter=0    
            for item in Clusters_R:

                if item.rtid in tidListL and item.rtid in tidListR:

                    if tempL.ltid != item.rtid and tempR.rtid != item.rtid:

                        counter+=1

                    # Maintaining 2 separate lists and doing separate comparisons may speed up code
                    # Doubling MARGIN below did not change results for small data set
                    elif (tempL.ltid != item.rtid or tempL.l_start - item.r_start > MARGIN) and (tempR.rtid != item.rtid or tempR.r_start - item.r_start > MARGIN):

                        counter+=1
                    else:
                        break

                Clusters_R = Clusters_R[counter:]
                        
        # This will insure no repeat comparisons and no self-comparisons
 
        compareCluster(tempL, Clusters_L, Claimed, compareGraph, matchGraph, offset, "LL", OCArray)
        compareCluster(tempL, Clusters_R, Claimed, compareGraph, matchGraph, offset, "LR", OCArray)
            
        Clusters_L.append(tempL)
        
        
        compareCluster(tempR, Clusters_L, Claimed, compareGraph, matchGraph, offset, "RL", OCArray)
        compareCluster(tempR, Clusters_R, Claimed, compareGraph, matchGraph, offset, "RR", OCArray)

        Clusters_R.append(tempR)
     
    fpL.seek(0)

    hashedVM = {}

    OCArray_C = OCArray_C + OCArray
    #print len(OCArray)
    #for x in OCArray:
        #OCArray_C.append(x)
        #print "At end. Appending", x
        #print "At end. Extracting", x
    OCArray = []

    for line in fp4:
 
        line_split = line.split()
        clNum = int(line_split[0])
        clList = line_split[1:]
        hashedVM[clNum] = clList
        #print clNum
        
    writeVariants(OCArray_C, fp3, fp5, hashedVM, 0)
    
    counter = -1
    TDStore = []
    TDArtefacts = []
    varCount = len(OCArray_C)
    OCArray_C = []
    store = 0
    
    for line in fpL:

        counter+= 1
        

        if line != "\n" and len(line) > 0:

         cluster = Cluster(line)
         #print "Unclaimed", counter
            
         if not Claimed.has_key(cluster.mapNum):
             
            store =1
            temp = OverlappingCluster()
            temp.bp1tid = cluster.ltid
            temp.bp2tid = cluster.rtid
            temp.bp1_start = cluster.l_start
            temp.bp1_end = cluster.l_end
            temp.bp2_start = cluster.r_start
            temp.bp2_end = cluster.r_end
            temp.clusterNums.append(cluster.mapNum)

            #$ These insertions will be written twice but bedtools comparison will not fault that. Clean up for other comparisons and diploid comparisons.
            if cluster.r_orient == "2":

                temp.typeO = "INS_R"
                temp.bp2_start = temp.bp1_start
                temp.bp2_end = temp.bp1_end

            # TDs will be unclaimed if no overlap occurs, which is fine -- no double counting
            elif cluster.l_orient == 1 and cluster.r_orient ==0 and cluster.l_start < cluster.r_start and cluster.ltid == cluster.rtid:
                temp.typeO = "TD"

                # For medium-sized TDs, avoid double counting from overlapping mate clusters. For small, only below will figure.
                TDStore.append(temp)

	    # $ remove -- test for comparison with other tools	
	    #elif cluster.l_orient == cluster.r_orient and cluster.ltid == cluster.rtid:
		#temp.typeO = "INV"
 
            elif cluster.l_orient == 0 and cluster.r_orient ==1 and cluster.l_start < cluster.r_end and cluster.ltid == cluster.rtid:
                temp.typeO = "DEL"
                #print cluster
                
            # artefact TD cluster from overlapping reads for v small TDs, addressed later   
            elif cluster.l_start > cluster.r_end and cluster.l_orient == 0 and cluster.r_orient == 1 and cluster.ltid == cluster.rtid:
                store=0
                temp.typeO = "TD"
                bpList = [temp.bp1_start, temp.bp1_end, temp.bp2_start, temp.bp2_end]
                tidList = [temp.bp1tid, temp.bp2tid]
                swapBPs(bpList, tidList)
                [temp.bp1_start, temp.bp1_end, temp.bp2_start, temp.bp2_end] = bpList
                [temp.bp1tid, temp.bp2tid] =  tidList
                
                TDArtefacts.append(temp)
                
            else:
                temp.typeO = "Unknown"

            #print store, temp

            if store:

                OCArray_C.append(temp)

    # May change results by a little but, but not much
    for elem in TDArtefacts:

        storeTD = 1

        for TD in TDStore:

            if abs(elem.bp1_start - TD.bp1_start) < MARGIN and abs(elem.bp2_start - TD.bp2_start) < MARGIN and TD.bp1tid == elem.bp1tid:
                storeTD = 0
                break

        if storeTD:
            OCArray_C.append(elem)

    writeVariants(OCArray_C, fp3, fp5, hashedVM, varCount)
    fpL.close()
    fp2.close()
    fp3.close()
    fp4.close()
    fp5.close()


        


# $ To speed up or improve

# Speed up by checking if cl1_end - cl2_start > 2(mu + 3 sigma IL) then don't call isOverlap fxn at all

# Remove columns from All_Clusters
# Keep different variant possibilities within variant, i.e. check with sequential if stms if bp1 matches, then if bp2 also matches, else new variant etc. till 3 bps formed for INS or max for others.

# Revise variant bounds when 3rd cluster matches

# can keep 3 dynamically sorted lists (each entry goes in sorted) for 3 bp's in variant and compare with incoming cluster to speed up further

# can call some variants complete (3 bps etc.) and move them to another list. No new clusters here will change designation of variant. May leave unclaimed clusters? Not significantly.
# Can do cluster comparisons before variant comparisons to skip comparing with variants whose clusters have been compared against
# Sub-sort cluster list and use all 4 as criteria to refresh variant list and cluster lists?
# streamline iF conditions for cluster comparison etc?

# Exception 2: IMP: Allow for variant branches to be formed. A variant can have 2 clusters matching, and 2 different third clusters yielding two different variants, e.g. with small insertions.
# Temporary/present fix: If 2 clusters match same cluster (after checking for if share variant), then don't compare them.
# OR perhaps, just postprocess to remove 2 bp insertions (and others) that are subsets of 3 bp ins (etc), within SLOP, at end -- this way there will be no extras or less!
# Exception 1: Variant complete, but 2 complete variants with same bp's will be IDd in such a case. Cl3 and Cl4 can miss each other in small insertion (draw it out with 4 clusters coming 1 by 1)
# If call all these diploid and true, b/c bp's do not match, will gain some in right calls but lose some in these wrong diploid calls.
# Other cases have been thought through.
# Diff w/ v9: order of cluster comparison different, so FPs can be different even when break used in both. A few less FPs and few less TPs here with small data set.
#$ Add condition to check if contained fully within
        






                







                

        

