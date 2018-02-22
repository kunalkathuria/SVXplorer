#!/usr/bin/env python

# Form variants from clusters output by formPEClusters.py, saved in allClusters.txt, i.e. "cluster consolidation."

import argparse
import logging

class Vertex:
    def __init__(self, node):
        self.id = node
        self.adjacent = {}
    def __str__(self):
        return str(self.id) + ' adjacent: ' + str([x.id for x in self.adjacent])
    def add_neighbor(self, neighbor, weight=0):
        self.adjacent[neighbor] = weight
    def get_connections(self):
        return self.adjacent.keys()
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

class clusterI(object):
    # "imported" cluster object read from allClusters.txt
    def __init__(self, data):
        data_split = data.split()
        self.l_start = int(data_split[4])
        self.l_end = int(data_split[5])
        self.r_start = int(data_split[7])
        self.r_end = int(data_split[8])
        self.l_orient = int(data_split[2][0])
        self.r_orient = int(data_split[2][1])
        self.typeC = data_split[2]
        self.lTID = data_split[3]
        self.rTID = data_split[6]
        self.mapNum = int(data_split[0])
        self.isSmall = int(data_split[9])
    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.mapNum, self.l_start,self.l_end, self.r_start, self.r_end, self.typeC, self.lTID, self.rTID, self.l_orient, self.r_orient)

class consCluster(object):
    # consolidated cluster or complex SV
    def __init__(self):
        self.clusterNums = [] # original cluster numbers that support this variant
        self.support = "PE"
        self.complete = 0 # All clusters that can match with this variant have arrived and matched
        self.variantNum = None
        # at most 4 clusters will overlap to produce a maximum of 3 breakpts in our regime
        self.bp1_start = -1
        self.bp2_start = -1
        self.bp3_start = -1
        self.bp1_end = -1
        self.bp2_end = -1
        self.bp3_end = -1
        self.count = 2 # 2 clusters overlap to initiate consolidated cluster
        self.SVType = None
        # it is an insertion cluster if overlapping point has had both orientation of reads
        self.bp1_hasAlmtF = 0
        self.bp1_hasAlmtR = 0
        self.bp1_orient = -1
        self.bp2_orient = -1
        self.bp3_orient = -1
        self.bp1TID = -1
        self.bp2TID = -1
        self.bp3TID = -1

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.SVType, self.bp1TID, self.bp1_start, self.bp1_end, self.bp2TID , self.bp2_start, self.bp2_end, self.bp3TID, self.bp3_start, self.bp3_end, self.support, self.count)

def isOverlapping(CorV, clus1, clus2, LR, slop):
    if CorV == "C":
        if LR=="LL":
            if not (clus2.l_start > clus1.l_end or clus2.l_end < clus1.l_start):
                return 1
        elif LR=="LR":
            if not (clus2.r_start > clus1.l_end or clus2.r_end < clus1.l_start):
                return 1
        elif LR=="RL":
            if not (clus2.l_start > clus1.r_end or clus2.l_end < clus1.r_start):
                return 1
        elif LR=="RR":
            if not (clus2.r_start > clus1.r_end or clus2.r_end < clus1.r_start):
                return 1

    elif CorV == "V":
        if LR=="L1":
            if not (clus1.l_start > clus2.bp1_end or clus1.l_end < clus2.bp1_start):
                return 1
        elif LR=="L2":
            if not (clus1.l_start > clus2.bp2_end or clus1.l_end < clus2.bp2_start):
                return 1
        elif LR=="L3":
            if not (clus1.l_start > clus2.bp3_end or clus1.l_end < clus2.bp3_start):
                return 1
        elif LR=="R1":
            if not (clus1.r_start > clus2.bp1_end or clus1.r_end < clus2.bp1_start):
                return 1
        elif LR=="R2":
            if not (clus1.r_start > clus2.bp2_end or clus1.r_end < clus2.bp2_start):
                return 1
        elif LR=="R3":
            if not (clus1.r_start > clus2.bp3_end or clus1.r_end < clus2.bp3_start):
                return 1
    return 0

def setBPs(clus1, clus2, LR):
    if LR == "LL":
        sorted_bps = sorted([clus2.l_start, clus2.l_end, clus1.l_start, clus1.l_end])
        return sorted_bps[1], sorted_bps[2]
    elif LR == "LR":
        sorted_bps = sorted([clus1.l_start, clus1.l_end, clus2.r_start, clus2.r_end])
        return sorted_bps[1], sorted_bps[2]
    elif LR == "RL":
        sorted_bps = sorted([clus1.r_start, clus1.r_end, clus2.l_start, clus2.l_end])
        return sorted_bps[1], sorted_bps[2]
    elif LR == "RR":
        sorted_bps = sorted([clus2.r_start, clus2.r_end, clus1.r_start, clus1.r_end])
        return sorted_bps[1], sorted_bps[2]
    return -1,-1

def formCutPasteINS(newVariant, cluster1, clusterP, LR):
    newSVFlag=1
    # if middle breakpoint is overlapping, then has to be a cut insertion
    lBP_end = min(newVariant.bp2_end, newVariant.bp3_end)
    rBP_start = max(newVariant.bp2_start, newVariant.bp3_start)

    if (LR == "RL" or LR == "LR") and newVariant.bp1TID == newVariant.bp2TID == newVariant.bp3TID:
        if lBP_end <  newVariant.bp1_start < newVariant.bp1_end < rBP_start:
            if newVariant.SVType == "INS":
                newVariant.SVType = "INS_C"
            elif newVariant.SVType == "INS_I":
                newVariant.SVType = "INS_C_I"
            newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, \
                newVariant.bp2_end = newVariant.bp2_start, newVariant.bp2_end, \
                newVariant.bp1_start, newVariant.bp1_end
            newVariant.bp2_orient, newVariant.bp1_orient = \
                newVariant.bp2_orient, newVariant.bp1_orient

            # make bp1 < bp3: this is conventional and consistent
            if newVariant.bp1_start > newVariant.bp3_start:
                newVariant.bp1_start, newVariant.bp1_end, newVariant.bp3_start,\
                    newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
                    newVariant.bp1_start, newVariant.bp1_end
                newVariant.bp1_orient, newVariant.bp2_orient = newVariant.bp2_orient,\
                    newVariant.bp1_orient
    # if chr where overlap occurs also is the site for 1 of the other mates/half clusters
    elif (LR == "RR" or LR == "LL") and newVariant.bp1TID == newVariant.bp2TID and \
        newVariant.bp2TID != newVariant.bp3TID:
        if cluster1.lTID == cluster1.rTID and \
            not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
            return 0
        if clusterP.lTID == clusterP.rTID and \
            not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
            return 0
        logging.debug('Orientation matches for INS_C')
        # breakpoint locations are confirmed
        if newVariant.SVType == "INS":
            newVariant.SVType = "INS_C_P"
        elif newVariant.SVType == "INS_I":
            newVariant.SVType = "INS_C_I_P"
        # swap bp1 and bp3, as bp1 is paste location by convention
        newVariant.bp1_start, newVariant.bp1_end, newVariant.bp3_start, \
            newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
            newVariant.bp1_start, newVariant.bp1_end
        newVariant.bp1TID, newVariant.bp3TID = newVariant.bp3TID, newVariant.bp1TID
        newVariant.bp1_orient, newVariant.bp3_orient = \
            newVariant.bp3_orient, newVariant.bp1_orient

        # make bp2 < bp3: this is conventional and consistent
        if newVariant.bp2_start > newVariant.bp3_start:
            newVariant.bp2_start, newVariant.bp2_end, newVariant.bp3_start,\
                newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
                newVariant.bp2_start, newVariant.bp2_end
            newVariant.bp3_orient, newVariant.bp2_orient = newVariant.bp2_orient,\
                newVariant.bp3_orient
    elif  (LR == "RR" or LR == "LL") and newVariant.bp1TID == newVariant.bp3TID and \
        newVariant.bp3TID != newVariant.bp2TID:
        if cluster1.lTID == cluster1.rTID and \
            not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
            return 0
        if clusterP.lTID == clusterP.rTID and \
            not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
            return 0
        logging.debug('Orientation matches for INS_C')
        if newVariant.SVType == "INS":
            newVariant.SVType = "INS_C_P"
        elif newVariant.SVType == "INS_I":
            newVariant.SVType = "INS_C_I_P"
        # swap bp1 and bp2, as bp1 is paste location by convention
        newVariant.bp1TID, newVariant.bp2TID = newVariant.bp2TID, newVariant.bp1TID
        newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, \
            newVariant.bp2_end = newVariant.bp2_start, newVariant.bp2_end, \
            newVariant.bp1_start, newVariant.bp1_end
        newVariant.bp1_orient, newVariant.bp2_orient = \
            newVariant.bp2_orient, newVariant.bp1_orient

        # make bp2 < bp3: this is conventional and consistent
        if newVariant.bp2_start > newVariant.bp3_start:
            newVariant.bp2_start, newVariant.bp2_end, newVariant.bp3_start,\
                newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
                newVariant.bp2_start, newVariant.bp2_end
            newVariant.bp3_orient, newVariant.bp2_orient = newVariant.bp2_orient,\
                newVariant.bp3_orient
    # bp2 < bp3 for regular copy-paste INS on same chr by convention
    elif (LR == "LL" or LR == "RR") and newVariant.bp2_start > newVariant.bp3_start:
        newVariant.bp2_start, newVariant.bp2_end, newVariant.bp3_start,\
        newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
        newVariant.bp2_start, newVariant.bp2_end
        newVariant.bp2_orient, newVariant.bp3_orient = newVariant.bp3_orient, newVariant.bp1_orient

    # note: if bp1 (overlapping bp) is on different chr from other 2 then could be cut or copy
    # in this case, don't change SV type till further evidence seen

def writeVariants(varList, fpAV, fpCVM, hashedVM, offset):
    for g, item in enumerate(varList):
        num = offset + 1 + g
        fpCVM.write("%s" %(num))
        if item.SVType == "INS" or item.SVType == "INS_I" or \
            item.SVType == "INS_C_P" or item.SVType == "INS_C_I_P":
                # just to be safe, but should be ensured anyway
                if item.bp2TID == item.bp3TID and item.bp3_start < item.bp2_start:
                    item.bp2_start, item.bp2_end, item.bp3_start, item.bp3_end=\
                        item.bp3_start, item.bp3_end, item.bp2_start, item.bp2_end
        elif item.SVType.find("TD") != -1:
            if item.bp1TID == item.bp2TID and item.bp1_start > item.bp2_start:
                item.bp1_start, item.bp1_end, item.bp2_start, item.bp2_end=\
                    item.bp2_start, item.bp2_end, item.bp1_start, item.bp1_end

        fpAV.write("%s\t%s\n" %(num,item))
        for elem in item.clusterNums:
            try:
                for elem2 in hashedVM[elem]:
                    #print "Writing", elem2, "from cluster", elem
                    fpCVM.write("\t%s" %elem2)
            except:
                print "Exception writing in Variant Map from cluster", elem
                exit(1)
        fpCVM.write("\n")

def compareCluster(cluster1, clusters, claimedCls, graph_c, graph_m, LR, consolidatedCls, slop, RDL_Factor, RDL):
    anyMatch = 0
    cl1 = 'c' + str(cluster1.mapNum)
    logging.debug('Start loop to compare this cluster to all clusters in list')
    for clusterP in clusters:
        if cluster1.mapNum == clusterP.mapNum:
            continue
        cl2 = 'c' + str(clusterP.mapNum)
        logging.debug("Comparing cluster %s and cluster %s as %s", cl1, cl2, LR)
        # If 2 clusters have been compared do not compare them again
        if graph_c.get_vertex(cl1) != None:
            if cl2 in graph_c.get_vertex(cl1).get_connections():
                logging.debug('Continue as these 2 clusters were compared before')
                continue
        # If 2 clusters share variant, then comparing them is redundant.
        # Any subtleties can be addressed when comparing to existing variant, e.g. variant branches.
        skip = 0
        if graph_m.get_vertex(cl1) != None and graph_m.get_vertex(cl2) != None:
            for variant in graph_m.get_vertex(cl1).get_connections():
                if variant in graph_m.get_vertex(cl2).get_connections():
                    skip = 1
                    break
        if skip:
            logging.debug('Continue as these 2 clusters already part of same consolidated cluster')
            continue

        # determine which sides of clusters overlap
        logging.debug('Determining signature of 2-cluster overlap: left with left (LL) etc.')
        LLOverlap=0
        LROverlap=0
        RLOverlap=0
        RROverlap=0
        if cluster1.lTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "LL", slop):
            LLOverlap = 1
        if cluster1.lTID == clusterP.rTID and isOverlapping("C", cluster1, clusterP, "LR", slop):
            LROverlap = 1
        if cluster1.rTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "RL", slop):
            RLOverlap = 1
        if cluster1.rTID == clusterP.rTID and isOverlapping("C", cluster1, clusterP, "RR", slop):
            RROverlap = 1
        if (LR == "LL" and not LLOverlap) or (LR == "LR" and not LROverlap) \
            or (LR =="RL" and not RLOverlap) or (LR == "RR" and not RROverlap):
            continue
        if not LLOverlap and not LROverlap and not RLOverlap and not RROverlap:
            continue

        # initialize all variables
        graph_c.add_edge(cl1, cl2, 1)
        newSVFlag = 0
        newVariant = consCluster()
        newVariant.SVType = None

        if cluster1.r_orient == 2 and clusterP.r_orient == 2 and LLOverlap \
            and cluster1.l_orient != clusterP.l_orient:
            logging.debug('Tagged as De Novo INS')

            newVariant.SVType = "INS_U"
            newSVFlag=1
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LL")
            newVariant.bp1TID = cluster1.lTID

        elif LLOverlap and RROverlap and cluster1.l_orient != clusterP.l_orient and \
            cluster1.r_orient != clusterP.r_orient and cluster1.r_orient == cluster1.l_orient \
            and cluster1.lTID == cluster1.rTID:
            logging.debug('Tagged as INV')

            newVariant.SVType = "INV"
            newSVFlag=1
            # bp1 is the left alignment
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LL")
            newVariant.bp2_start, newVariant.bp2_end = setBPs(cluster1, clusterP, "RR")
            newVariant.bp1TID, newVariant.bp2TID  = cluster1.lTID, cluster1.rTID

        # small insertions
        elif LROverlap and RROverlap and not LLOverlap and not RLOverlap and cluster1.isSmall \
            and not clusterP.isSmall:
            logging.debug('Tagged as Small INS')
            # strictly, this should be TD_I_I but we are grouping them into INS_I
            if clusterP.l_orient == clusterP.r_orient:
                newVariant.SVType = "INS_I"
            elif not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID) \
                or (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                newVariant.SVType = "INS"
            else:
                # TD or INS
                newVariant.SVType = "TD_I"
            newSVFlag=1
            # only bp_1 numbering is crucial
            # not nec--can make more precise by evaluating which of 3 bp pairs yields least bp width
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LR")
            newVariant.bp2_start, newVariant.bp2_end = clusterP.l_start, clusterP.l_end
            newVariant.bp1_hasAlmtF = not clusterP.r_orient
            newVariant.bp1_hasAlmtR = clusterP.r_orient
            newVariant.bp1TID, newVariant.bp2TID = cluster1.lTID, clusterP.lTID

        elif RLOverlap and RROverlap and not LLOverlap and not LROverlap and \
            clusterP.isSmall and not cluster1.isSmall:

            if cluster1.l_orient == cluster1.r_orient:
                newVariant.SVType = "INS_I"
            elif not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID)\
                or (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                newVariant.SVType = "INS"
            else:
                newVariant.SVType = "TD_I"
            newSVFlag=1
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "RL")
            newVariant.bp2_start, newVariant.bp2_end = cluster1.l_start, cluster1.l_end
            newVariant.bp1_hasAlmtF = not cluster1.r_orient
            newVariant.bp1_hasAlmtR = cluster1.r_orient
            newVariant.bp1TID = cluster1.rTID
            newVariant.bp2TID = cluster1.lTID

        elif LLOverlap and RLOverlap and not LROverlap and not RROverlap and \
            cluster1.isSmall and not clusterP.isSmall:

            if clusterP.l_orient == clusterP.r_orient:
                newVariant.SVType = "INS_I"
            elif not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID)\
                or (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                newVariant.SVType = "INS"
            else:
                newVariant.SVType = "TD_I"
            newSVFlag=1
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "RL")
            newVariant.bp2_start, newVariant.bp2_end = clusterP.r_start, clusterP.r_end
            newVariant.bp1_hasAlmtF = not clusterP.l_orient
            newVariant.bp1_hasAlmtR = clusterP.l_orient
            newVariant.bp1TID = cluster1.lTID
            newVariant.bp2TID = clusterP.rTID

        elif LLOverlap and LROverlap and not RLOverlap and not RROverlap and \
            clusterP.isSmall and not cluster1.isSmall:

            if cluster1.l_orient == cluster1.r_orient:
                newVariant.SVType = "INS_I"
            elif not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID)\
                or (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                newVariant.SVType = "INS"
            else:
                newVariant.SVType = "TD_I"
            newSVFlag=1
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LR")
            newVariant.bp2_start, newVariant.bp2_end = cluster1.r_start, cluster1.r_end
            newVariant.bp1_hasAlmtF = not cluster1.l_orient
            newVariant.bp1_hasAlmtR = cluster1.l_orient
            newVariant.bp1TID = cluster1.lTID
            newVariant.bp2TID = cluster1.rTID

        # Large insertions
        elif LLOverlap and cluster1.l_orient != clusterP.l_orient and (cluster1.lTID == \
            cluster1.rTID or cluster1.lTID == clusterP.rTID or cluster1.rTID == clusterP.rTID):
            logging.debug('Large INS check 1: left with left overlap')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.l_end < cluster1.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.l_end < clusterP.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2')
                continue
            if cluster1.l_orient == cluster1.r_orient and clusterP.l_orient == clusterP.r_orient:
                newVariant.SVType = "INS_I"
            elif cluster1.l_orient != cluster1.r_orient and clusterP.l_orient != clusterP.r_orient:
                newVariant.SVType = "INS"
                if cluster1.rTID == clusterP.rTID and cluster1.r_end < clusterP.r_start and \
                    not (cluster1.r_orient == 1 and clusterP.r_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:1")
                    continue
                if cluster1.rTID == clusterP.rTID and clusterP.r_end < cluster1.r_start and \
                    not (clusterP.r_orient == 1 and cluster1.r_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:2")
                    continue

            if newVariant.SVType == "INS" or newVariant.SVType == "INS_I":
                newSVFlag=1
                # set breakpoints
                newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LL")
                newVariant.bp2_start, newVariant.bp2_end = cluster1.r_start, cluster1.r_end
                newVariant.bp3_start, newVariant.bp3_end = clusterP.r_start, clusterP.r_end
                newVariant.bp1TID, newVariant.bp2TID, newVariant.bp3TID = \
                cluster1.lTID, cluster1.rTID, clusterP.rTID
                newVariant.bp2_orient, newVariant.bp3_orient = cluster1.r_orient, clusterP.r_orient
                # check for cut-paste
                logging.debug('Large INS check 1: check if could be cut-paste')
                isValidCP = formCutPasteINS(newVariant, cluster1, clusterP, "LL")
                if isValidCP == 0:
                    continue

        elif RROverlap and cluster1.r_orient != clusterP.r_orient and \
            (cluster1.rTID == cluster1.lTID or cluster1.rTID == clusterP.lTID or cluster1.lTID == clusterP.lTID):
            logging.debug('Large INS check 2: right with right overlap')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.r_start > cluster1.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.r_start > clusterP.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2')
                continue
            if cluster1.l_orient == cluster1.r_orient and clusterP.l_orient == clusterP.r_orient:
                newVariant.SVType = "INS_I"
            elif cluster1.l_orient != cluster1.r_orient and clusterP.l_orient != clusterP.r_orient:
                newVariant.SVType = "INS"
                if cluster1.lTID == clusterP.lTID and cluster1.l_end < clusterP.l_start and \
                    not (cluster1.l_orient == 1 and clusterP.l_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:1")
                    continue
                if cluster1.lTID == clusterP.lTID and clusterP.l_end < cluster1.l_start and \
                    not (clusterP.l_orient == 1 and cluster1.l_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:2")
                    continue

            if newVariant.SVType == "INS" or newVariant.SVType == "INS_I":
                newSVFlag=1
                # set breakpoints
                newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "RR")
                newVariant.bp2_start, newVariant.bp2_end = cluster1.l_start, cluster1.l_end
                newVariant.bp3_start, newVariant.bp3_end = clusterP.l_start, clusterP.l_end
                newVariant.bp1TID, newVariant.bp2TID, newVariant.bp3TID = \
                cluster1.rTID, cluster1.lTID, clusterP.lTID
                newVariant.bp2_orient, newVariant.bp3_orient = cluster1.l_orient, clusterP.l_orient
                # check for cut-paste
                logging.debug('Large INS check 2: check if could be cut-paste')
                isValidCP = formCutPasteINS(newVariant, cluster1, clusterP, "RR")
                if isValidCP == 0:
                    continue
        
        elif LROverlap and not RLOverlap and (cluster1.rTID == cluster1.lTID or \
            clusterP.lTID == cluster1.lTID or cluster1.rTID == clusterP.lTID) and \
            cluster1.l_orient != clusterP.r_orient:
            logging.debug('Large INS check 3: left mate of cluster 1 overlapping with right of 2')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.r_end < cluster1.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.l_start > clusterP.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2')
                continue
            if cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID and \
                not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                continue
            newVariant.SVType = "INS_C"
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LR")
            newVariant.bp2_start, newVariant.bp2_end = cluster1.r_start, cluster1.r_end
            newVariant.bp3_start, newVariant.bp3_end = clusterP.l_start, clusterP.l_end
            newVariant.bp1TID, newVariant.bp2TID, newVariant.bp3TID \
                = cluster1.lTID, cluster1.rTID, clusterP.lTID
            newVariant.bp2_orient, newVariant.bp3_orient = cluster1.r_orient, clusterP.l_orient
            newSVFlag = 1
            # one mate same orientation, other opposite, then inverted cut insertion
            if (cluster1.l_orient == cluster1.r_orient and clusterP.l_orient != clusterP.r_orient) \
                or (clusterP.r_orient == clusterP.l_orient and cluster1.l_orient != cluster1.r_orient):
                newVariant.SVType = "INS_C_I"
            #regular insertion if both mates show inversion
            if newVariant.bp1TID != newVariant.bp2TID and newVariant.bp1TID != newVariant.bp3TID \
                and cluster1.l_orient == cluster1.r_orient and clusterP.r_orient == clusterP.l_orient:
                newVariant.SVType = "INS_I"
            # as always if paste location is on diff chr, it is regular INS
            elif newVariant.bp1TID != newVariant.bp2TID and newVariant.bp1TID != newVariant.bp3TID:
                    newVariant.SVType = "INS"
            else:
                logging.debug('Large INS check 3: check for cut-paste')
                formCutPasteINS(newVariant, cluster1, clusterP, "LR")

                logging.debug('Large INS check 3: mark as cut-paste if one of source breakpoints is on same chr as overlap pt')
                if newVariant.bp1TID == newVariant.bp2TID and newVariant.bp3TID != newVariant.bp1TID:
                    if cluster1.lTID == cluster1.rTID and \
                        not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
                        continue
                    if clusterP.lTID == clusterP.rTID and \
                        not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                        continue
                    logging.debug('Orientation matches for INS_C')
                    # breakpoint locations are confirmed
                    if newVariant.SVType == "INS_C":
                        newVariant.SVType = "INS_C_P"
                    elif newVariant.SVType == "INS_C_I":
                        newVariant.SVType = "INS_C_I_P"

                    # Put 2 and 3 on same chromosome, as 1 is pasted location in our convention
                    newVariant.bp1_start, newVariant.bp1_end, newVariant.bp3_start, newVariant.bp3_end=\
                        newVariant.bp3_start, newVariant.bp3_end, newVariant.bp1_start, newVariant.bp1_end
                    newVariant.bp1TID, newVariant.bp3TID = newVariant.bp3TID, newVariant.bp1TID
                    newVariant.bp1_orient, newVariant.bp3_orient = newVariant.bp3_orient, newVariant.bp1_orient

                elif newVariant.bp1TID == newVariant.bp3TID and newVariant.bp2TID != newVariant.bp1TID:
                    if cluster1.lTID == cluster1.rTID and \
                        not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
                        continue
                    if clusterP.lTID == clusterP.rTID and \
                        not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                        continue
                    logging.debug('Orientation matches for INS_C')
                    if newVariant.SVType == "INS_C":
                        newVariant.SVType = "INS_C_P"
                    elif newVariant.SVType == "INS_C_I":
                        newVariant.SVType = "INS_C_I_P"

                    # Put 2 and 3 on same chromosome, as 1 is pasted location in our convention
                    newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, newVariant.bp2_end=\
                        newVariant.bp2_start, newVariant.bp2_end, newVariant.bp1_start, newVariant.bp1_end
                    newVariant.bp1TID, newVariant.bp2TID = newVariant.bp2TID, newVariant.bp1TID
                    newVariant.bp1_orient, newVariant.bp2_orient = newVariant.bp2_orient, newVariant.bp1_orient

            # make bp2 < bp3 for C_p and regular INS: this is conventional and consistent
            if (newVariant.SVType == "INS_C_P" or newVariant.SVType == "INS_C_I_P" or \
                newVariant.SVType == "INS" or newVariant.SVType == "INS_I") and \
                newVariant.bp2_start > newVariant.bp3_start:

                newVariant.bp2_start, newVariant.bp2_end, newVariant.bp3_start,\
                    newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
                    newVariant.bp2_start, newVariant.bp2_end
                newVariant.bp3_orient, newVariant.bp2_orient = newVariant.bp2_orient,\
                    newVariant.bp3_orient

        elif RLOverlap and not LROverlap and (cluster1.rTID == cluster1.lTID or \
            clusterP.rTID == cluster1.rTID or cluster1.lTID == clusterP.rTID) and \
            cluster1.r_orient != clusterP.l_orient:
            logging.debug('Large INS check 4: right mate of cluster 1 overlapping with left mate of cluster 2')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.l_start > cluster1.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.r_end < clusterP.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2 %d %d',cluster1.r_end,clusterP.r_end)
                continue
            if cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID and \
                not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
                continue
            newVariant.SVType = "INS_C"
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "RL")
            newVariant.bp2_start, newVariant.bp2_end = cluster1.l_start, cluster1.l_end
            newVariant.bp3_start, newVariant.bp3_end = clusterP.r_start, clusterP.r_end
            newVariant.bp1TID, newVariant.bp2TID, newVariant.bp3TID \
                = cluster1.rTID, cluster1.lTID, clusterP.rTID
            newVariant.bp2_orient, newVariant.bp3_orient = cluster1.l_orient, clusterP.r_orient
            newSVFlag = 1
            # one mate same orientation, other opposite, then inverted cut insertion
            if (cluster1.l_orient == cluster1.r_orient and clusterP.l_orient != clusterP.r_orient) \
                or (clusterP.r_orient == clusterP.l_orient and cluster1.l_orient != cluster1.r_orient):
                newVariant.SVType = "INS_C_I"
            #regular insertion if both mates show inversion
            if newVariant.bp1TID != newVariant.bp2TID and newVariant.bp1TID != newVariant.bp3TID \
                and cluster1.l_orient == cluster1.r_orient and clusterP.r_orient == clusterP.l_orient:
                newVariant.SVType = "INS_I"
            # as always if paste location is on diff chr, it is regular INS
            elif newVariant.bp1TID != newVariant.bp2TID and newVariant.bp1TID != newVariant.bp3TID:
                    newVariant.SVType = "INS"
            else:
                formCutPasteINS(newVariant, cluster1, clusterP, "RL")
                
                logging.debug('Large INS check 4: mark as cut-paste if one of source breakpoints is on same chr as overlap pt')
                if newVariant.bp1TID == newVariant.bp2TID and newVariant.bp3TID != newVariant.bp1TID:
                    logging.debug('bp3 chr diff')
                    if cluster1.lTID == cluster1.rTID and \
                        not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
                        continue
                    if clusterP.lTID == clusterP.rTID and \
                        not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                        continue
                    # breakpoint locations are confirmed
                    if newVariant.SVType == "INS_C":
                        newVariant.SVType = "INS_C_P"
                    elif newVariant.SVType == "INS_C_I":
                        newVariant.SVType = "INS_C_I_P"

                    # Put 2 and 3 on same chromosome, as 1 is pasted location in our convention
                    newVariant.bp1_start, newVariant.bp1_end, newVariant.bp3_start, newVariant.bp3_end=\
                        newVariant.bp3_start, newVariant.bp3_end, newVariant.bp1_start, newVariant.bp1_end
                    newVariant.bp1TID, newVariant.bp3TID = newVariant.bp3TID, newVariant.bp1TID
                    newVariant.bp1_orient, newVariant.bp3_orient = newVariant.bp3_orient, newVariant.bp1_orient

                elif newVariant.bp1TID == newVariant.bp3TID and newVariant.bp2TID != newVariant.bp1TID:
                    logging.debug('bp2 chr diff')
                    if cluster1.lTID == cluster1.rTID and \
                        not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
                        continue
                    if clusterP.lTID == clusterP.rTID and \
                        not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                        continue

                    if newVariant.SVType == "INS_C":
                        newVariant.SVType = "INS_C_P"
                    elif newVariant.SVType == "INS_C_I":
                        newVariant.SVType = "INS_C_I_P"
                    logging.debug('%d %d %d %d', newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, newVariant.bp2_end)
                    # Put 2 and 3 on same chromosome, as 1 is pasted location in our convention
                    newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, newVariant.bp2_end=\
                        newVariant.bp2_start, newVariant.bp2_end, newVariant.bp1_start, newVariant.bp1_end
                    newVariant.bp1TID, newVariant.bp2TID = newVariant.bp2TID, newVariant.bp1TID
                    newVariant.bp1_orient, newVariant.bp2_orient = newVariant.bp2_orient, newVariant.bp1_orient

            # make bp2 < bp3: this is conventional and consistent
            if (newVariant.SVType == "INS_C_P" or newVariant.SVType == "INS_C_I_P" or \
                newVariant.SVType == "INS" or newVariant.SVType == "INS_I") and \
                newVariant.bp2_start > newVariant.bp3_start:
                logging.debug('bp2 > bp3')
                newVariant.bp2_start, newVariant.bp2_end, newVariant.bp3_start,\
                    newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
                    newVariant.bp2_start, newVariant.bp2_end
                newVariant.bp3_orient, newVariant.bp2_orient = newVariant.bp2_orient,\
                    newVariant.bp3_orient

        # if clusters matched to form new SV
        if newSVFlag:
            logging.debug('Writing new complex variant:%s', newVariant.SVType)
            newVariant.clusterNums.append(cluster1.mapNum)
            newVariant.clusterNums.append(clusterP.mapNum)
            # should never occur after start, so okay to refresh variant buffer
            if len(consolidatedCls) > 0:
                newVariant.variantNum = consolidatedCls[-1].variantNum
            else:
                newVariant.variantNum = 1
            consolidatedCls.append(newVariant)
            claimedCls.add(clusterP.mapNum)
            anyMatch = 1
            v1 = 'v' + str(newVariant.variantNum)
            graph_m.add_edge(cl2,v1,1)
            graph_m.add_edge(cl1,v1,1)
            graph_c.add_edge(cl2,v1,1)
            graph_c.add_edge(cl1,v1,1)
            #print "Match:", v1, " from", cluster1, "and", clusterP

    if anyMatch:
        logging.debug('Adding cluster to claimed list')
        claimedCls.add(cluster1.mapNum)

def compareVariant(cluster1, varList, claimedCls, graph_c, graph_m, LR, maxClCombGap, slop):
    anyMatch = 0
    cl1 = 'c' + str(cluster1.mapNum)
    logging.debug('Start loop to compare this cluster to all complex variants in list')
    for g,elem in enumerate(varList):
        match = 0
        v1 = 'v' + str(elem.variantNum)

        # If cluster has been compared to variant, do not compare
        logging.debug('Do not compare cluster to variant if compared before')
        if graph_c.get_vertex(cl1) != None and graph_c.get_vertex(v1) != None:
            if cl1 in graph_c.get_vertex(v1).get_connections():
                continue
        # if cluster compared to all clusters currently in variant, no need to compare again.
        # implement variant branching later to resolve unlikely exceptions to this, if any.
        skip = 0
        logging.debug('Do not compare cluster to variant if compared to all members of variant individually already')
        if graph_c.get_vertex(cl1) !=None and graph_m.get_vertex(v1) != None:
            skip = 1
            for connection in graph_m.get_vertex(v1).get_connections():
                if connection not in graph_c.get_vertex(cl1).get_connections():
                    skip = 0
                    break
        if skip:
            continue

        # don't compare if exceeds left-sorted comparison bounds.
        # will appear again for right bound comparison later.
        logging.debug('Do not compare to variant if gap exceeds necessary comparison bounds')
        if LR == "L" and (elem.bp1TID != cluster1.lTID or abs(elem.bp1_start - cluster1.l_start) \
            > maxClCombGap) and (elem.bp2TID != cluster1.lTID or abs(elem.bp2_start - cluster1.l_start) \
            > maxClCombGap) and (elem.bp3TID != cluster1.lTID or abs(elem.bp3_start - cluster1.l_start) \
            > maxClCombGap):
            continue
        # analogous to above for right -sorted
        if LR == "R" and (elem.bp1TID != cluster1.rTID or abs(elem.bp1_start - cluster1.r_start) \
            > maxClCombGap) and (elem.bp2TID != cluster1.rTID or abs(elem.bp2_start - cluster1.r_start) \
            > maxClCombGap) and (elem.bp3TID != cluster1.rTID or abs(elem.bp3_start - cluster1.r_start) \
            > maxClCombGap):
            continue
        graph_c.add_edge(cl1,v1,1)

        # check for this cluster's signature and location match with existing variants
        logging.debug('Check conditions for match: cluster against variant if variant is TD_I')
        if elem.SVType == "TD_I":
            # small-medium TD's: second small cluster overlap on other side possible with TD's
            # but not with insertions.
            if cluster1.isSmall and cluster1.lTID == elem.bp2TID and cluster1.rTID == elem.bp2TID and \
                isOverlapping("V", cluster1, elem, "L2", slop) and isOverlapping("V", cluster1, elem, "R2", slop):

                elem.SVType = "TD"
                elem.count+=1
                match = 1
                anyMatch = 1
                elem.complete = 1
            # all insertions
            elif cluster1.lTID == elem.bp1TID and isOverlapping("V", cluster1, elem, "L1", slop):
                if cluster1.l_orient != cluster1.r_orient and \
                    ((not cluster1.l_orient and elem.bp1_hasAlmtR) or \
                    (cluster1.l_orient and elem.bp1_hasAlmtF)):

                    elem.bp1_hasAlmtF = 1
                    elem.bp1_hasAlmtR = 1
                    elem.SVType = "INS"
                    elem.bp3_start, elem.bp3_end = cluster1.r_start, cluster1.r_end
                    elem.bp3TID = cluster1.rTID
                    elem.count+=1
                    match = 1
                    anyMatch = 1
            elif cluster1.rTID == elem.bp1TID and isOverlapping("V", cluster1, elem, "R1", slop):
                if cluster1.l_orient != cluster1.r_orient and \
                    ((not cluster1.r_orient and elem.bp1_hasAlmtR) or \
                    (cluster1.r_orient and elem.bp1_hasAlmtF)):

                    elem.bp1_hasAlmtF = 1
                    elem.bp1_hasAlmtR = 1
                    elem.SVType = "INS"
                    elem.bp3_start, elem.bp3_end = cluster1.l_start, cluster1.l_end
                    elem.bp3TID = cluster1.lTID
                    elem.count+=1
                    match = 1
                    anyMatch = 1
            # if cluster indicating translocation deletion arrives here
            elif cluster1.l_orient == 0 and cluster1.l_orient != cluster1.r_orient and \
                cluster1.rTID == cluster1.lTID == elem.bp2TID and \
                isOverlapping("V", cluster1, elem, "L2", slop) and not isOverlapping("V", cluster1, elem, "R1", slop)\
                and elem.bp2_end < cluster1.r_start < cluster1.r_end < elem.bp1_start:

                    match = 1
                    elem.count+=1
                    anyMatch = 1
                    elem.bp3_start, elem.bp3_end = cluster1.r_start, cluster1.r_end
                    elem.bp3TID = cluster1.rTID
                    elem.SVType = "INS_C"
            elif cluster1.l_orient == 0 and cluster1.l_orient != cluster1.r_orient and \
                cluster1.lTID == cluster1.rTID == elem.bp2TID and \
                not isOverlapping("V", cluster1, elem, "L1", slop) and isOverlapping("V", cluster1, elem, "R2", slop)\
                and elem.bp1_end < cluster1.l_start < cluster1.l_end < elem.bp2_start:

                    match = 1
                    elem.count+=1
                    anyMatch = 1
                    elem.bp3_start, elem.bp3_end = cluster1.l_start, cluster1.l_end
                    elem.bp3TID = cluster1.lTID
                    elem.SVType = "INS_C"
        elif elem.SVType == "INS_C" or elem.SVType == "INS_C_I" or \
            elem.SVType == "INS_C_P" or elem.SVType == "INS_C_I_P":
            logging.debug('Check conditions for match: cluster against variant if variant is in INS_C family')

            if cluster1.lTID == elem.bp1TID and cluster1.rTID == elem.bp3TID and elem.bp1_orient != -1 \
                and elem.bp3_orient != -1 and isOverlapping("V",cluster1,elem,"L1", slop) and \
                cluster1.l_orient != elem.bp1_orient and isOverlapping("V",cluster1,elem,"R3", slop) \
                and cluster1.r_orient != elem.bp3_orient:

                match = 1
                elem.count+=1
                anyMatch = 1
                elem.bp1_orient = -1
                elem.bp3_orient = -1
                # bp has 2 overlapping reads now
                #print "INS_C 1", elem.count
            elif cluster1.lTID == elem.bp3TID and cluster1.rTID == elem.bp1TID and elem.bp1_orient != -1 \
                and elem.bp3_orient != -1 and isOverlapping("V",cluster1,elem,"R1", slop) and \
                cluster1.r_orient != elem.bp1_orient and isOverlapping("V",cluster1,elem,"L3", slop) \
                and cluster1.l_orient != elem.bp3_orient:

                match = 1
                elem.count+=1
                anyMatch = 1
                elem.bp1_orient = -1
                elem.bp3_orient = -1
                #print "INS_C 2", elem.count
            elif cluster1.lTID == elem.bp1TID and cluster1.rTID == elem.bp2TID and elem.bp1_orient != -1 \
                and elem.bp2_orient != -1 and isOverlapping("V",cluster1,elem,"L1", slop) and \
                cluster1.l_orient != elem.bp1_orient and isOverlapping("V",cluster1,elem,"R2", slop) \
                and cluster1.r_orient != elem.bp2_orient:

                match = 1
                elem.count+=1
                anyMatch = 1
                elem.bp1_orient = -1
                elem.bp3_orient = -1
                # bp has 2 overlapping reads now
                #print "INS_C 1", elem.count
            elif cluster1.lTID == elem.bp2TID and cluster1.rTID == elem.bp1TID and elem.bp1_orient != -1 \
                and elem.bp2_orient != -1 and isOverlapping("V",cluster1,elem,"R1", slop) and \
                cluster1.r_orient != elem.bp1_orient and isOverlapping("V",cluster1,elem,"L2", slop) \
                and cluster1.l_orient != elem.bp2_orient:

                match = 1
                elem.count+=1
                anyMatch = 1
                elem.bp1_orient = -1
                elem.bp3_orient = -1
            # the _P subscript denotes that the breakpoints are confirmed.
            # bp1 is indeed the pasted location for this INS_C
            # small cluster overlap to confirm paste location-- overlaps both bp 1 and 2
            elif isOverlapping("V", cluster1, elem, "L1", slop) and cluster1.lTID == elem.bp1TID and \
                isOverlapping("V", cluster1, elem, "R1", slop) and cluster1.rTID == elem.bp1TID and \
                cluster1.l_orient != cluster1.r_orient:

                match = 1
                elem.count+=1
                anyMatch = 1
                if elem.SVType == "INS_C":
                    elem.SVType = "INS_C_P"
                elif elem.SVType =="INS_C_I":
                    elem.SVType = "INS_C_I_P"
                elem.complete = 1
                #print "INS_C_P 3", elem.count
            elif not (elem.SVType == "INS_C_P" or elem.SVType == "INS_C_I_P") and \
                isOverlapping("V", cluster1, elem, "L3", slop) and cluster1.lTID == elem.bp3TID and \
                isOverlapping("V", cluster1, elem, "R3", slop) and cluster1.rTID == elem.bp3TID and \
                cluster1.l_orient != cluster1.r_orient:

                #print "INS_C_P 4", elem.count
                match = 1
                elem.count+=1
                anyMatch = 1
                # make bp1 the paste location since now it is known
                elem.bp1_start, elem.bp1_end, elem.bp3_start, elem.bp3_end=\
                    elem.bp3_start, elem.bp3_end, elem.bp1_start, elem.bp1_end
                #elem.bp1TID, elem.bp3TID = elem.bp3TID, elem.bp1TID #not needed
                elem.bp1_orient, elem.bp3_orient = elem.bp3_orient, elem.bp1_orient
                elem.complete = 1
                if elem.SVType == "INS_C":
                    elem.SVType = "INS_C_P"
                elif elem.SVType =="INS_C_I":
                    elem.SVType = "INS_C_I_P"
        elif elem.SVType == "TD":
            logging.debug('Check conditions for match: cluster against variant if variant is TD')            
            if cluster1.isSmall and (isOverlapping("V", cluster1, elem, "L1", slop) and \
                cluster1.lTID == elem.bp1TID) and (isOverlapping("V", cluster1, elem, "R2", slop) and \
                cluster1.rTID == elem.bp2TID):

                match = 1
                elem.count+=1
                anyMatch = 1
                elem.complete = 1
        elif elem.SVType == "INS" or elem.SVType == "INS_I":
            logging.debug('Check conditions for match: cluster against variant if variant is regular INS')        
            # 2-bp-thus-far INS's
            if elem.bp3_start == -1:
                if cluster1.l_orient == 0 and cluster1.l_orient != cluster1.r_orient and \
                    cluster1.lTID == elem.bp2TID  and isOverlapping("V", cluster1, elem, "L2", slop) and \
                    (cluster1.rTID != elem.bp1TID or not isOverlapping("V", cluster1, elem, "R1", slop))\
                    and (cluster1.lTID != elem.bp1TID or\
                    elem.bp2_end < cluster1.r_start < cluster1.r_end):

                    match = 1
                    elem.count+=1
                    anyMatch = 1
                    elem.bp3_start, elem.bp3_end = cluster1.r_start, cluster1.r_end
                    elem.bp3TID = cluster1.rTID
                    if elem.SVType == "INS":
                        elem.SVType = "INS_C"
                    elif elem.SVType == "INS_I":
                        elem.SVType = "INS_C_I"
                elif cluster1.l_orient == 0 and cluster1.l_orient != cluster1.r_orient and \
                    cluster1.rTID == elem.bp2TID and (cluster1.lTID != elem.bp1TID or \
                    not isOverlapping("V", cluster1, elem, "L1", slop)) and \
                    isOverlapping("V", cluster1, elem, "R2", slop) and (cluster1.lTID != elem.bp1TID\
                    or cluster1.l_start < cluster1.l_end < elem.bp2_start):

                    match = 1
                    elem.count+=1
                    anyMatch = 1
                    elem.bp3_start, elem.bp3_end = cluster1.l_start, cluster1.l_end
                    elem.bp3TID = cluster1.lTID
                    if elem.SVType == "INS":
                        elem.SVType = "INS_C"
                    elif elem.SVType == "INS_I":
                        elem.SVType = "INS_C_I"
                elif (cluster1.lTID == elem.bp1TID and isOverlapping("V", cluster1, elem, "L1", slop)\
                    and elem.bp1_hasAlmtF and cluster1.l_orient == 1)\
                    or (cluster1.rTID == elem.bp1TID and isOverlapping("V", cluster1, elem, "R1", slop)\
                    and elem.bp1_hasAlmtR and cluster1.l_orient == 0):

                    if (elem.SVType == "INS" and cluster1.l_orient != cluster1.r_orient) or\
                        (elem.SVType == "INS_I" and cluster1.l_orient == cluster1.r_orient):
                            match = 1
                            elem.count+=1
                            anyMatch = 1
                            if elem.bp1_hasAlmtF and cluster1.l_orient == 1:
                                elem.bp3_start, elem.bp3_end = cluster1.r_start, cluster1.r_end
                                elem.bp3TID = cluster1.rTID
                            elif elem.bp1_hasAlmtR and cluster1.l_orient == 0:
                                elem.bp3_start, elem.bp3_end = cluster1.l_start, cluster1.l_end
                                elem.bp3TID = cluster1.lTID
            elif elem.bp2TID == cluster1.lTID == cluster1.rTID and isOverlapping("V", cluster1, \
                elem, "L2", slop) and isOverlapping("V", cluster1, elem, "R3", slop) and \
                cluster1.l_orient != elem.bp2_orient and cluster1.r_orient != elem.bp3_orient:

                    match = 1
                    elem.count+=1
                    anyMatch = 1
                    #print "INS -> INS_C 1", elem.count
                    # breakpoints are confirmed for translocation only if INS involves 2 diff chr
                    if elem.bp1TID == elem.bp2TID:
                        if elem.SVType == "INS_I" and cluster1.l_orient == cluster1.r_orient:
                            elem.SVType = "INS_C_I_P"
                            # swap bp1, bp3; then bp2, bp3 for both ustream and dstream
                            elem.bp1_start, elem.bp1_end, elem.bp3_start, elem.bp3_end=\
                                elem.bp3_start, elem.bp3_end, elem.bp1_start, elem.bp1_end
                            elem.bp1TID, elem.bp3TID = elem.bp3TID, elem.bp1TID
                            elem.bp1_orient, elem.bp3_orient = elem.bp3_orient, elem.bp1_orient
                            elem.bp2_start, elem.bp2_end, elem.bp3_start, elem.bp3_end=\
                                elem.bp3_start, elem.bp3_end, elem.bp2_start, elem.bp2_end
                            elem.bp2TID, elem.bp3TID = elem.bp3TID, elem.bp2TID
                            elem.bp2_orient, elem.bp3_orient = elem.bp3_orient, elem.bp2_orient
                        elif cluster1.l_orient == 0 and cluster1.r_orient == 1:
                            if elem.SVType == "INS":
                                elem.SVType = "INS_C"
                            if elem.SVType == "INS_I":
                                elem.SVType = "INS_C_I_P"
                    else:
                        if elem.SVType == "INS_I":
                            elem.SVType = "INS_C_I_P"
                        elif elem.SVType == "INS":
                            elem.SVType = "INS_C_P"
            # small cluster check
            elif cluster1.lTID == elem.bp1TID and cluster1.rTID == elem.bp1TID and \
                isOverlapping("V", cluster1, elem, "L1",slop) and isOverlapping("V", cluster1, elem, "R1", slop):

                match = 1
                elem.count+=1
                anyMatch = 1
        logging.debug('Append cluster to complex variant if supports it')
        if match:
            #print "Cluster", cluster1, "and variant", g, "matched~"
            if cluster1.mapNum not in elem.clusterNums:
                elem.clusterNums.append(cluster1.mapNum)
            graph_m.add_edge(cl1,v1,1)
    logging.debug('Write cluster to claimed list if matched with existing variants')
    if anyMatch:
        claimedCls.add(cluster1.mapNum)

def refreshCCList(consolidatedCls, tidListL, tidListR, maxClCombGap, cluster_L, cluster_R, \
    refRate, consolidatedCls_C):

    if len(consolidatedCls) % refRate == 0:
        for h,item in enumerate(consolidatedCls):
            varTIDs = set([item.bp1TID, item.bp2TID, item.bp3TID])
            # if ref IDs in variant have occurred in both cluster lists
            if varTIDs.issubset(tidListL) and varTIDs.issubset(tidListR):
                if cluster_L.lTID not in varTIDs and cluster_R.rTID not in varTIDs and h < len(consolidatedCls) -1:
                    consolidatedCls_C.append(item)
                    del consolidatedCls[h]
                elif (cluster_L.l_start - item.bp1_start > maxClCombGap or cluster_L.lTID != item.bp1TID) and \
                    (cluster_R.rTID != item.bp1TID or cluster_R.r_start - item.bp1_start > maxClCombGap):
                    if (cluster_L.l_start - item.bp2_start > maxClCombGap or cluster_L.lTID != item.bp2TID) and \
                        (cluster_R.rTID != item.bp2TID or cluster_R.r_start - item.bp2_start > maxClCombGap):
                        if (cluster_L.l_start - item.bp3_start > maxClCombGap or cluster_L.lTID != item.bp3TID) and \
                            (cluster_R.rTID != item.bp3TID or cluster_R.r_start - item.bp3_start > maxClCombGap):
                            if h < len(consolidatedCls) -1:
                                consolidatedCls_C.append(item)
                                del consolidatedCls[h]

def refreshClusterBuffer(clusters, LR, tidListL, tidListR, cluster_L, cluster_R,  maxClCombGap):
    counter=0
    for item in clusters:
        if LR == "L":
            relevTID = item.lTID
            relevStart = item.l_start
        elif LR == "R":
            relevTID = item.rTID
            relevStart = item.r_start
        if relevTID in tidListL and relevTID in tidListR:
            if cluster_L.lTID != relevTID and cluster_R.rTID != relevTID:
                counter+=1
            elif (cluster_L.lTID != relevTID or cluster_L.l_start - relevStart > maxClCombGap) and \
                    (cluster_R.rTID != relevTID or cluster_R.r_start - relevStart > maxClCombGap):
                counter+=1
            else:
                break
    clusters = clusters[counter:]

def consolidatePEClusters(workDir, statFile, clusterFileLS, clusterFileRS,
                          clusterMapFile, maxClCombGap, slop, refRate):
    variantNum = 1
    RDL_Factor=1.2 # default recommended
    fStat = open(statFile,"r")
    RDL = int(fStat.readline().split()[0])
    disc_thresh = int(fStat.readlines()[5].split()[0])
    # clusters sorted by left TID and position and right TID and position for faster comparison
    fClustersLS = open(clusterFileLS,"r")
    fClustersRS = open(clusterFileRS,"r")
    fClusterMap = open(clusterMapFile, "r")
    fVariantsPE = open(workDir+"/allVariants.pe.txt","w")
    fVariantMapPE = open(workDir+"/variantMap.pe.txt", "w")
    clusters_LS = []
    clusters_RS = []
    consolidatedCls = []
    consolidatedCls_C = []
    # maintain list of reference IDs that have been read in clusters
    tidListL = set([-1]) #"-1" represents unset 3rd breakpoint in variant
    tidListR = set([-1])
    claimedCls = set() # clusters that have matched with other clusters or variants
    newBlockL = 0
    newBlockR = 0
    matchGraph = Graph() # stores all cluster and variant matches thus far
    # compareGraph graph stores all comparisons between variants and clusters
    # to prevent repeat comparisons
    compareGraph = Graph()
    loopCount = 0

    # form all possible complex variants from PE clusters
    logging.debug('Started comparison of clusters')
    for lineLC in fClustersLS:
        # Both files are same size so will reach end at same time
        lineRC = fClustersRS.readline()
        cluster_L = clusterI(lineLC)
        cluster_R = clusterI(lineRC)
        if cluster_L.lTID not in tidListL:
            tidListL.add(cluster_L.lTID)
        if cluster_R.rTID not in tidListR:
            tidListR.add(cluster_R.rTID)

        # refresh and match clusters to variants in complex variant buffer
        if len(consolidatedCls) > 0:
            # refresh complex variant buffer list for necessary comparisons only
            refreshCCList(consolidatedCls, tidListL, tidListR, maxClCombGap,
                          cluster_L, cluster_R, refRate, consolidatedCls_C)
            # may need both L,R cluster comparisons since variant buffer is small
            compareVariant(cluster_L, consolidatedCls, claimedCls, compareGraph,
                           matchGraph, "L", maxClCombGap, slop)
            compareVariant(cluster_R, consolidatedCls, claimedCls, compareGraph,
                           matchGraph, "R", maxClCombGap, slop)

        # refresh left and right cluster buffers and mutually match
        if len(clusters_LS) > 0:
            refreshClusterBuffer(clusters_LS, "L", tidListL, tidListR,
                                 cluster_L, cluster_R, maxClCombGap)
        if len(clusters_RS) > 0:
            refreshClusterBuffer(clusters_RS, "R", tidListL, tidListR,
                                 cluster_L, cluster_R, maxClCombGap)
        compareCluster(cluster_L, clusters_LS, claimedCls, compareGraph,
                       matchGraph, "LL", consolidatedCls, slop, RDL_Factor, RDL)
        compareCluster(cluster_L, clusters_RS, claimedCls, compareGraph,
                       matchGraph, "LR", consolidatedCls, slop, RDL_Factor, RDL)
        clusters_LS.append(cluster_L)
        compareCluster(cluster_R, clusters_LS, claimedCls, compareGraph,
                       matchGraph, "RL", consolidatedCls, slop, RDL_Factor, RDL)
        compareCluster(cluster_R, clusters_RS, claimedCls, compareGraph,
                       matchGraph, "RR", consolidatedCls, slop, RDL_Factor, RDL)
        clusters_RS.append(cluster_R)
    logging.debug('Finished comparison of clusters')

    # write list of complex variants and complex variant map to 2 sep files
    hashedVM = {}
    consolidatedCls_C = consolidatedCls_C + consolidatedCls
    consolidatedCls = []
    for line in fClusterMap:
        line_split = line.split()
        clNum = int(line_split[0])
        clList = line_split[1:]
        hashedVM[clNum] = clList
        #print clNum
    logging.debug('Started writing complex variants to output files')
    writeVariants(consolidatedCls_C, fVariantsPE, fVariantMapPE, hashedVM, 0)
    logging.debug('Finished writing complex variants to output files')

    # write unmatched clusters as appropriate variants
    TDStore = []
    TDArtefacts = []
    varCount = len(consolidatedCls_C)
    consolidatedCls_C = []
    fClaimed = open(workDir+"/claimedClusters.txt","w")
    fClustersLS.seek(0)
    logging.debug('Started recording unclaimed clusters')
    for line in fClustersLS:
        clusterC = clusterI(line)
        if not clusterC.mapNum in claimedCls:
            #print "Unclaimed", clusterC.mapNum
            store =1
            newSimpleSV = consCluster()
            newSimpleSV.bp1TID = clusterC.lTID
            newSimpleSV.bp2TID = clusterC.rTID
            newSimpleSV.bp1_start = clusterC.l_start
            newSimpleSV.bp1_end = clusterC.l_end
            newSimpleSV.bp2_start = clusterC.r_start
            newSimpleSV.bp2_end = clusterC.r_end
            newSimpleSV.count = 1
            newSimpleSV.clusterNums.append(clusterC.mapNum)
            # In case did not match with other half cluster for INS_U (de novo INS)
            if clusterC.r_orient == "2":
                newSimpleSV.SVType = "INS_U"
                newSimpleSV.bp2_start = newSimpleSV.bp1_start
                newSimpleSV.bp2_end = newSimpleSV.bp1_end
            elif clusterC.l_orient == 1 and clusterC.r_orient == 0 and \
                clusterC.l_start < clusterC.r_start and \
                clusterC.lTID == clusterC.rTID:
                newSimpleSV.SVType = "TD"
                # avoid double counting TDs -- see artefacts below
                TDStore.append(newSimpleSV)
            elif clusterC.l_orient == clusterC.r_orient and clusterC.lTID == clusterC.rTID:
                newSimpleSV.SVType = "INV"
            elif clusterC.l_orient == 0 and clusterC.r_orient ==1 and \
                clusterC.l_start < clusterC.r_end and \
                clusterC.lTID == clusterC.rTID and \
                clusterC.isSmall != 1:
                newSimpleSV.SVType = "DEL"
            # artefact FR TD cluster
            elif clusterC.l_start < clusterC.r_end and \
                clusterC.l_orient == 0 and \
                clusterC.r_orient == 1 and clusterC.lTID == clusterC.rTID and \
                clusterC.isSmall == 1:
                store=0
                newSimpleSV.SVType = "TD"
                TDArtefacts.append(newSimpleSV)
            else:
                newSimpleSV.SVType = "Unknown"
            # store if appropriate
            if store:
                consolidatedCls_C.append(newSimpleSV)
        else:
            fClaimed.write("%s\t%s\n" %(clusterC.mapNum,str(clusterC.l_orient)+\
                str(clusterC.r_orient)))
    logging.debug('Finished recording unclaimed clusters')

    # $$$ Denovo insertions should not called unless there is SR support
    for elem in TDArtefacts:
        storeTD = 1
        for TD in TDStore:
            if TD.bp1TID == elem.bp1TID and \
               abs(elem.bp1_start-TD.bp2_start) < disc_thresh and \
               abs(elem.bp2_end-TD.bp1_end) < disc_thresh:
                storeTD = 0
                break
        if storeTD:
            consolidatedCls_C.append(elem)
        else:
            newSimpleSV.SVType = "INS_U"
            newSimpleSV.bp1_end = newSimpleSV.bp2_end # set "both" breakpts same
            newSimpleSV.bp2_start = newSimpleSV.bp1_start
            consolidatedCls_C.append(elem)

    logging.debug('Started writing unclaimed variants')
    writeVariants(consolidatedCls_C, fVariantsPE, fVariantMapPE, hashedVM,
                  varCount)
    logging.debug('Finished writing unclaimed variants')
    fClustersLS.close()
    fClustersRS.close()
    fVariantsPE.close()
    fClusterMap.close()
    fVariantMapPE.close()

if __name__ == "__main__":

    # parse arguments
    PARSER = argparse.ArgumentParser(description='Consolidate clusters, formed by formPEClusters.py, into variants', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    PARSER.add_argument('workDir', help='Work directory')
    PARSER.add_argument('statFile', help='File containing BAM statistics, typically bamStats.txt')
    PARSER.add_argument('clusterFileLS', help='File containing clusters sorted by left chr and position, typically allClusters.ls.txt')
    PARSER.add_argument('clusterFileRS', help='File containing clusters sorted by right chr and position, typically allClusters.rs.txt')
    PARSER.add_argument('clusterMapFile', help='File containing list of clusters and their supporting fragments, typically clusterMap.txt')
    PARSER.add_argument('-d', action='store_true', dest='debug',
        help='print debug information')
    PARSER.add_argument('-r', default=550, dest='maxClCombGap', type=int,
        help='Maximum gap between start position of cluster breakpoints to consider them for matching into 1 variant')
    PARSER.add_argument('-s', default=0,dest='slop', type=int,
        help='Additional slop added to dynamically calculated cluster breakpoint margins, if desired. Default recommended.')
    PARSER.add_argument('-f', default=5, dest='refRate', type=int,
        help='Parameter controlling size of buffer variant list in cluster consolidation, optimized for speed.')
    ARGS = PARSER.parse_args()

    LEVEL = logging.INFO
    if ARGS.debug:
        LEVEL = logging.DEBUG

    logging.basicConfig(level=LEVEL,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    consolidatePEClusters(ARGS.workDir, ARGS.statFile, ARGS.clusterFileLS, ARGS.clusterFileRS, ARGS.clusterMapFile, ARGS.maxClCombGap, ARGS.slop, ARGS.refRate)

    logging.shutdown()
