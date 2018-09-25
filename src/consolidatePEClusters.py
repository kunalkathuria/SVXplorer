#!/usr/bin/env python

# Form variants from clusters output by formPEClusters.py, saved in allClusters.txt, i.e. "cluster consolidation."

import argparse
import logging
from interlap import InterLap
from collections import OrderedDict

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
        if int(data_split[9]) == 1:
            self.isSmall = True
        else:
            self.isSmall = False

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
        self.orient = "22"

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
    # if middle breakpoint is overlapping, then has to be a cut insertion
    lBP_start = min(newVariant.bp2_start, newVariant.bp3_start)
    rBP_end = max(newVariant.bp2_end, newVariant.bp3_end)

    if (LR == "RL" or LR == "LR") and newVariant.bp1TID == newVariant.bp2TID == newVariant.bp3TID:
        if lBP_start <  newVariant.bp1_start < newVariant.bp1_end < rBP_end:
            if newVariant.SVType == "INS":
                newVariant.SVType = "INS_C"
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
        newVariant.bp1TID != newVariant.bp3TID:
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
        newVariant.bp1TID != newVariant.bp2TID:
        if cluster1.lTID == cluster1.rTID and \
            not (cluster1.l_orient == 0 and cluster1.r_orient == 1):
            return 0
        if clusterP.lTID == clusterP.rTID and \
            not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
            return 0
        logging.debug('Orientation matches for INS_C')
        if newVariant.SVType == "INS":
            newVariant.SVType = "INS_C_P"
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

def writeVariants(varHash, fpAV, fpCVM, hashedVM, offset):
    for g, varNum in enumerate(varHash):
        num = offset + 1 + g
        fpCVM.write("%s" %(num))
        item = varHash[varNum]
        if item.SVType == "INS" or item.SVType == "INS_C_P":
                # just to be safe, but should be ensured anyway
                if item.bp2TID == item.bp3TID and item.bp3_start < item.bp2_start:
                    item.bp2_start, item.bp2_end, item.bp3_start, item.bp3_end=\
                        item.bp3_start, item.bp3_end, item.bp2_start, item.bp2_end
        elif item.SVType.find("TD") != -1:
            if item.bp1TID == item.bp2TID and item.bp1_start > item.bp2_start:
                item.bp1_start, item.bp1_end, item.bp2_start, item.bp2_end=\
                    item.bp2_start, item.bp2_end, item.bp1_start, item.bp1_end

        suppCount = 0
        for elem in item.clusterNums:
            try:
                for elem2 in hashedVM[elem]:
                    suppCount+=1
                    fpCVM.write("\t%s" %elem2)
            except:
                print "Exception writing in Variant Map from cluster", elem
                exit(1)
        fpAV.write("%s\t%s\t%s\t%s\t%s\n" %(num,item,suppCount,"0",item.orient))
        fpCVM.write("\n")

def compareCluster(cluster1, clusters, claimedCls, consolidatedCls, 
                    slop, RDL_Factor, RDL, as_relative_thresh, interVariant, dontCompareSet):
    anyMatch = 0
    cl1 = 'c' + str(cluster1.mapNum)
    logging.debug('Start loop to compare this cluster to all clusters in list')
    for clusterP in clusters:
        cl2 = 'c' + str(clusterP.mapNum)
        if cluster1.mapNum == clusterP.mapNum or (dontCompareSet is not None and \
            ((cluster1.mapNum, clusterP.mapNum) in dontCompareSet or (clusterP.mapNum, cluster1.mapNum) in dontCompareSet)):
            logging.debug("Not comparing cluster %s and cluster %s", cl1, cl2)
            continue
        logging.debug("Comparing cluster %s and cluster %s", cl1, cl2)

        # determine which sides of clusters overlap
        logging.debug('Determining signature of 2-cluster overlap: left with left (LL) etc.')
        LLOverlap, LLOverlap2, LROverlap, RLOverlap, RROverlap, LLOverlap3, RLOverlap2, LROverlap2 = 0,0,0,0,0,0,0,0

        if cluster1.rTID != "None" and clusterP.rTID != "None":
            if cluster1.lTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "LL", slop):
                LLOverlap = 1
            if cluster1.lTID == clusterP.rTID and isOverlapping("C", cluster1, clusterP, "LR", slop):
                LROverlap = 1
            if cluster1.rTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "RL", slop):
                RLOverlap = 1
            if cluster1.rTID == clusterP.rTID and isOverlapping("C", cluster1, clusterP, "RR", slop):
                RROverlap = 1
        elif cluster1.rTID == "None" and clusterP.rTID == "None" and \
            cluster1.lTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "LL", slop):
            LLOverlap2 = 1
        elif cluster1.rTID == "None" and clusterP.rTID != "None" and clusterP.isSmall and \
            cluster1.lTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "LL", slop):
            LLOverlap3 = 1
        elif cluster1.rTID == "None" and clusterP.rTID != "None" and clusterP.isSmall and\
            cluster1.lTID == clusterP.rTID and isOverlapping("C", cluster1, clusterP, "LR", slop):
            LROverlap2 = 1
        elif cluster1.rTID != "None" and clusterP.rTID == "None" and cluster1.isSmall and \
            cluster1.lTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "LL", slop):
            LLOverlap3 = 1
        elif cluster1.rTID != "None" and clusterP.rTID == "None" and cluster1.isSmall and \
            cluster1.rTID == clusterP.lTID and isOverlapping("C", cluster1, clusterP, "RL", slop):
            RLOverlap2 = 1

        if not LLOverlap and not LROverlap and not RLOverlap and not RROverlap and not LLOverlap2 and\
            not RLOverlap2 and not LROverlap2 and not LLOverlap3:
            logging.debug('Continue as no overlap between any breakpoints')
            continue

        # initialize all variables
        newSVFlag = 0
        newVariant = consCluster()
        newVariant.SVType = None

        if (LLOverlap2 and cluster1.l_orient != clusterP.l_orient) or LLOverlap3:
            logging.debug('Tagged as De Novo INS')
            newVariant.SVType = "DN_INS"
            newSVFlag=1
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LL")
            newVariant.bp1TID = cluster1.lTID
        elif LROverlap2:
            logging.debug('Tagged as De Novo INS')
            newVariant.SVType = "DN_INS"
            newSVFlag=1
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "LR")
            newVariant.bp1TID = cluster1.lTID
        elif RLOverlap2:
            logging.debug('Tagged as De Novo INS')
            newVariant.SVType = "DN_INS"
            newSVFlag=1
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "RL")
            newVariant.bp1TID = clusterP.lTID
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
            and not clusterP.isSmall and clusterP.l_orient != clusterP.r_orient:
            logging.debug('Tagged as Small INS')
            if clusterP.l_orient != clusterP.r_orient:
                if not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID) \
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
            clusterP.isSmall and not cluster1.isSmall and cluster1.l_orient != cluster1.r_orient:

            if cluster1.l_orient != cluster1.r_orient:
                if not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID)\
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
            cluster1.isSmall and not clusterP.isSmall and clusterP.l_orient != clusterP.r_orient:

            if clusterP.l_orient != clusterP.r_orient:
                if not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID)\
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
            clusterP.isSmall and not cluster1.isSmall and cluster1.l_orient != cluster1.r_orient:

            if cluster1.l_orient != cluster1.r_orient:
                if not (cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID)\
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
            cluster1.rTID or cluster1.lTID == clusterP.rTID or cluster1.rTID == clusterP.rTID) and \
            clusterP.l_orient != clusterP.r_orient and cluster1.l_orient != cluster1.r_orient:
            logging.debug('Large INS check 1: left with left overlap')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.l_end < cluster1.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.l_end < clusterP.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2')
                continue
            if cluster1.l_orient != cluster1.r_orient and clusterP.l_orient != clusterP.r_orient:
                newVariant.SVType = "INS"
                
            if newVariant.SVType == "INS":
                if cluster1.rTID == clusterP.rTID and cluster1.r_end < clusterP.r_start and \
                    not (cluster1.r_orient == 1 and clusterP.r_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:1")
                    continue
                if cluster1.rTID == clusterP.rTID and clusterP.r_end < cluster1.r_start and \
                    not (clusterP.r_orient == 1 and cluster1.r_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:2")
                    continue

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
            (cluster1.rTID == cluster1.lTID or cluster1.rTID == clusterP.lTID or cluster1.lTID == clusterP.lTID) and \
            clusterP.l_orient != clusterP.r_orient and cluster1.l_orient != cluster1.r_orient:
            logging.debug('Large INS check 2: right with right overlap')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.r_start > cluster1.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.r_start > clusterP.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2')
                continue
            if cluster1.l_orient != cluster1.r_orient and clusterP.l_orient != clusterP.r_orient:
                newVariant.SVType = "INS"
                
            if newVariant.SVType == "INS": 
                if cluster1.lTID == clusterP.lTID and cluster1.l_end < clusterP.l_start and \
                    not (cluster1.l_orient == 1 and clusterP.l_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:1")
                    continue
                if cluster1.lTID == clusterP.lTID and clusterP.l_end < cluster1.l_start and \
                    not (clusterP.l_orient == 1 and cluster1.l_orient == 0):
                    logging.debug("Orientation mismatch for INS/INS_C:2")
                    continue

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
            cluster1.l_orient != clusterP.r_orient and \
            clusterP.l_orient != clusterP.r_orient and cluster1.l_orient != cluster1.r_orient:
            logging.debug('Large INS check 3: left mate of cluster 1 overlapping with right of 2')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.r_end < cluster1.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.l_start > clusterP.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2')
                continue

            if cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID and \
                not (cluster1.l_orient == 0 and cluster1.r_orient == 1) and \
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
            # as always if paste location is on diff chr, it is regular INS
            if newVariant.bp1TID != newVariant.bp2TID and newVariant.bp1TID != newVariant.bp3TID \
                and cluster1.l_orient != cluster1.r_orient and clusterP.r_orient != clusterP.l_orient:
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

                    # Put 2 and 3 on same chromosome, as 1 is pasted location in our convention
                    newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, newVariant.bp2_end=\
                        newVariant.bp2_start, newVariant.bp2_end, newVariant.bp1_start, newVariant.bp1_end
                    newVariant.bp1TID, newVariant.bp2TID = newVariant.bp2TID, newVariant.bp1TID
                    newVariant.bp1_orient, newVariant.bp2_orient = newVariant.bp2_orient, newVariant.bp1_orient

            if newVariant.SVType == "INS" and \
                cluster1.rTID == clusterP.lTID and cluster1.r_end < clusterP.l_start and \
                not (cluster1.r_orient == 1 and clusterP.l_orient == 0):
                logging.debug("Orientation mismatch for INS/INS_C:1")
                continue
            elif newVariant.SVType == "INS" and \
                cluster1.rTID == clusterP.lTID and clusterP.l_end < cluster1.r_start and \
                not (clusterP.l_orient == 1 and cluster1.r_orient == 0):
                logging.debug("Orientation mismatch for INS/INS_C:2")
                continue

            # make bp2 < bp3 for C_p and regular INS: this is conventional and consistent
            if (newVariant.SVType == "INS_C_P" or newVariant.SVType == "INS") and \
                newVariant.bp2_start > newVariant.bp3_start:

                newVariant.bp2_start, newVariant.bp2_end, newVariant.bp3_start,\
                    newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
                    newVariant.bp2_start, newVariant.bp2_end
                newVariant.bp3_orient, newVariant.bp2_orient = newVariant.bp2_orient,\
                    newVariant.bp3_orient

        elif RLOverlap and not LROverlap and (cluster1.rTID == cluster1.lTID or \
            clusterP.rTID == cluster1.rTID or cluster1.lTID == clusterP.rTID) and \
            cluster1.r_orient != clusterP.l_orient and \
            clusterP.l_orient != clusterP.r_orient and cluster1.l_orient != cluster1.r_orient:
            logging.debug('Large INS check 4: right mate of cluster 1 overlapping with left mate of cluster 2')

            # if mate between the reads that overlap, then not a bona fide match
            if cluster1.lTID == cluster1.rTID and not (clusterP.l_start > cluster1.l_start + RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:1')
                continue
            if clusterP.lTID == clusterP.rTID and not (cluster1.r_end < clusterP.r_end - RDL_Factor*RDL):
                logging.debug('Not safe distance between overlap point and other 2 mate almts:2 %d %d',cluster1.r_end,clusterP.r_end)
                continue

            if cluster1.lTID == cluster1.rTID == clusterP.lTID == clusterP.rTID and \
                not (cluster1.l_orient == 0 and cluster1.r_orient == 1) and \
                not (clusterP.l_orient == 0 and clusterP.r_orient == 1):
                continue
            newVariant.SVType = "INS_C"
            newVariant.bp1_start, newVariant.bp1_end = setBPs(cluster1, clusterP, "RL")
            newVariant.bp2_start, newVariant.bp2_end = cluster1.l_start, cluster1.l_end
            newVariant.bp3_start, newVariant.bp3_end = clusterP.r_start, clusterP.r_end
            newVariant.bp1TID, newVariant.bp2TID, newVariant.bp3TID \
                = cluster1.rTID, cluster1.lTID, clusterP.rTID
            newVariant.bp2_orient, newVariant.bp3_orient = cluster1.l_orient, clusterP.r_orient
            newSVFlag = 1
            # as always if paste location is on diff chr, it is regular INS
            if newVariant.bp1TID != newVariant.bp2TID and newVariant.bp1TID != newVariant.bp3TID \
                and cluster1.l_orient != cluster1.r_orient and clusterP.r_orient != clusterP.l_orient:
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
                    logging.debug('%d %d %d %d', newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, newVariant.bp2_end)
                    # Put 2 and 3 on same chromosome, as 1 is pasted location in our convention
                    newVariant.bp1_start, newVariant.bp1_end, newVariant.bp2_start, newVariant.bp2_end=\
                        newVariant.bp2_start, newVariant.bp2_end, newVariant.bp1_start, newVariant.bp1_end
                    newVariant.bp1TID, newVariant.bp2TID = newVariant.bp2TID, newVariant.bp1TID
                    newVariant.bp1_orient, newVariant.bp2_orient = newVariant.bp2_orient, newVariant.bp1_orient

            if newVariant.SVType == "INS" and \
                cluster1.lTID == clusterP.rTID and cluster1.l_end < clusterP.r_start and \
                not (cluster1.l_orient == 1 and clusterP.r_orient == 0):
                logging.debug("Orientation mismatch for INS/INS_C:1")
                continue
            elif newVariant.SVType == "INS" and \
                cluster1.lTID == clusterP.rTID and clusterP.r_end < cluster1.l_start and \
                not (clusterP.r_orient == 1 and cluster1.l_orient == 0):
                logging.debug("Orientation mismatch for INS/INS_C:2")
                continue

            # make bp2 < bp3: this is conventional and consistent
            if (newVariant.SVType == "INS_C_P" or newVariant.SVType == "INS") and \
                newVariant.bp2_start > newVariant.bp3_start:
                logging.debug('bp2 > bp3')
                newVariant.bp2_start, newVariant.bp2_end, newVariant.bp3_start,\
                    newVariant.bp3_end = newVariant.bp3_start, newVariant.bp3_end, \
                    newVariant.bp2_start, newVariant.bp2_end
                newVariant.bp3_orient, newVariant.bp2_orient = newVariant.bp2_orient,\
                    newVariant.bp3_orient

        # if clusters matched to form new SV
        if newSVFlag and newVariant.SVType != None:
            newVariant.clusterNums.append(cluster1.mapNum)
            newVariant.clusterNums.append(clusterP.mapNum)
            # should never occur after start, so okay to refresh variant buffer
            if len(consolidatedCls) > 0:
                newVariant.variantNum = 1 + consolidatedCls.keys()[-1]
            else:
                logging.debug("CC list empty, about to write new variant")
                newVariant.variantNum = 1
           
            consolidatedCls[newVariant.variantNum] = newVariant
            claimedCls.add(clusterP.mapNum)
            interVariant.add((newVariant.bp1_start, newVariant.bp1_end, newVariant.variantNum, newVariant.bp1TID))
            interVariant.add((newVariant.bp2_start, newVariant.bp2_end, newVariant.variantNum, newVariant.bp2TID))
            if newVariant.bp3TID != -1:
                interVariant.add((newVariant.bp3_start, newVariant.bp3_end, newVariant.variantNum, newVariant.bp3TID))
            anyMatch = 1
            v1 = 'v' + str(newVariant.variantNum)
            logging.debug("Writing new complex variant %s, type: %s", v1, newVariant.SVType)
            if as_relative_thresh > 1:
                logging.debug('Breaking loop after match for primary almts')
                break
    if anyMatch:
        logging.debug('Adding cluster to claimed list')
        claimedCls.add(cluster1.mapNum)

def compareVariant(cluster1, varList, claimedCls, slop, as_relative_thresh,
                  consolidatedCls, dontCompareSet):
    anyMatch = 0
    cl1 = cluster1.mapNum
    logging.debug('Start loop to compare cluster %s to all complex variants in list', cl1)
    for g,elem in enumerate(varList):
        match = 0
        logging.debug("Trying to compare cluster %s and variant %s", cl1, elem.variantNum)
        # if using primary almts only, then no need to form multiple variants with same clusters
        if as_relative_thresh > 1:
            for clusterDCNum in elem.clusterNums:
                dontCompareSet.add((cl1, clusterDCNum))
        # check for this cluster's signature and location match with existing variants
        if elem.SVType == "DN_INS" and (cluster1.isSmall and cluster1.lTID == elem.bp1TID and \
            (isOverlapping("V", cluster1, elem, "L1", slop) or isOverlapping("V", cluster1, elem, "R1", slop)) or \
            (cluster1.r_orient == 2 and cluster1.lTID == elem.bp1TID and isOverlapping("V", cluster1, elem, "L1", slop))):
            elem.count+=1
            match = 1
            anyMatch = 1
        elif elem.SVType == "TD_I":
            logging.debug('Check conditions for match: cluster against variant for TD_I')
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
        elif elem.SVType == "INS_C" or elem.SVType == "INS_C_P":
            logging.debug('Check conditions for match: cluster against variant for INS_C family')

            if cluster1.lTID == elem.bp1TID and cluster1.rTID == elem.bp3TID and elem.bp1_orient != -1 \
                and elem.bp3_orient != -1 and isOverlapping("V",cluster1,elem,"L1", slop) and \
                cluster1.l_orient != elem.bp1_orient and isOverlapping("V",cluster1,elem,"R3", slop) \
                and cluster1.r_orient != elem.bp3_orient:

                match = 1
                elem.count+=1
                anyMatch = 1
                elem.bp1_orient = -1
                elem.bp3_orient = -1
                elem.bp1_start, elem.bp1_end = sorted([elem.bp1_start, elem.bp1_end, cluster1.l_start, cluster1.l_end])[1:3]
                elem.bp3_start, elem.bp3_end = sorted([elem.bp3_start, elem.bp3_end, cluster1.r_start, cluster1.r_end])[1:3]
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
                elem.bp1_start, elem.bp1_end = sorted([elem.bp1_start, elem.bp1_end, cluster1.r_start, cluster1.r_end])[1:3]
                elem.bp3_start, elem.bp3_end = sorted([elem.bp3_start, elem.bp3_end, cluster1.l_start, cluster1.l_end])[1:3]
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
                elem.bp1_start, elem.bp1_end = sorted([elem.bp1_start, elem.bp1_end, cluster1.l_start, cluster1.l_end])[1:3]
                elem.bp2_start, elem.bp2_end = sorted([elem.bp2_start, elem.bp2_end, cluster1.r_start, cluster1.r_end])[1:3]
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
                elem.bp1_start, elem.bp1_end = sorted([elem.bp1_start, elem.bp1_end, cluster1.r_start, cluster1.r_end])[1:3]
                elem.bp2_start, elem.bp2_end = sorted([elem.bp2_start, elem.bp2_end, cluster1.l_start, cluster1.l_end])[1:3]
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
                elem.complete = 1
                #print "INS_C_P 3", elem.count
            elif (not elem.SVType == "INS_C_P") and \
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
        elif elem.SVType == "TD":
            logging.debug('Check conditions for match: cluster against variant for TD')            
            if cluster1.isSmall and (isOverlapping("V", cluster1, elem, "L1", slop) and \
                cluster1.lTID == elem.bp1TID) and (isOverlapping("V", cluster1, elem, "R2", slop) and \
                cluster1.rTID == elem.bp2TID):

                match = 1
                elem.count+=1
                anyMatch = 1
                elem.complete = 1
        elif elem.SVType == "INS":
            logging.debug('Check conditions for match: cluster against variant for INS')        
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
                    elem.bp1_start, elem.bp1_end = sorted([elem.bp1_start, elem.bp1_end, cluster1.r_start, cluster1.r_end])[1:3]
                    elem.bp2_start, elem.bp2_end = sorted([elem.bp2_start, elem.bp2_end, cluster1.l_start, cluster1.l_end])[1:3]
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
                    elem.bp1_start, elem.bp1_end = sorted([elem.bp1_start, elem.bp1_end, cluster1.l_start, cluster1.l_end])[1:3]
                    elem.bp2_start, elem.bp2_end = sorted([elem.bp2_start, elem.bp2_end, cluster1.r_start, cluster1.r_end])[1:3]
            elif elem.bp2TID == cluster1.lTID == cluster1.rTID and isOverlapping("V", cluster1, \
                elem, "L2", slop) and isOverlapping("V", cluster1, elem, "R3", slop) and \
                cluster1.l_orient != elem.bp2_orient and cluster1.r_orient != elem.bp3_orient:

                    match = 1
                    elem.count+=1
                    anyMatch = 1
                    #print "INS -> INS_C 1", elem.count
                    # breakpoints are confirmed for translocation only if INS involves 2 diff chr
                    if elem.bp1TID == elem.bp2TID:
                        if cluster1.l_orient == 0 and cluster1.r_orient == 1:
                            if elem.SVType == "INS":
                                elem.SVType = "INS_C"
                    else:
                        if elem.SVType == "INS":
                            elem.SVType = "INS_C_P"
                    elem.bp3_start, elem.bp3_end = sorted([elem.bp3_start, elem.bp3_end, cluster1.r_start, cluster1.r_end])[1:3]
                    elem.bp2_start, elem.bp2_end = sorted([elem.bp2_start, elem.bp2_end, cluster1.l_start, cluster1.l_end])[1:3]
            # small cluster check
            elif cluster1.lTID == elem.bp1TID and cluster1.rTID == elem.bp1TID and \
                isOverlapping("V", cluster1, elem, "L1",slop) and isOverlapping("V", cluster1, elem, "R1", slop):

                match = 1
                elem.count+=1
                anyMatch = 1
        if match:
            logging.debug('Append cluster %s to complex variant %s', cluster1.mapNum, elem.variantNum)
            if cluster1.mapNum not in elem.clusterNums:
                elem.clusterNums.append(cluster1.mapNum)
            # if using secondary almts, don't compare to clusters in variant again if match
            # done for primaries in comparison stage itself
            if as_relative_thresh <= 1:    
                for clusterDCNum in elem.clusterNums:
                    dontCompareSet.add((cl1, clusterDCNum))
            logging.debug('Breaking loop after match for primary almts')
            if as_relative_thresh > 1:
                break
    if anyMatch:
        logging.debug('Write cluster %s to claimed list', cl1)
        claimedCls.add(cluster1.mapNum)

def consolidatePEClusters(workDir, statFile, clusterFile,
                          clusterMapFile, slop, refRate, as_relative_thresh):
    RDL_Factor=1.2 # default recommended
    fStat = open(statFile,"r")
    RDL = int(fStat.readline().split()[0])
    disc_thresh = int(fStat.readlines()[5].split()[0])
    # clusters sorted by left TID and position and right TID and position for faster comparison
    fClusters = open(clusterFile,"r")
    fClusterMap = open(clusterMapFile, "r")
    fVariantsPE = open(workDir+"/allVariants.pe.txt","w")
    fVariantMapPE = open(workDir+"/variantMap.pe.txt", "w")
    consolidatedCls = OrderedDict()
    claimedCls = set() # clusters that have matched with other clusters or variants
    interCluster = InterLap()
    interVariant = InterLap()
    clusterHash = {}

    # form all possible complex variants from PE clusters
    logging.debug('Started comparison of clusters')
    for lineC in fClusters:
        lineC_split = lineC.split()
        mapNum = int(lineC_split[0])
        clusterHash[mapNum] = clusterI(lineC)
    fClusters.seek(0)

    for lineC in fClusters:
        cluster = clusterI(lineC)
        dontCompareClSet = None
        if len(consolidatedCls) > 0:
            dontCompareClSet = set()
            variants_M = []
            comparedSet = set()
            for variantM in list(interVariant.find((cluster.l_start, cluster.l_end))):
                varNum = variantM[2]
                if varNum not in comparedSet and variantM[3] == cluster.lTID:
                    variants_M.append(consolidatedCls[varNum])
                    comparedSet.add(varNum)
            #except de novo INS candidates        
            if cluster.r_start != -1:
                for variantM in list(interVariant.find((cluster.r_start, cluster.r_end))):
                    varNum = variantM[2]
                    if varNum not in comparedSet and variantM[3] == cluster.rTID:
                        variants_M.append(consolidatedCls[varNum])
                        comparedSet.add(varNum)
            compareVariant(cluster, variants_M, claimedCls, slop, as_relative_thresh,
                           consolidatedCls, dontCompareClSet)
        
        clusters_M = []
        comparedSet = set()
        for clusterM in list(interCluster.find((cluster.l_start, cluster.l_end))):
            mapNum = clusterM[2]
            if mapNum not in comparedSet and clusterM[3] == cluster.lTID:
                clusters_M.append(clusterHash[mapNum])
                comparedSet.add(mapNum)
        if cluster.r_start != -1:        
            for clusterM in list(interCluster.find((cluster.r_start, cluster.r_end))):
                mapNum = clusterM[2]
                if mapNum not in comparedSet and clusterM[3] == cluster.rTID:
                    clusters_M.append(clusterHash[mapNum])
        interCluster.add((cluster.l_start, cluster.l_end, cluster.mapNum, cluster.lTID))
        if cluster.r_start != -1:
            interCluster.add((cluster.r_start, cluster.r_end, cluster.mapNum, cluster.rTID))

        compareCluster(cluster, clusters_M, claimedCls, consolidatedCls, 
                      slop, RDL_Factor, RDL, as_relative_thresh, interVariant, dontCompareClSet)
    logging.debug('Finished comparison of clusters')

    # write list of complex variants and complex variant map to 2 sep files
    fVariantsPE.write("VariantNum\tType\tchr1\tstart1\tstop1\tchr2\tstart2\tstop2\tchr3\tstart3\tstop3\tSupportBy\tNPEClusterSupp\tNFragPESupp\tNFragSRSupp\tOrientation\n") 
    hashedVM = {}
    for line in fClusterMap:
        line_split = line.split()
        clNum = int(line_split[0])
        clList = line_split[1:]
        hashedVM[clNum] = clList
        #print clNum
    logging.debug('Started writing complex variants to output files')
    writeVariants(consolidatedCls, fVariantsPE, fVariantMapPE, hashedVM, 0)
    logging.debug('Finished writing complex variants to output files')

    # write unmatched clusters as appropriate variants
    varNum = 0
    TDStore = []
    TDArtefacts = []
    varCount = len(consolidatedCls)
    consolidatedCls = {}
    fClaimed = open(workDir+"/claimedClusters.txt","w")
    fClusters.seek(0)
    logging.debug('Started recording unclaimed clusters')
    for line in fClusters:
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
            # In case did not match with other half cluster for DN_INS (de novo INS)
            if clusterC.r_orient == 2:
                newSimpleSV.SVType = "DN_INS"
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
                newSimpleSV.orient = str(clusterC.l_orient) + str(clusterC.r_orient)
            #crossover TD cluster
            elif clusterC.l_orient == 0 and clusterC.r_orient ==1 and \
                clusterC.lTID == clusterC.rTID and \
                clusterC.l_start > clusterC.r_end and \
                clusterC.isSmall == 1:
                newSimpleSV.SVType = "TD"
                TDArtefacts.append(newSimpleSV)
            elif clusterC.l_orient == 0 and clusterC.r_orient ==1 and \
                clusterC.l_start < clusterC.r_end and \
                clusterC.lTID == clusterC.rTID and \
                clusterC.isSmall != 1:
                newSimpleSV.SVType = "DEL"
            elif clusterC.lTID != clusterC.rTID and clusterC.l_orient == 0 and \
                clusterC.r_orient == 1:
                newSimpleSV.SVType = "INS_halfFR"
            elif clusterC.lTID != clusterC.rTID and clusterC.l_orient == 1 and \
                clusterC.r_orient == 0:
                newSimpleSV.SVType = "INS_halfRF"
            elif clusterC.lTID != clusterC.rTID and clusterC.l_orient == clusterC.r_orient:
                newSimpleSV.SVType = "INS_half_I"
            else:
                newSimpleSV.SVType = "Unknown"
            # store if appropriate
            if store:
                consolidatedCls[varNum] = newSimpleSV
                varNum+=1
        else:
            fClaimed.write("%s\t%s\n" %(clusterC.mapNum,str(clusterC.l_orient)+\
                str(clusterC.r_orient)))
    logging.debug('Finished recording unclaimed clusters')

    # Unnec now it seems: $$$ Denovo insertions should not called unless there is SR support
    #check if TD already stored via cluster of different signature
    for elem in TDArtefacts:
        storeTD = 1
        for TD in TDStore:
            if TD.bp1TID == elem.bp1TID and \
               abs(elem.bp1_start-TD.bp2_start) < disc_thresh and \
               abs(elem.bp2_end-TD.bp1_end) < disc_thresh:
                storeTD = 0
                break
        # store as de novo INS if was not due to existing TD    
        if storeTD:
            consolidatedCls[varNum] = elem
            varNum+=1

    logging.debug('Started writing unclaimed variants')
    writeVariants(consolidatedCls, fVariantsPE, fVariantMapPE, hashedVM,
                  varCount)
    logging.debug('Finished writing unclaimed variants')
    fClusters.close()
    fVariantsPE.close()
    fClusterMap.close()
    fVariantMapPE.close()

if __name__ == "__main__":

    # parse arguments
    PARSER = argparse.ArgumentParser(description='Consolidate clusters, formed by formPEClusters.py, into variants', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    PARSER.add_argument('workDir', help='Work directory')
    PARSER.add_argument('statFile', help='File containing BAM statistics, typically bamStats.txt')
    PARSER.add_argument('clusterFile', help='File containing discordant clusters')
    PARSER.add_argument('clusterMapFile', help='File containing list of clusters and their supporting fragments, typically clusterMap.txt')
    PARSER.add_argument('-d', action='store_true', dest='debug',
        help='print debug information')
    PARSER.add_argument('-r', default=550, dest='maxClCombGap', type=int,
        help='Maximum gap between start position of cluster breakpoints to consider them for matching into 1 variant')
    PARSER.add_argument('-s', default=0,dest='slop', type=int,
        help='Additional slop added to dynamically calculated cluster breakpoint margins, if desired. Default recommended.')
    PARSER.add_argument('-f', default=5, dest='refRate', type=int,
        help='Parameter controlling size of buffer variant list in cluster consolidation, optimized for speed.')
    PARSER.add_argument('-a', default=2, dest='as_relative_thresh', type=int,
        help='Threshold of relative alignment score of given alignment to primary alignment of same name\
        --set to less than 1 (.95 recommended) if using secondary almts')
    ARGS = PARSER.parse_args()

    LEVEL = logging.INFO
    if ARGS.debug:
        LEVEL = logging.DEBUG

    logging.basicConfig(level=LEVEL,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    consolidatePEClusters(ARGS.workDir, ARGS.statFile, ARGS.clusterFile, ARGS.clusterMapFile, ARGS.slop, ARGS.refRate, ARGS.as_relative_thresh)

    logging.shutdown()
