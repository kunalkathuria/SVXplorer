#!/usr/bin/env python

# Form PE clusters from discordant PE alignments written by writeDiscordantFragments.py, via a maximal clique formulation
import pysam
import networkx as nx
import argparse as ap
import logging
from sklearn.cluster import KMeans

#global variables
IL_BinDistHash = {}
#margin to increase max_cluster_length due to split reads affecting insert size
SR_GRACE_MARGIN = 40

class fragment(object):
    def __init__(self):
        self.lpos = None
        self.rpos = None
        self.mq = None
    def __str__(self):
        return "%s\t%s\t%s\t" %(self.lpos, self.rpos, self.mq)

class cluster(object):
    # cluster of aligned fragments
    def __init__(self):
        self.l_bound=None
        self.r_bound=None
        self.count=1
        self.cType=None #orientation, e.g. "01" represents "FR" etc.
        self.lTID = -1
        self.rTID = -1
        self.l_min = -1
        self.lmax = -1
        self.r_min = -1
        self.r_max = -1
        # final cluster locations
        self.l_start = -1
        self.l_end = -1
        self.r_start = -1
        self.r_end = -1
        self.clSmall = -1
        self.fragNum = -1
        self.used = 0

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.count, self.cType, self.lTID, self.l_start, self.l_end, self.rTID, self.r_start, self.r_end, self.clSmall)

def readBamStats(statFile):
    stats = []
    with open(statFile,"r") as fp:
        for line in fp:
            stats.append(line)
    stats = map(float,stats)
    return stats[0], stats[1], stats[5], stats[6], stats[7]

def readDistHash(fp):
    global IL_BinDistHash
    total_entries = 0
    for line in fp:
        ls1 = float(line.split()[1])
        IL_BinDistHash[int(line.split()[0])]=ls1
        total_entries+= ls1
    return total_entries

def calcDistPen(f1_lPos, f2_lPos, rdl, dist_end, mean_IL, dist_penalty):
        
    if dist_penalty <= dist_end:
        if dist_end == mean_IL:
            if abs(f2_lPos-f1_lPos)+2*rdl-dist_end <= 0:
                distPen = 1
            else:
                distPen = 0
        else:
            distPen = 1 - (abs(f2_lPos-f1_lPos)+2*rdl- SR_GRACE_MARGIN -dist_end)*1.0/(dist_end - mean_IL)
    else:
        distPen = 1 - (abs(f2_lPos-f1_lPos)+2*rdl- SR_GRACE_MARGIN -dist_end)*1.0/(dist_penalty - dist_end)
    if distPen > 1:
        distPen = 1
    return distPen

def calcEdgeWeight(f1_lPos, f1_rPos, f2_lPos, f2_rPos, IL_BinTotalEntries, c_type, dist_penalty, dist_end, rdl, IL_BinDistHash, mean_IL):
    # first handle one-mapped-read cases
    if f1_rPos == -1 and f2_rPos == -1:
        # calculate weight component 1
        weight = calcDistPen(f1_lPos, f2_lPos, rdl, dist_end, mean_IL, dist_penalty)
        return weight
    elif f1_rPos == -1 and f2_rPos != -1:
        return 0
    elif f1_rPos != -1 and f2_rPos == -1:
        return 0

    # calculate difference in insert lengths between the 2 almts
    bin_size = 10
    if c_type == "01" or c_type == "10":
        ildist = abs(bin_size*int((f2_lPos - f1_lPos + f1_rPos - f2_rPos)/(1.0*bin_size)))
    elif c_type == "00" or c_type == "11":
        ildist = abs(bin_size*int((f2_lPos - f1_lPos + f2_rPos - f1_rPos)/(1.0*bin_size)))

    #larger IL distance for liberal approach
    if ildist > SR_GRACE_MARGIN:
        ildist_L = ildist - abs(bin_size*int(SR_GRACE_MARGIN/(1.0*bin_size)))
    else:
        ildist_L = ildist
    # calculate weight component 1
    distPen = calcDistPen(f1_lPos, f2_lPos, rdl, dist_end, mean_IL, dist_penalty)

    # multiply to weight component 2 to get final weight
    # RF clusters from TDs need not have overlapping almts
    if ildist in IL_BinDistHash and distPen > 0 and (c_type == "10" or (f1_lPos < f2_rPos and f2_lPos < f1_rPos)):
        weight = distPen*IL_BinDistHash[abs(ildist)]/(1.0*IL_BinTotalEntries)
    elif ildist_L in IL_BinDistHash and distPen > 0 and (c_type == "10" or (f1_lPos < f2_rPos and f2_lPos < f1_rPos)):
        weight = distPen*IL_BinDistHash[abs(ildist_L)]/(1.0*IL_BinTotalEntries)
    else:
        weight = 0
    return weight

def calculateMargin(clusterC, max_cluster_length, disc_thresh, bp_margin):

    l_orient = int(clusterC.cType[0])
    r_orient = int(clusterC.cType[1])
    cl_margin_l = max_cluster_length - (clusterC.lmax - clusterC.l_min)
    cl_margin_r = max_cluster_length - (clusterC.r_max - clusterC.r_min)
    if r_orient != 2:
        cl_margin = min(cl_margin_l, cl_margin_r)
    else:
        cl_margin = cl_margin_l

    cl_margin = int(cl_margin)
    if cl_margin < 2*bp_margin:
        # use small margin
        cl_margin = 2*bp_margin
    # small margin on other side of "alignment tip"
    reverse_margin = bp_margin

    # apply margins
    if l_orient == 0:
        clusterC.l_end = clusterC.l_bound + cl_margin
        clusterC.l_start = clusterC.l_bound - reverse_margin
    elif l_orient == 1:
        clusterC.l_start = clusterC.l_bound - cl_margin
        clusterC.l_end = clusterC.l_bound + reverse_margin
    if r_orient == 0:
        clusterC.r_end = clusterC.r_bound + cl_margin
        clusterC.r_start = clusterC.r_bound - reverse_margin
    elif r_orient == 1:
        clusterC.r_start = clusterC.r_bound - cl_margin
        clusterC.r_end = clusterC.r_bound + reverse_margin
    if clusterC.l_start < 0:
        clusterC.l_start = 0
    #if diff chr, "right" almt pos can also be 0
    if r_orient != 2 and clusterC.r_start < 0:
        clusterC.r_start = 0
        clusterC.r_end = max(clusterC.r_start + 1, clusterC.r_end)

    # resolve breakpoint margin overlap if needed
    if clusterC.lTID == clusterC.rTID:
        if not l_orient and r_orient:
            clusterC.l_end = min(clusterC.l_end, clusterC.r_end-1)
            clusterC.r_start = max(clusterC.r_start, clusterC.l_start+1)
        elif not l_orient and not r_orient:
            clusterC.r_start = max(clusterC.lmax+1, clusterC.r_start)
        elif l_orient and r_orient:
            clusterC.l_end = min(clusterC.r_min-1, clusterC.l_end)

def writeClusters(fragGraph, fragHash, fCliques, fClusters, fClusterMap,
                  preserveList, clusterNum, mean_IL, disc_thresh,
                  max_cluster_length, bp_margin, min_cluster_size, debug):

    max_clique_list = []
    # do not form cliques out of most recent fragments, as they will connect further
    for connected_comp in list(nx.connected_component_subgraphs(fragGraph)):
        logging.debug('Connected component: %s', sorted([x for x in connected_comp.nodes()]))
        preserveComp = 0
        if len(preserveList) > 0:
            for node in connected_comp:
                if node in preserveList:
                    preserveComp = 1
                    break
        if len(connected_comp) >= min_cluster_size and not preserveComp:
            cliquesCC = nx.find_cliques(connected_comp)
            max_clique_list.extend(list(cliquesCC))
            logging.debug("Removed %d nodes from graph after forming cliques", len(connected_comp))
            fragGraph.remove_nodes_from(connected_comp)
        elif len(connected_comp) < min_cluster_size and not preserveComp:
            fragGraph.remove_nodes_from(connected_comp)

    max_clique_list_s = sorted(max_clique_list, key=len, reverse=True)
    # delete cliques with too few member fragments
    for k, clique in enumerate(max_clique_list_s):
        if len(clique) < min_cluster_size:
            break
    if len(max_clique_list_s) > 0 and k != len(max_clique_list_s) - 1:
        del max_clique_list_s[k:]

    # process and write clusters from cliques
    while len(max_clique_list_s) > 0:
        clique = max(max_clique_list_s, key=len)
        if len(clique) >= min_cluster_size:
            cliqueRev = clique
            sample = []
            for item in clique:
                if sample == []:
                    sample = [[fragHash[item].l_bound, fragHash[item].r_bound]]
                else:
                    sample.append([fragHash[item].l_bound, fragHash[item].r_bound])

            if len(sample) > 1:
                kmeans = KMeans(n_clusters=2)
                kmeans = kmeans.fit(sample)
                label = kmeans.predict(sample)
                labelList = list(label)
                count0 = labelList.count(0)
                count1 = len(labelList) - count0
                if 1.0*count1/count0 < (1/6.) or 1.0*count0/count1 < (1/6.):
                    if count0 < count1:
                        removeBit = 0
                    else:
                        removeBit = 1
                    removeList = set()
                    for k,item in enumerate(clique):
                        if labelList[k] == removeBit:
                            removeList.add(k)
                    cliqueRev = [v for i, v in enumerate(clique) if i not in removeList]

            pickedFrags = {}
            goodFrags = []
            count = 0
            first = 1

            for item in cliqueRev:
                # do not put diff almts of same *fragment* in same cluster
                item_s = item.split("_")[0]
                if item_s not in pickedFrags and not fragHash[item].used:
                    if first == 1: 
                        lTID = fragHash[item].lTID
                        rTID = fragHash[item].rTID
                        lbp = fragHash[item].l_bound
                        rbp = fragHash[item].r_bound
                        l_min = fragHash[item].l_bound
                        lmax = fragHash[item].l_bound + 1
                        r_min = fragHash[item].r_bound
                        r_max = fragHash[item].r_bound + 1
                        first = 0

                    pickedFrags[item_s] = 1
                    count+=1
                    goodFrags.append(item)
                    fragHash[item].used = 1
                    if fragHash[item].cType[0] == "0" and fragHash[item].l_bound > lbp:
                        lbp = fragHash[item].l_bound
                    elif fragHash[item].cType[0] == "1" and fragHash[item].l_bound < lbp:
                        lbp = fragHash[item].l_bound
                    if fragHash[item].cType[1] == "0" and fragHash[item].r_bound > rbp:
                        rbp = fragHash[item].r_bound
                    elif fragHash[item].cType[1] == "1" and fragHash[item].r_bound < rbp:
                        rbp = fragHash[item].r_bound
                    if fragHash[item].l_bound < l_min:
                        l_min = fragHash[item].l_bound
                    elif fragHash[item].l_bound > lmax:
                        lmax = fragHash[item].l_bound
                    if fragHash[item].r_bound < r_min:
                        r_min = fragHash[item].r_bound
                    elif fragHash[item].r_bound > r_max:
                        r_max = fragHash[item].r_bound
            # if enough good support
            if count >= min_cluster_size:
                newCl = cluster()
                newCl.count = count
                newCl.lTID = lTID
                newCl.rTID = rTID
                newCl.l_bound = lbp
                newCl.r_bound = rbp
                newCl.cType = fragHash[goodFrags[0]].cType
                newCl.clSmall = fragHash[goodFrags[0]].clSmall
                newCl.l_min = l_min
                newCl.r_min = r_min
                newCl.lmax = lmax
                newCl.r_max = r_max
                calculateMargin(newCl, max_cluster_length, disc_thresh, bp_margin)
                if debug:
                    fCliques.write("%s\t%s\n" %("@Cluster"+str(clusterNum),newCl))
                fClusters.write("%s\t%s\n" %(clusterNum, newCl))
                fClusterMap.write("%s" %clusterNum)
                for item in goodFrags:
                    # write all supporting fragments to clique info file as well
                    if debug:
                        fCliques.write("%s\t%s\t%s\n" %(item, fragHash[item].l_bound, fragHash[item].r_bound))
                    item_s = item.split("_")[0]
                    fClusterMap.write("\t%s" %item_s)
                fClusterMap.write("\n")
                clusterNum+=1

                # remove each picked element from universe of fragments and resort list
                usedSet = set(goodFrags)
                max_clique_list_T = []
                max_clique_list_s.remove(clique)
                for clNum, cliqueSet in enumerate(max_clique_list_s):
                    cliqueSet_R = [x for x in cliqueSet if x not in usedSet]
                    max_clique_list_T.append(cliqueSet_R)
                max_clique_list_s = max_clique_list_T
            else:
                for item in cliqueRev:
                    if fragHash[item].used:
                        fragHash[item].used = 0
                max_clique_list_s.remove(clique)
        else:
            max_clique_list_s.remove(clique)

    return clusterNum

def doSubsample(fp):
    subsample = 0
    # write special "whether to subsample" routine here if desired
    return subsample

def runSubsample(almt, block_hash):
    found = 0
    process = 0
    for z in range(almt.l_bound - block_gap_thresh, almt.l_bound + block_gap_thresh):
        for q in range(almt.r_bound - block_gap_thresh, almt.r_bound + block_gap_thresh):
            if (z,q) in block_hash and block_hash[(z,q)] > block_thresh:
                found = 1
                block_hash[(z,q)]+=1
                break
            elif (z,q) in block_hash:
                found = 1
                process = 1
                block_hash[(z,q)]+=1
                break
        if found:
            break
    return found, process

def processNewFrag(fragList, almt, IL_BinTotalEntries, fragmentGraph,
                   edge_weight_thresh, dist_penalty, dist_end, rdl,
                   IL_BinDistHash, mean_IL, debug):
    """Include the new fragment into the graph.
    """

    logging.debug("Processing fragment: %s against %d fragments", almt.fragNum, len(fragList))
    for storedAlmt in fragList:
        if storedAlmt.lTID == almt.lTID and \
           storedAlmt.rTID == almt.rTID and \
           storedAlmt.cType == almt.cType and \
           storedAlmt.clSmall == almt.clSmall:
            logging.debug("Comparing %s~%s-%s to %s~%s-%s", almt.fragNum, almt.l_bound, almt.r_bound, storedAlmt.fragNum, storedAlmt.l_bound, storedAlmt.r_bound)
            #compare 2 frags
            f1_lPos = storedAlmt.l_bound
            f1_rPos = storedAlmt.r_bound
            f2_lPos = almt.l_bound
            f2_rPos = almt.r_bound
            edge_weight = calcEdgeWeight(f1_lPos, f1_rPos, f2_lPos, f2_rPos, IL_BinTotalEntries, almt.cType, dist_penalty, dist_end, rdl, IL_BinDistHash, mean_IL)
            if edge_weight > 0: logging.debug("Edge weight is: %f", edge_weight)

            # only add node if calculation threshold to get edge_weight_thresh
            # has been reached, else determine later
            if edge_weight > edge_weight_thresh:
                if storedAlmt.fragNum != almt.fragNum:
                    fragmentGraph.add_edge(storedAlmt.fragNum, almt.fragNum)

def refreshFragList(fragList, almt, fragmentGraph, fragHash,
                    fCliques, fClusters, fClusterMap, max_cluster_length,
                    clusterNum, newClusterBlock, mean_IL, disc_thresh,
                    bp_margin, debug, min_cluster_size):

    if len(fragList) > 1:
        if almt.lTID != fragList[-2].lTID or almt.rTID != fragList[-2].rTID:
            newClusterBlock = 0
            fragList = [fragList[-1]]
            logging.debug('Refreshed list for new chromosomes')
            clusterNum = writeClusters(fragmentGraph, fragHash, fCliques, fClusters, fClusterMap, [], clusterNum, mean_IL, disc_thresh, max_cluster_length, bp_margin, min_cluster_size, debug)
            fragmentGraph.clear()
            fragmentGraph.add_node(fragList[0].fragNum)
        elif almt.lTID == fragList[newClusterBlock].lTID and \
                (almt.l_bound - fragList[newClusterBlock].l_bound) > 2*max_cluster_length:
            nodeList = {}
            for frag in fragList[newClusterBlock:]:
                nodeList[frag.fragNum] = 1
            clusterNum = writeClusters(fragmentGraph, fragHash, fCliques, fClusters, fClusterMap, nodeList, clusterNum, mean_IL, disc_thresh, max_cluster_length, bp_margin, min_cluster_size, debug)
            del fragList[0:newClusterBlock]
            logging.debug('Deleted %d elems from list', newClusterBlock)
            newClusterBlock = len(fragList)-1

    return fragList, clusterNum, newClusterBlock, fragHash

def formPEClusters(workDir, statFile, IL_BinFile, min_cluster_size,
                   disc_enhancer, bp_margin, subsample, debug):
    # read the stats
    rdl, mean_IL, disc_thresh, dist_penalty, dist_end = readBamStats(statFile)
    max_cluster_length = mean_IL + disc_enhancer*disc_thresh - 2*rdl + SR_GRACE_MARGIN
    logging.info('max_cluster_length is %f', max_cluster_length)

    # subsampling routine variables -- should be rarely required if at all
    block_gap_thresh = 25
    block_thresh = 3*min_cluster_size
    block_hash = {}

    edge_weight_thresh = 0.00
    global IL_BinDistHash
    fragList = []
    fBin=open(IL_BinFile,"r")
    IL_BinTotalEntries = readDistHash(fBin)
    fragHash = {}
    secCounter = {}
    fragmentGraph = nx.Graph()
    clusterNum = 1
    newClusterBlock = 0

    fDiscAlmts=open(workDir+"/allDiscordants.txt","r")
    fClusters=open(workDir+"/allClusters.txt","w")
    fClusterMap=open(workDir+"/clusterMap.txt","w")
    fCliques = None
    if debug:
        fCliques = open(workDir+"/clusterCliques.txt","w")

    logging.info('Started PE cluster formation')
    for line_num,line in enumerate(fDiscAlmts):
        if line_num % 1000 == 0:
            logging.info("Processed %s alignments", line_num)

        parsed = line.strip().split()
        almt = cluster()
        almt.l_bound = int(parsed[2])
        almt.r_bound = int(parsed[4])
        almt.cType = parsed[5]
        almt.lTID = parsed[1]
        almt.rTID = parsed[3]
        almt.clSmall = int(parsed[6])
        almt.fragNum = parsed[0]

        #ignore artefact seen often
        if almt.lTID == almt.rTID and almt.l_bound == almt.r_bound and (almt.cType == "00" or almt.cType == "11"):
            logging.debug("Fragment %s appears to be an alignment artefact", almt.fragNum)
            continue

        #criss crossing read from small TD may appear like this 
        if almt.l_bound > almt.r_bound and almt.lTID == almt.rTID and almt.cType == "01":
            logging.debug("Fragment %s appears to be from a TD with criss-crossing alignments", almt.fragNum)
            almt.l_bound, almt.r_bound = almt.r_bound, almt.l_bound
            almt.cType = "10"
            almt.clSmall = -1

        #subsample if alignments close, for speed -- should be rarely needed
        if subsample:
            logging.info('Subsampling since the coverage is extremely high')
            process = 0
            found = 0
            found, process = runSubsample(almt, block_hash)
            if found and not process:
                continue
            elif not found:
                block_hash[(almt.l_bound, almt.r_bound)] = 1

        fragmentGraph.add_node(almt.fragNum)
        if almt.fragNum not in fragHash:
            fragHash[almt.fragNum] = almt
            secCounter[almt.fragNum] = 1
        # add sec almts of same fragment with underscore
        else:
            secCounter[almt.fragNum]+=1
            almt.fragNum=str(almt.fragNum) + "_" + str(secCounter[almt.fragNum])
            fragHash[almt.fragNum] = almt
        fragList.append(almt)

        # process new PE alignment
        processNewFrag(fragList, almt, IL_BinTotalEntries, fragmentGraph, edge_weight_thresh, dist_penalty, dist_end, rdl, IL_BinDistHash,mean_IL, debug)

        # refresh fragment list
        fragList, clusterNum, newClusterBlock, fragHash = refreshFragList(fragList, almt, fragmentGraph, fragHash, fCliques,  fClusters, fClusterMap, max_cluster_length, clusterNum, newClusterBlock, mean_IL, disc_thresh, bp_margin, debug, min_cluster_size)
    logging.info('Finished PE cluster formation')

    logging.info('Started writing remaining clusters')
    clusterNum = writeClusters(fragmentGraph, fragHash, fCliques, fClusters, fClusterMap, [], clusterNum, mean_IL, disc_thresh, max_cluster_length, bp_margin, min_cluster_size, debug)
    logging.info('Finished writing remaining clusters')

    if debug:
        fCliques.close()
    fClusters.close()
    fClusterMap.close()

if __name__== "__main__":
    # parse arguments
    PARSER = ap.ArgumentParser(description="""
            Form PE clusters from discordant PE alignments written by
            writeDiscordantFragments.py, via a maximal clique formulation""",
            formatter_class=ap.ArgumentDefaultsHelpFormatter)
    PARSER.add_argument('workDir', help='Work directory')
    PARSER.add_argument('statFile',
        help='File containing BAM statistics typically named bamStats.txt')
    PARSER.add_argument('IL_BinFile',
        help='File containing insert length statistics typically named binDist.txt')
    PARSER.add_argument('-d', action='store_true', dest='debug', help='print debug information')
    PARSER.add_argument('-m', default=3, dest='min_cluster_size', type=int,
        help='Minimum number of almts per cluster')
    PARSER.add_argument('-x', default=1.67, dest='disc_enhancer', type=float,
        help='Factor to multiply "end" of IL distribution by a "safety" factor so as \
        to not miss larger sampling ILs in certain distributions')
    PARSER.add_argument('-p', default=25, dest='bp_margin', type=int,
        help='Breakpoint margin for each breakpoint if its calculation yields zero \
        or negative value')
    PARSER.add_argument('-s', action='store_true',
        help='Subsample alignments. Set to 1 ONLY if code is slow, e.g. taking hours on large file')
    ARGS = PARSER.parse_args()

    LEVEL = logging.INFO
    if ARGS.debug:
        LEVEL = logging.DEBUG

    logging.basicConfig(level=LEVEL,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    formPEClusters(ARGS.workDir, ARGS.statFile, ARGS.IL_BinFile,
                   ARGS.min_cluster_size, ARGS.disc_enhancer,
                   ARGS.bp_margin, ARGS.s, ARGS.debug)

    logging.shutdown()
