import pysam as ps

MIN_PROCESS_SIZE = 3
MIN_VAR_SIZE = 100


def markBadRegionsMQ(wdir, bamfile):  # TODO: unused function
    # read the bam informations
    stats = []
    with open("%s/bamStats.txt" % wdir, "r") as fp:
        for line in fp:
            stats.append(line)
        stats = map(float, stats)
    # print stats

    # read allDiscordants.txt 
    fragments = {}
    with open("%s/allDiscordants.txt" % wdir, "r") as f:
        for line in f:
            tokens = line.strip().split("\t")
            fragments[tokens[0]] = (tokens[1], int(tokens[2]), tokens[3], int(tokens[4]), tokens[5])

    # read the cluster->fragment information
    clusters = {}
    with open("%s/clusterMap.txt" % wdir, "r") as f:
        for line in f:
            tokens = line.strip().split("\t")
            clusters[tokens[0]] = tokens[1:]

    # now read the clusters and process them
    samfile = ps.AlignmentFile(bamfile, "rb")
    fClean = open("%s/allClusters.postClean.txt" % wdir, "w")

    with open("%s/allClusters.txt" % wdir, "r") as f:
        for line in f:
            tokens = line.strip().split("\t")

            # I require a minimum of 3 reads to support the variant
            if int(tokens[1]) <= MIN_PROCESS_SIZE: continue

            # if both the breakpoints are on the same chromosome; it is either 
            # a deletion or a duplication; then the variant 
            # should be at least 100 bps long    
            if tokens[3] == tokens[6] and tokens[2] in ["01", "10"] and \
                    (int(tokens[8]) - int(tokens[4])) < MIN_VAR_SIZE: continue

            index = tokens[0]
            cmap = clusters[index]
            left = []
            right = []
            for findex in cmap:
                frag = fragments[findex]
                left.append(frag[1])
                right.append(frag[3])
            left.sort()
            right.sort()

            # l_is_reverse = fragments[cmap[0]][4][0] == "1"
            # r_is_reverse = fragments[cmap[0]][4][1] == "0"

            # print tokens
            # print min(left), max(left)
            # print min(right), max(right)

            supporting = contradicting = 0
            smq = []
            cmq = []
            reads = set()
            for pc in samfile.pileup(tokens[3], min(left), max(left), stepper="all", truncate=True):
                for pr in pc.pileups:
                    aln = pr.alignment
                    if aln.is_proper_pair: continue
                    if aln.mate_is_unmapped: continue
                    if aln.mapping_quality == 0: continue

                    if aln.query_name not in reads:
                        reads.add(aln.query_name)
                        s = aln.next_reference_start
                        if aln.mate_is_reverse == False: s += stats[0]
                        if aln.next_reference_name != tokens[6] or s < min(right) or s > max(right):
                            # print "1",aln.query_name, aln.flag, aln.next_reference_name, s
                            contradicting += 1
                            cmq.append(aln.mapping_quality)
                        else:
                            supporting += 1
                            smq.append(aln.mapping_quality)

            reads = set()
            for pc in samfile.pileup(tokens[6], min(right), max(right), stepper="all", truncate=True):
                for pr in pc.pileups:
                    aln = pr.alignment
                    if aln.is_proper_pair: continue
                    if aln.mate_is_unmapped: continue
                    if aln.mapping_quality == 0: continue

                    if aln.query_name not in reads:
                        reads.add(aln.query_name)
                        s = aln.next_reference_start
                        if aln.mate_is_reverse == False: s += stats[0]
                        if aln.next_reference_name != tokens[3] or s < min(left) or s > max(left):
                            # print "2",aln.query_name, aln.flag, aln.next_reference_name, s
                            contradicting += 1
                            cmq.append(aln.mapping_quality)
                        else:
                            supporting += 1
                            smq.append(aln.mapping_quality)

            tokens.append(str(supporting))
            try:
                savgmq = sum(smq) * 1. / supporting
            except ZeroDivisionError:
                savgmq = 0
            tokens.append(str(savgmq))
            tokens.append(str(contradicting))
            try:
                cavgmq = sum(cmq) * 1. / contradicting
            except ZeroDivisionError:
                cavgmq = 0
            tokens.append(str(cavgmq))
            fClean.write("%s\n" % ("\t".join(tokens[:10])))

    fClean.close()
