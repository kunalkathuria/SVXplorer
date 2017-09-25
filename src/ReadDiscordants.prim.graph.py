#KK 08/31/17
#Accepts list of name-sorted discordant read pairs (2 files), filters out all relevant discordant reads and stores them in file with attributes for subsequent cluster formation
#Picks alignment if both reads mapped to same chromosome

import math
import sys
import pysam
import time
from itertools import izip

CALC_THRESH = 500000
PCT_THRESH = 0
FILE1 = "./results/aln1s.bam" # read 1 discordants
FILE2 = "./results/aln2s.bam"# read 2 discordants
BAMFILE = sys.argv[1] # position-sorted same BAM file.
ignoreR_FILE = sys.argv[2]
ignoreChr = sys.argv[3]
MAP_THRESH = int(sys.argv[4])
AS_THRESH = float(sys.argv[5])
AS_CALC_THRESH = .999
SIG_THRESH=.9985
ignoreTIDList = []
chrHash = {}

def calcMeanSig(file2):
   bamfile = pysam.Samfile(file2,"rb")
   summedIL = 0
   summedQL = 0
   counter = 0
   stdev = 0
   max_IL = 0
   loopCount = 0
   IL_list = []

   while counter < CALC_THRESH and loopCount < 2*CALC_THRESH:

        try:

            m1 = bamfile.next()

        except StopIteration:
            break

	#print m1.get_tag("XS")
	#last condition ensures split read will not be included in analysis because "next" may not be mate then
        if m1.is_proper_pair and not m1.is_secondary and not m1.is_supplementary and m1.get_tag("AS") > AS_CALC_THRESH*m1.infer_query_length() and m1.template_length > 0:

	     IL = m1.template_length
	     IL_list.append(IL)	
	     
	     if IL > max_IL:
		max_IL = IL
             summedIL+= IL
	     summedQL+= m1.infer_query_length()
             counter+= 1
	loopCount+=1
   if counter == 0 or counter ==1:
	print "Division by 0. Please check order of name-sorted and position-sorted files supplied."
	return
   meanIL = summedIL/counter
   meanQL = summedQL/counter
   IL_list = sorted(IL_list)
  
   PENALTY_PERC = .90
   #"generalized" 3 sigma distance if not normal distribution, SIG_THRESH = 99.85
   DIST_END = IL_list[int(SIG_THRESH*len(IL_list)) - 1]
   DIST_PEN = IL_list[int(PENALTY_PERC*len(IL_list)) - 1]
   #print "DISC stats:", SIG_THRESH, DISC_D, len(IL_list), IL_list[-50:]
   DISC_D = DIST_END - meanIL
   bamfile.close()

   BIN_SIZE = 10
   BIN_DIST = BIN_SIZE
   DIST_HASH = {}
   minIL = IL_list[0]
   print minIL, IL_list[-1]
   for item in IL_list:
        if item - minIL < BIN_SIZE:
                if BIN_DIST not in DIST_HASH:
                        DIST_HASH[BIN_DIST]=0
                DIST_HASH[BIN_DIST]+=1
        else:
                BIN_DIST+=BIN_SIZE
                minIL = item

   BINDIST_HASH = {}
   for item in DIST_HASH:
        for item2 in DIST_HASH:
                #print item, item2
                temp = item2 - item
                if temp >= 0:
                        if temp not in BINDIST_HASH:
                                BINDIST_HASH[temp] = 0
                        BINDIST_HASH[temp]+= DIST_HASH[item] + DIST_HASH[item2]

   fh = open("./results/bindist.txt","w")
   for item in BINDIST_HASH:
        fh.write("%s\t%s\n" %(item,BINDIST_HASH[item]))

   #"generalized" 3 sigma distance if not normal distribution, SIG_THRESH = 99.85
   DISC_D = IL_list[int(SIG_THRESH*len(IL_list)) - 1]
   DISC_D = DISC_D - meanIL
   bamfile.close()
   
   bamfile = pysam.Samfile(file2,"rb")
   loopCount=0
   counter=0
   while counter < CALC_THRESH and loopCount < 2*CALC_THRESH:

        try:

            m1 = bamfile.next()

        except StopIteration:
            break

	if m1.is_proper_pair and not m1.is_secondary and not m1.is_supplementary and m1.get_tag("AS") > AS_CALC_THRESH*m1.infer_query_length() and m1.template_length > 0:

	     stdev+= (m1.template_length - meanIL)**2
	     counter+=1
	loopCount+=1

   bamfile.close()
   bamfile = pysam.AlignmentFile(file2,"rb")

   cov = 0
   width = 20
   counter2=0

   for pileupcolumn in bamfile.pileup():
   	cov+=pileupcolumn.n
	counter2+=1
	if counter2 > CALC_THRESH:
		break

   bamfile.close()
   stdev = stdev/counter
   stdev = stdev**(.5)
   cov=cov/counter2

   if DISC_D < 0:
        DISC_D = 3*stdev

   #print meanQL, meanIL, stdev, cov, max_IL, DISC_D
   fp = open("./results/bam_stats.txt","w")
   fp.write("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n" %(meanQL, meanIL, stdev,cov, max_IL, DISC_D, DIST_PEN, DIST_END))
   fp.close()
   return meanQL, meanIL, stdev, DISC_D

[RDL, MEAN_D, SIG_D, DISC_dist] = calcMeanSig(BAMFILE)
ignore_buffer =0* RDL

def formHash():

        prevTID = "*" # invalid value
        fo=open(ignoreR_FILE,"r")
        print ignoreR_FILE

        for line in fo:
                line_s = line.split()
                currentTID = line_s[0]
                if currentTID != prevTID:

                        chrHash[currentTID] = {}

                for x in range(int(line_s[1])-ignore_buffer, int(line_s[2])+ignore_buffer):
                        chrHash[currentTID][x] = 1

                prevTID = currentTID

def ignoreRead(chr_l, loc_l, chr_r, loc_r):

        if chr_l in chrHash and loc_l in chrHash[chr_l]:
                return 1

        if chr_r in chrHash and loc_r in chrHash[chr_r]:
                return 1

        return 0

def recordChr():

	f= open(ignoreChr, "r")
	for line in f:
		ignoreTIDList.append(line.split()[0])

class Cluster(object):

    def __init__(self):

       #Whether l_bound is start or end of left-aligned read (in reference) depends on orientation of read in reference
       self.l_bound=None
       self.r_bound=None
       self.l_bound_O=None # debug variable. first pair in cluster: all positions assessed relative to this one for future additions to cluster
       self.r_bound_O=None
       self.C_type=None # LR orientation: F=0, R=1
       self.Cmap_type = None # 0 if both reads map, 1 if one read maps only
       self.ltid = None
       self.rtid = None
       self.mapqual = -1
       self.clsmall = -1
       self.ascore1 = -1
       self.ascore2 = -1

    def __str__(self):
        return "%s %s %s %s %s %s %s %s %s" % (self.ltid, self.l_bound, self.rtid, self.r_bound, self.C_type, self.ascore1, self.ascore2, self.clsmall, self.mapqual)

def FormDiscordant(al1, al2, DList1):

    		OUTER_D = MEAN_D # non-discordant default value. safe.

		[al1_reference_name, al1_reference_start, al1_reference_end, al1_is_unmapped, al1_mapping_quality, al1_is_reverse, Al_S_1, al1_infer_query_length] = [al1.reference_name, al1.reference_start, al1.reference_end, al1.is_unmapped, al1.mapping_quality, al1.is_reverse, float(al1.get_tag("AS")), al1.infer_query_length()]
		
		[al2_reference_name, al2_reference_start, al2_reference_end, al2_is_unmapped, al2_mapping_quality, al2_is_reverse, Al_S_2, al2_infer_query_length] = [al2.reference_name, al2.reference_start, al2.reference_end, al2.is_unmapped, al2.mapping_quality, al2.is_reverse, float(al2.get_tag("AS")), al2.infer_query_length()]
	
                al1_match = 0
		al2_match = 0
                l_orient = 2 # improper value
                r_orient = 2
	
		
		if al1_reference_start == None and al1_reference_end == None and al2_reference_start == None and al2_reference_end == None:
			return

		if al1_reference_name in ignoreTIDList or al2_reference_name in ignoreTIDList:
                        return
	
		try:
                        if ignoreR_FILE != "none" and ignoreRead(al1_reference_name, al1_reference_start, al2_reference_name, al2_reference_start):

                                return
                except:
                        pass

                if Al_S_1 > AS_THRESH*al1_infer_query_length:
			al1_match = 1

                if Al_S_2 > AS_THRESH*al2_infer_query_length:
                    	al2_match = 1

		#want both alignments on same chr
		if (al1_reference_name != al2_reference_name) or (not al1_match and not al2_match):
			return

                if not al1_is_unmapped and not al2_is_unmapped and al1_reference_start!= None and al1_reference_end != None and al1_reference_start > -1 and al2_reference_start!=None and al2_reference_end !=None and al2_reference_start > -1 and al1_match and al2_match:

		    # if primary alignment below mapping quality threshold, don't use fragment at all
		    if (al1_mapping_quality < MAP_THRESH or al2_mapping_quality < MAP_THRESH):
			return
		
 		    if isinstance(al1_mapping_quality, int) and isinstance(al2_mapping_quality, int):
			    map_qual = min(al1_mapping_quality, al2_mapping_quality)

		    temp = Cluster()
                    temp.Cmap_type = 0
		    temp.mapqual = map_qual

                    # All right if TID different or orientation different. Won't make any difference in those situations.
                    if (al1_reference_start <= al2_reference_start):
			OUTER_D = abs(al1.reference_end - al2.reference_start) + al1_infer_query_length + al2_infer_query_length
                        temp.ltid = al1_reference_name
                        temp.rtid = al2_reference_name
                        l_orient = al1_is_reverse
                        r_orient = al2_is_reverse

                    else:
			OUTER_D = abs(al2.reference_end - al1.reference_start) + al1_infer_query_length + al2_infer_query_length
                        temp.ltid = al2_reference_name
                        temp.rtid = al1_reference_name
                        l_orient = al2_is_reverse
                        r_orient = al1_is_reverse
		
		    temp.C_type = str(int(l_orient)) + str(int(r_orient))

		else:
			return

                # If discordant (due to mapping distance or orientation or mapping TID)
                # Check just to be safe even though discordant flag set in creating BAM files
                if (temp.Cmap_type == 0 and abs(OUTER_D - MEAN_D) > DISC_dist) or (temp.Cmap_type==0 and CLType !="01"):
			#print "Marked disc", al1, al2
			temp.ascore1 = Al_S_1
			temp.ascore2 = Al_S_2

			# "small" means discordant due to IL being smaller than threshold	
			if temp.ltid == temp.rtid and OUTER_D - MEAN_D < -1*DISC_dist:
				temp.clsmall = 1

                        if temp.Cmap_type == 0:
                           
			    # "Left" TID by convention is the one with the lower mapped coordinate if TIDs are different
                            if (al1_reference_start <= al2_reference_start) and (not al1_is_reverse) and (not al2_is_reverse):
                        
                                temp.l_bound = al1_reference_end                
                                temp.r_bound = al2_reference_end
      
                            elif (al1_reference_start <= al2_reference_start) and (not al1_is_reverse) and (al2_is_reverse):

                                temp.l_bound = al1_reference_end                
                                temp.r_bound = al2_reference_start
              			
                            elif (al1_reference_start <= al2_reference_start) and (al1_is_reverse) and (not al2_is_reverse):

                                temp.l_bound = al1_reference_start               
                                temp.r_bound = al2_reference_end
     
                            elif (al1_reference_start <= al2_reference_start) and (al1_is_reverse) and (al2_is_reverse):

                                temp.l_bound = al1_reference_start                
                                temp.r_bound = al2_reference_start
             			
                            elif (al1_reference_start > al2_reference_start) and (not al1_is_reverse) and (not al2_is_reverse):
                                
                                temp.l_bound = al2_reference_end              
                                temp.r_bound = al1_reference_end

                            elif (al1_reference_start > al2_reference_start) and (not al1_is_reverse) and (al2_is_reverse):

                                temp.l_bound = al2_reference_start               
                                temp.r_bound = al1_reference_end
         			
                            elif (al1_reference_start > al2_reference_start) and (al1_is_reverse) and (not al2_is_reverse):

                                temp.l_bound = al2_reference_end                
                                temp.r_bound = al1_reference_start
                                
                            elif (al1_reference_start > al2_reference_start) and (al1_is_reverse) and (al2_is_reverse):

                                temp.l_bound = al2_reference_start                
                                temp.r_bound = al1_reference_start
            				
                            temp.l_bound_O = temp.l_bound
                            temp.r_bound_O = temp.r_bound
		
			    	
                        if temp.l_bound == None and temp.r_bound == None:
                            print "WARNING STORING BAD DATA ", temp


                        if temp.Cmap_type==0:
                            	
				DList1.append(temp)
 
                # If even 1 alignment is concordant, don't use any mappings from read
                # Being safe, as discordant flag is anyway set when creating the input BAM files
                else:

                    return
			       
    
def ReadNextReadAlignments(bamname):

    bamfile = pysam.AlignmentFile(bamname)

    for alignment in bamfile:
	if not alignment.is_secondary:
            yield alignment.qname, alignment                        
    
if __name__ == "__main__":

    f1 = open("./results/All_Discords_P.txt","w")
   
    print "Ignoring:", ignoreChr, ignoreR_FILE
 
    if ignoreChr != "none":
	recordChr()
 	print "Chromosomes", ignoreTIDList, "will be ignored."

    if ignoreR_FILE != "none":
        print "Forming hash table for bad regions..."
        formHash()
	print "Done."

    currentFrag = 1

    # Read discordant alignments and write all possible discordant pairs to file. 

    print "Entering read loop..."

    start_for = time.clock()
    
    for (q1,aln1),(q2,aln2) in izip(ReadNextReadAlignments(FILE1),ReadNextReadAlignments(FILE2)):

	if (currentFrag % 100000) == 0:
		print "Fragment", currentFrag, "analyzed."

        assert q1[:-2] == q2[:-2]

        DList1 = []
        
        start_fd = time.clock()
        FormDiscordant(aln1, aln2, DList1)

        for item in DList1:
             f1.write("%s %s\n" %(currentFrag, item))

        end_fd = time.clock()
    
        currentFrag = currentFrag + 1


    end_for = time.clock()

    print "Done with loop"

    f1.close()
