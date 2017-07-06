#extract primary alignments in VARSECer format
import sys

if __name__== "__main__":

	f2=open(sys.argv[1],"r")
	f1=open("../results/text/All_Discords_P.txt","r")
	f3=open("../results/text/All_Discords_P_N.txt","w")

	skip = 0
	primary = 0
	prev_frag = "*"
	ignore_frag = "$"
	tid_matched_L = 0
	tid_matched_R = 0

	for counter,line1 in enumerate(f1):

		print counter
		line1_s = line1.split()
		current_frag = line1_s[0]
	
		if ignore_frag == current_frag:

			continue

		if current_frag !=  prev_frag:

			primary = 1
		else:
			primary = 0

		prev_frag = current_frag
		
		skip = 0
		f2.seek(0)

		for line2 in f2:

			line2_s = line2.split()

			# if matched with a previous TID but not this one, then no need to go further				
			if tid_matched_L and tid_matched_R and line1_s[1] != line2_s[0] and line1_s[3] != line2_s[0]:
				break

			if line1_s[1] == line2_s[0] and int(line2_s[1]) <= int(line1_s[2]) <= int(line2_s[2]):

				skip = 1

			if line1_s[3] == line2_s[0] and int(line2_s[1]) <= int(line1_s[4]) <= int(line2_s[2]):
	
				skip = 1

			if line1_s[1] == line2_s[0]:
                                tid_matched_L = 1

                        if line1_s[3] == line2_s[0]:
                                tid_matched_R = 1

			if skip:
			
				break

		if not skip:
			f3.write("%s" %line1)
		elif skip and primary:
			ignore_frag = current_frag
			

