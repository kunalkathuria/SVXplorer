import sys

if __name__=="__main__":

	counter=0
	f=open(sys.argv[1], "r")
	fw=open(str(sys.argv[1]+".vcf"), "w")
	fw.write("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n" %("##fileformat=VCF4.3", "##source=SVCaller-0.0.1","##INFO=<ID=END, Number=1, Type=Integer, Description=\"end point of SV\">","##INFO=<ID=SVTYPE, Number=1, Type=String, Description=\"SV Type\">", "##INFO=<ID=CM, Number=1, Type=String, Description=\"Comment denoting likely kind of variant\">","##INFO=<ID=ISINV, Number=1, Type=Flag, Description=\"Whether on inverted or positive strand\">", "##INFO=<ID=CHR2, Number=1, Type=Integer, Description=\"For BNDs the reference ID of the 'END' breakpoint if different from that of start 'POS'\">", "##INFO=<ID=POS2, Number=1, Type=Integer, Description=\"source location start of copy-paste INS\">", "##INFO=<ID=END2, Number=1, Type=Integer, Description=\"source location end of copy-paste INS\">", "##INFO=<ID=GROUPID, Number=1, Type=Integer, Description=\"GROUPID correlating 2 translocation or non-tandem-duplication events\">"))

        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

	for line in f:
		counter+=1
		line_split = line.split()
		precise=line_split[7].find("SR")
		support="SUPPORT=" + line_split[-1] + ";"
		cipos = str(int(line_split[2])-int(line_split[1]))
		svlen = str(int(line_split[4]) - int(line_split[1]))
		if precise == -1:
			precise = "IMPRECISE"
		else:
			precise="PRECISE"

		#$Use SR variant flag to remove IMPRECISE
		if line_split[6] == "BND" and line_split[7] != "INS_C":

			#l1 = min(int(line_split[1]), int(line_split[4]))
			#l2 = max(int(line_split[1]), int(line_split[4]))	

			chr2=""
			if line_split[0] != line_split[3]:
				chr2="CHR2="+ line_split[3] + ";"

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,"N", "<BND>",".","PASS", "SVTYPE=BND;END=" + line_split[5] + ";" + svlen + ";CIPOS=0," + cipos + ";CIEND=-" + cipos + ",0;CM="+line_split[6]+";" + chr2 + support + precise))

		elif line_split[6] == "BND" and line_split[7] == "INS_C":

			#$write exact imprecision and third breakend based on number of supporting clusters for INS_C	
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,"N", "<BND>",".","PASS", "SVTYPE=BND;END=" + line_split[2]+ ";CIPOS=-100,100;CIEND=-100,100;CM=" + line_split[6]+";" + support + precise))
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[3], line_split[2],counter,"N", "<BND>",".","PASS", "SVTYPE=BND;END=" + line_split[5]+ ";CIPOS=-100,100;CIEND=-100,100;CM=" + line_split[6]+";" + support + precise))
			
		elif line_split[6] == "DEL":

			 fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,"N", "<DEL>",".","PASS", "SVTYPE=DEL;END=" + line_split[5] + ";SVLEN=-" + svlen + ";CIPOS=0," + cipos + ";CIEND=-" + cipos + ",0;CM="+line_split[6]+ ";" + support + precise))

		elif line_split[6] == "TD" or line_split[6] == "TD_INV":

			isinv=""
                        if line_split[6]=="TD_INV":
                                isinv="ISINV"

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[2],counter,"N", "<DUP:TANDEM>",".","PASS", "SVTYPE=DUP;END=" + line_split[4] + ";SVLEN=" + svlen + ";CIPOS=-" + cipos + ",0;CIEND=0," + cipos + ";CM="+line_split[6]+";" + support + precise))

		elif line_split[6] == "INV":

			ciend = int(line_split[5]) - int(line_split[4])
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,"N", "<INV>",".","PASS", "SVTYPE=INV;END=" + line_split[5] + ";CIPOS=-" + str(int(cipos)/2.0) +"," + str(int(cipos)/2.0) + ";CIEND=-" + str(int(ciend)/2.0) +"," + str(int(ciend)/2.0) + ";CM="+line_split[6]+ ";" + support + precise))

		elif line_split[6] == "INS" or line_split[6]=="INS_I":

			svlen=str(int(line_split[2])-int(line_split[1]))
			cipos = 2*(int(line_split[5])-int(line_split[4]))
			isinv=""
			if line_split[6]=="INS_I":
				isinv="ISINV"

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,"N", "<DUP>",".","PASS", "SVTYPE=DUP;END=" + line_split[2] + ";SVLEN=" + svlen + ";CIPOS=0," + str(cipos) + ";CIEND=-" + str(cipos) +",0;GROUPID=" + str(counter) + ";CM="+line_split[6]+";" + isinv + ";" + support + precise))
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[3], line_split[4],counter,"N", "<INS:ME>",".","PASS", "SVTYPE=INS;END=" + str(int(line_split[4]) + 1) + ";SVLEN=" + svlen + ";CIPOS=0," + str(cipos/2.0) + ";CIEND=0," + str(cipos/2.0) + ";GROUPID=" + str(counter) + ";CM="+line_split[6]+ ";" + isinv + ";" + support + precise))

		#$for INS_C write separately with comment of bp may be switched, provide alternate	
		elif line_split[6] == "INS_C_P" or line_split[6] == "INS_C_I_P":
	
			cipos = 2*(int(line_split[5])-int(line_split[4]))	
			isinv=""
                        if line_split[6] == "INS_C_I_P":
                                isinv="ISINV"
			len_tr = int(line_split[2]) - int(line_split[1])
                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,"N", "<DEL:ME>",".","PASS", "SVTYPE=DEL;END=" + line_split[2] + ";CIPOS=0," + str(cipos) + ";CIEND=-" + str(cipos) + ";GROUPID="+str(counter)+";SVLEN=-" + str(len_tr) + ";" + isinv + ";" + support + precise))
                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[3], line_split[4],counter+1,"N", "<INS:ME>",".","PASS", "SVTYPE=INS;END=" + str(int(line_split[4]) + 1) +";CIPOS=0," + str(cipos/2.0) + ";CIEND=-1," + str(cipos/2.0) + ";GROUPID="+str(counter) + ";SVLEN=" + str(len_tr) + ";" + isinv + ";" + support + precise))
			counter+=1


