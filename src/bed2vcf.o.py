import sys

if __name__=="__main__":

	counter=0
	f=open(sys.argv[1], "r")
	fw=open(str(sys.argv[1]+".vcf"), "w")
	fw.write("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n" %("##fileformat=VCF4.2", "##source=SVCaller-0.0.1","##INFO=<ID=END, Number=1, Type=Integer, Description=\"end point of SV\">","##INFO=<ID=SVTYPE, Number=1, Type=String, Description=\"SV Type\">", "##INFO=<ID=CM, Number=1, Type=String, Description=\"Specific comment\">","##INFO=<ID=ISINV, Number=1, Type=Flag, Description=\"Whether on inverted or positive strand\">", "##INFO=<ID=CHR2, Number=1, Type=Integer, Description=\"source location chromosome of copy-paste INS\">", "##INFO=<ID=POS2, Number=1, Type=Integer, Description=\"source location start of copy-paste INS\">", "##INFO=<ID=END2, Number=1, Type=Integer, Description=\"source location end of copy-paste INS\">"))

        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

	for line in f:
		counter+=1
		line_split = line.split()

		#$Use SR variant flag to remove IMPRECISE
		if line_split[6] == "BND":

			l1 = min(int(line_split[1]), int(line_split[4]))
			l2 = max(int(line_split[1]), int(line_split[4]))	

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], l1,counter,".", "<BND>",".","BND", "SVTYPE=BND;"+"CHR2=" + line_split[3] + ";END=" + str(l2) + "CM="+line_split[7]+";IMPRECISE"))

		elif line_split[6] == "DEL":

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,"N", "<DEL>",".","PASS", "SVTYPE=DEL;END="+line_split[5]+";SVLEN=-"+str(int(line_split[5])-int(line_split[1]))+";CIPOS=-1,"+str(int(line_split[2])-int(line_split[1]))+";CIEND=-"+str(int(line_split[2])-int(line_split[1]))+",1;IMPRECISE;PE=.;SR=.","GT","./."))

		elif line_split[6] == "TD" or line_split[6] == "TD_INV":

			isinv="INV=FALSE"
                        if line_split[6]=="TD_INV":
                                isinv="ISINV"

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[2],counter,".", "<DUP:TANDEM>",".","PASS", "SVTYPE=DUP;END=" + line_split[4]+";SVLEN="+str(int(line_split[4])-int(line_split[2]))+";CIPOS=-"+str(int(line_split[2])-int(line_split[1]))+",1;CIEND=-1,"+str(int(line_split[2])-int(line_split[1]))+";IMPRECISE;"+isinv))

		elif line_split[6] == "INV":

                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,".", "<INV>",".","PASS", "SVTYPE=INV;END="+line_split[4]+";SVLEN="+str(int(line_split[4])-int(line_split[2]))+";CIPOS=-1,"+str(int(line_split[2])-int(line_split[1]))+";CIEND=-"+str(int(line_split[2])-int(line_split[1]))+";IMPRECISE"))

		elif line_split[6] == "INS" or line_split[6]=="INS_I":

			isinv="INV=FALSE"
			if line_split[6]=="INS_I":
				isinv="ISINV"
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[3], line_split[4],counter,".", "<DUP:ISP>",".","PASS", "SVTYPE=DUP;CHR2="+line_split[0]+";POS2="+line_split[1]+";END2="+line_split[2] +";SVLEN="+str(int(line_split[2])-int(line_split[1])) +";CIPOS=-1,"+str(int(line_split[5])-int(line_split[4]))+";CIEND=-"+str(int(line_split[5])-int(line_split[4]))+",1;IMPRECISE;"+isinv))

		#$for INS_C write separately with comment of bp may be switched, provide alternate	
		elif line_split[6] == "INS_C" or line_split[6] == "INS_C_P" or line_split[6] == "INS_C_I" or line_split[6] == "INS_C_I_P":
		
			isinv="INV=FALSE"
                        if line_split[6]=="INS_C_I" or line_split[6] == "INS_C_I_P":
                                isinv="ISINV"
			len_tr = int(line_split[2]) - int(line_split[1])
                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], line_split[1],counter,".", "<DEL:TRA>",".","PASS", "SVTYPE=DEL;TRAID="+str(counter)+";SVLEN=-" + str(len_tr) +"END="+line_split[2]+ ";CIPOS=-50,50;CIEND=-50,50;IMPRECISE;"+ isinv))
                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[3], line_split[4],counter+1,".", "<DUP:TRA>",".","PASS", "SVTYPE=DUP;TRAID="+str(counter) + ";SVLEN=" + str(len_tr) + ";END="+str(int(line_split[2])-int(line_split[1])) + ";CHR2=" + line_split[0] + ";POS2=" + line_split[1] + ";END2=" + line_split[2] + ";CIPOS=-1,"+str(int(line_split[5])-int(line_split[4]))+";CIEND=-"+str(int(line_split[5])-int(line_split[4]))+",1;IMPRECISE;"+isinv))
			counter+=1


