import sys

if __name__=="__main__":

	counter=0
	f=open(sys.argv[1], "r")
	fw=open(sys.argv[1]+".vcf", "w")
	fw.write("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n" %("##fileformat=VCF4.2", "##source=SVCaller-0.0.1","##INFO=<ID=END, Number=1, Type=Integer, Description=\"end pofloat of SV\">","##INFO=<ID=SVTYPE, Number=1, Type=String, Description=\"SV Type\">", "##INFO=<ID=CM, Number=1, Type=String, Description=\"Specific comment\">","##INFO=<ID=ISINV, Number=1, Type=Flag, Description=\"Whether on inverted or positive intand\">", "##INFO=<ID=CHR2, Number=1, Type=Integer, Description=\"source location chromosome of copy-paste INS\">", "##INFO=<ID=POS2, Number=1, Type=Integer, Description=\"source location start of copy-paste INS\">", "##INFO=<ID=END2, Number=1, Type=Integer, Description=\"source location end of copy-paste INS\">"))

        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

	for line in f:
		counter+=1
		line_split = line.split()

		#$Use SR variant flag to remove IMPRECISE
		if line_split[6] == "BND":

			l1 = min(float(line_split[1]), float(line_split[4]))
			l2 = max(float(line_split[1]), float(line_split[4]))	

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], l1,counter,".", "<BND>",".","BND", "SVTYPE=BND;"+"CHR2=" + line_split[3] + ";END=" + str(int(l2)) + "CM="+line_split[7]+";IMPRECISE"))

		elif line_split[6] == "DEL":

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], int(float(line_split[1])),counter,"N", "<DEL>",".","PASS", "SVTYPE=DEL;END="+line_split[5]+";SVLEN=-"+str(int(float(line_split[5])-float(line_split[1])))+";CIPOS=-1,"+str(int(float(line_split[2])-float(line_split[1])))+";CIEND=-"+str(int(float(line_split[2])-float(line_split[1])))+",1;IMPRECISE;PE=.;SR=.","GT","./."))

		elif (line_split[6] == "TD" or line_split[6] == "TD_INV") and float(line_split[2]) > float(line_split[1]):

			isinv="INV=FALSE"
                        if line_split[6]=="TD_INV":
                                isinv="ISINV"

			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], int(float(line_split[2])),counter,".", "<DUP:TANDEM>",".","PASS", "SVTYPE=DUP;END=" + str(int(float(line_split[4])))+";SVLEN="+str(int(float(line_split[4])-float(line_split[2])))+";CIPOS=-50,50;CIEND=-50,50;IMPRECISE;"+isinv))

		elif line_split[6] == "INV":

                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], int(float(line_split[1])) ,counter,".", "<INV>",".","PASS", "SVTYPE=INV;END="+str(int(float(line_split[4])))+";SVLEN="+str(int(float(line_split[4])-float(line_split[2])))+";CIPOS=-50,50;CIEND=-50,50;IMPRECISE"))

		elif line_split[6] == "INS" or line_split[6]=="INS_I":

			isinv="INV=FALSE"
			if line_split[6]=="INS_I":
				isinv="ISINV"
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], int(float(line_split[1])),counter,".", "<DUP>",".","PASS", "SVTYPE=DUP;;END="+str(int(float(line_split[2]))) +";SVLEN="+str(int(float(line_split[2])-float(line_split[1]))) +";CIPOS=-50,50;CIEND=-50,50;IMPRECISE;"+isinv))

		#$for INS_C write separately with comment of bp may be switched, provide alternate	
		elif line_split[6] == "INS_C_P" or line_split[6] == "INS_C_I_P":
		
			isinv="INV=FALSE"
                        if line_split[6] == "INS_C_I_P":
                                isinv="ISINV"
			len_tr = int(float(line_split[2]) - float(line_split[1]))
                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], int(float(line_split[1])),counter,".", "<DEL>",".","PASS", "SVTYPE=DEL;TRAID="+str(counter)+";SVLEN=-" + str(len_tr) +"END="+str(int(float(line_split[2])))+ ";CIPOS=-50,50;CIEND=-50,50;IMPRECISE;"+ isinv))
                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line_split[0], int(float(line_split[1])),counter+1,".", "<DUP>",".","PASS", "SVTYPE=DUP;TRAID="+str(counter) + ";SVLEN=" + str(len_tr) + ";END=" + str(int(float(line_split[2]))) + ";CIPOS=-50,50;CIEND=-50,50;IMPRECISE;"+isinv))
			counter+=1


