#!/usr/bin/env python

# Write a BEDPE file of variants

import argparse
import sys
import logging

def writeBEDs(variantFile, passFile, outname, libINV): 
    passed = None
    if passFile != None:
        passed = set()
        with open(passFile, 'r') as f:
            for line in f:
                passed.add(line.strip())
    if passed != None:
        logging.info('%d variants passed threshold', len(passed))

    outfile = sys.stdout
    if outname != sys.stdout: outfile = open(outname, 'w')
    outfile.write("%s\n" %("chr1\tstart1\tstop1\tchr2\tstart2\tstop2\tnameID\tscore\tstrand1\tstrand2\tSVType\tSupportedBy\tconfPos\tconfEnd\tInverted\tGT\tNSupport\tBNDTag\tSVSubtype\tN_PE\tN_SR\tN_PE_Cl\tBNDAlt1\tBNDAlt2\tGROUPID\tcovLocRejection"))
    with open(variantFile, 'r') as fAV:
        header = fAV.readline()
        for counter, lineAV in enumerate(fAV):
            tokens = lineAV.split()

            svindx = tokens[0]
            if passed is not None and svindx not in passed:
                continue

            svtype = tokens[1]
            chrom1 = tokens[2]
            start1 = int(tokens[3])
            end1 = int(tokens[4])
            chrom2 = tokens[5]
            start2 = int(tokens[6])
            end2 = int(tokens[7])
            chrom3 = tokens[8]
            start3 = int(tokens[9])
            end3 = int(tokens[10])
            support_tag = tokens[11]
            cl_support = tokens[12]
            suppPE = tokens[13]
            suppSR = tokens[14]
            BNDAlt1, BNDAlt2 = ".","."
            GT = "."
            swap = "0"
            bnd = "0"
            support = "0"
            GROUPID = "."
            covInfo = "-1"
            if len(tokens) > 15:
                swap = tokens[15]
                bnd = tokens[16]
                support = tokens[17]
                GT = tokens[18]
                covInfo = tokens[19]

            if svtype == "DN_INS":
                output = tokens[2:8]
                name = "DN_INS"
            if svtype == "DEL":
                output = tokens[2:8]
                name = 'DEL'
            elif svtype == 'TD' or \
                 (svtype == 'TD_I' and support_tag.find("PE") != -1):
                output = tokens[2:8]
                name = 'TD'
            elif svtype == 'INS_I' and chrom1 == chrom2 and \
                 chrom3 == "-1" and support_tag.find("PE") != -1:
                bp1_s, bp1_e = min(start1, start2), min(end1, end2)
                bp2_s, bp2_e = max(start1, start2), max(end1, end2)
                output = [chrom1, bp1_s, bp1_e, chrom2, bp2_s, bp2_e]
                name = 'TD_INV'
            elif svtype == "INV_POSS" or (svtype == "INV" and \
                 ((libINV and support_tag.find("SR") != -1) or  cl_support == "2")):
                output = tokens[2:8]
                name = 'INV'
            elif svtype in ["BND", "Unknown", "INV", "INS_POSS", "TD_I", "INV_POSS", "INS_C",
                            "INS_C_I", "INS_halfFR", "INS_halfRF"] or \
                (svtype in ["INS", "INS_I", "INS_C_P", "INS_C_I_P"] and \
                (start3 == -1 or bnd == "1")):

                name = 'BND_2'
                out3 = "."
                if (svtype in ["INS", "INS_I"] and bnd == "1"):
                    if not svtype.startswith("INS_I"):
                        BNDAlt1, BNDAlt2 = "N[p[", "]p]N"
                        BNDAlt3, BNDAlt4 = "]p]N", "N[p["
                    else:
                        BNDAlt1, BNDAlt2 = "[p[N", "N]p]"
                        BNDAlt3, BNDAlt4 = "N]p]", "[p[N"

                    out1 = [chrom1, start1, end1, chrom2, start2, end2, 
                        "SV" + str(counter), ".", ".", ".", name, support_tag, ".", ".", 
                        ".", GT, support, bnd, svtype, suppPE, suppSR, cl_support, BNDAlt1, BNDAlt2, GROUPID, covInfo]
                    out2 = [chrom1, start1, end1, chrom3, start3, end3, 
                        "SV" + str(counter), ".", ".", ".", name, support_tag, ".", ".", 
                        ".", GT, support, bnd, svtype, suppPE, suppSR, cl_support, BNDAlt3, BNDAlt4, GROUPID, covInfo]
                elif svtype in ["INS_C", "INS_C_I"] or \
                    (svtype in ["INS_C_P", "INS_C_I_P"] and bnd == "1"):
                    GROUPID = "T" + str(counter)
                    if not svtype.startswith("INS_C_I"):
                        BNDAlt1, BNDAlt2 = "N[p[", "]p]N"
                    else:
                        BNDAlt1, BNDAlt2 = "[p[N", "N]p]"
                    BNDAlt3, BNDAlt4 = "N[p[", "]p]N" 

                    out1 = [chrom1, start1, end1, chrom2, start2, end2, 
                        "SV" + str(counter), ".", ".", ".", name, support_tag, ".", ".", 
                        ".", GT, support, bnd, svtype, suppPE, suppSR, cl_support, BNDAlt1, BNDAlt2, GROUPID, covInfo]
                    out2 = [chrom2, start2, end2, chrom3, start3, end3, 
                        "SV" + str(counter), ".", ".", ".", name, support_tag, ".", ".", 
                        ".", GT, support, bnd, svtype, suppPE, suppSR, cl_support, BNDAlt3, BNDAlt4, GROUPID, covInfo]
                    if int(cl_support) >= 3:
                        if not svtype.startswith("INS_C_I"):
                            BNDAlt5, BNDAlt6 = "]p]N", "N[p["
                        else:
                            BNDAlt5, BNDAlt6 = "N]p]", "[p[N"
                        out3 = [chrom1, start1, end1, chrom3, start3, end3,
                            "SV" + str(counter), ".", ".", ".", name, support_tag, ".", ".",
                            ".", GT, support, bnd, svtype, suppPE, suppSR, cl_support, BNDAlt5, BNDAlt6, GROUPID, covInfo]
                else:
                    output = tokens[2:8]
                    name = 'BND'
                    if svtype == "INS_halfFR":
                        BNDAlt1, BNDAlt2 = "N[p[", "]p]N"
                    elif svtype == "INS_halfRF":
                        BNDAlt1, BNDAlt2 = "]p]N", "N[p["
                    elif svtype == "INV":
                        BNDAlt1, BNDAlt2 = "N]p]", "[p[N"

            elif svtype.startswith("INS"):
                bp1_s, bp1_e = start1, end1
                bp2_s, bp2_e = start2, end2
                bp3_s, bp3_e = start3, end3

                if swap == "1" and bp1_s < bp3_s:
                    bp1_s, bp3_s = bp3_s, bp1_s
                    bp1_e, bp3_e = bp3_e, bp1_e
                elif swap == "1":
                    bp1_s, bp2_s = bp2_s, bp1_s
                    bp1_e, bp2_e = bp2_e, bp1_e
                if bp3_s < bp2_s:
                    bp2_s = bp3_s
                    bp3_e = bp2_e

                output = [chrom2, bp2_s, bp3_e, chrom1, bp1_s, bp1_e]
                name = svtype
            else:
                output = tokens[2:8]
                name = 'BND'

            if name != 'BND_2':
                output.extend(["SV" + str(counter), ".", ".", "."])
                output.append(name)
                output.append(support_tag)
                output.extend([".", ".", "."])
                output.extend([GT, support, bnd])

            if name == 'BND': 
                output.extend([svtype, suppPE, suppSR, cl_support, BNDAlt1, BNDAlt2, GROUPID, covInfo])
            elif name != 'BND_2': 
                output.extend(['.', suppPE, suppSR, cl_support, BNDAlt1, BNDAlt2, GROUPID, covInfo])
            
            if name != 'BND_2':
                print >> outfile, "\t".join(map(str, output))
            else:
                out1[10], out2[10] = 'BND', 'BND'
                print >> outfile, "\t".join(map(str, out1))
                print >> outfile, "\t".join(map(str, out2))
                if out3 != ".": 
                    out3[10] = 'BND'
                    print >> outfile, "\t".join(map(str, out3))

    if outname != sys.stdout: outfile.close()

if __name__ == "__main__":
    # parse arguments
    PARSER = argparse.ArgumentParser(description='Write variants in BEDPE format given a allVariant.*.txt file, and optionally a map of variants and fragments/reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    PARSER.add_argument('-d', action='store_true', dest='debug',
                        help='print debug information')
    PARSER.add_argument('-l', action='store_true', dest='libINV',
                        help='Allow liberal inversions')
    PARSER.add_argument('-p', dest='passed', default=None,
                        help='Only print the variants with index in file')
    PARSER.add_argument('-o', dest='out', default=sys.stdout,
                        help='Print the output to this file')
    PARSER.add_argument('variantFile', help='File with the variants')
    ARGS = PARSER.parse_args()

    LEVEL = logging.INFO
    if ARGS.debug:
        LEVEL = logging.DEBUG

    logging.basicConfig(level=LEVEL,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    writeBEDs(ARGS.variantFile, ARGS.passed, ARGS.out, ARGS.libINV)
