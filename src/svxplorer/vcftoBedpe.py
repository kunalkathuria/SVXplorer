import sys


def main():
    # TODO: add parser here, remove sys.argv
    const = int(sys.argv[2])
    search = 0
    f = open(sys.argv[1], "r")
    fw = open(str(sys.argv[1] + ".bedpe"), "w")
    first = 0

    for line in f:
        line_split = line.split()
        if line_split[0][0] != "#":
            search = 1
        if search:
            eflag = 0
            end = -1
            start = -1
            if len(line_split) > 7:
                eflag = line_split[7].find(";END")
                if eflag == -1:
                    eflag = line_split[7].find("END")
                    short_str = line_split[7][eflag:]
                    end = int(short_str[4: 1 + short_str[1:].find(";")])
                else:
                    short_str = line_split[7][eflag:]
                    end = int(short_str[5: 1 + short_str[1:].find(";")])
                start = int(line_split[1])
                # these SVs assumed to be on same chr
                if start > end:
                    start, end = end, start

        if (search and line_split[0][0] != "#"):
            fw.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                line_split[0], start - const / 2, start + const / 2, line_split[0], end - const / 2, end + const / 2))

    fw.close()


if __name__ == "__main__":
    main()
