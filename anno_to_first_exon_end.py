import sys

with open(sys.argv[1], "r") as input:
    for line in input.readlines():
        line = line.strip().split("\t")
        if line[0] == "protein_coding,protein_coding":
            gene = line[-3]
            enst = line[1]
            chr = line[2]
            strand = line[3]
            starts = line[9].strip(",").split(",")
            ends = line[10] # 1-based
            ends = ends.strip(",").split(",")
            if strand == "+":
                end = int(ends[0])
                print("{}\t{}\t{}\t{}\t{}\t{}".format(chr, end-1, end, enst, gene, strand))
            elif strand == "-":
                end = int(starts[-1]) + 1
                print("{}\t{}\t{}\t{}\t{}\t{}".format(chr, end-1, end, enst, gene, strand))

            