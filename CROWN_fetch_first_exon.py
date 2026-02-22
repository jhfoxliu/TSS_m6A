import os, sys
from Bio import SeqIO
import pandas as pd
import argparse
import pysam
import pandas as pd
from Bio.Seq import reverse_complement
import multiprocessing
from multiprocessing import Process,Pool
import signal
import time
from time import gmtime, strftime

if __name__ == "__main__":
    description = """"""
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i","--input_table",dest="input_table",required=True,help="input table file")
    group_required.add_argument("-b","--input_bam",dest="input_bam",required=True,help="input BAM file")
    group_required.add_argument("-o","--output_file",dest="output_file",required=True,help="output file")
    options = parser.parse_args()

    df = pd.read_csv(options.input_table, index_col=None, names=["chr", "TSS", "strand", "TSS_base", "TSS_cov", "Target", "A", "T", "C", "G", "N", "dist_to_5p", "EJC"], sep="\t", low_memory=False)
    df["cov"] = df["A"] + df["G"]
    df["level"] = df["A"] / (df["A"] + df["G"]) 
    df = df[df["cov"] >= 20]
    df = df[(df["C"] + df["T"])/(df["A"] + df["T"] + df["C"] + df["G"]) < 0.05]
    df.to_csv(options.input_table + ".all.cov20.csv")
    data = {}
    for idx, row in df.iterrows():
        key = (str(row["chr"]), int(row["TSS"]), str(row["strand"]))
        data[key] = int(row["TSS"]) # values = closest exon-exon junction corrdinate 1-based

    with pysam.AlignmentFile(options.input_bam, "rb") as bam_file:
        for key, value in data.items():
            chr, pos_1, strand = key
            pos_0 = pos_1 - 1
            for read in bam_file.fetch(contig=chr, start=pos_0, end=pos_1, multiple_iterators=True):
                if strand == "+" and read.is_read1 == True and read.is_reverse == False:
                    if read.reference_start != pos_0 or read.cigartuples[0][0] == 4:
                        continue
                    else:
                        pos_1_mapped_end = pos_0
                        for cigar in  read.cigartuples:
                            cigar_id = cigar[0]
                            length = cigar[1]
                            if cigar_id == 0:
                                pos_1_mapped_end += length
                            elif cigar_id == 2:
                                pos_1_mapped_end += length
                            elif cigar_id == 3:
                                break
                            else:
                                continue
                    if value < pos_1_mapped_end:
                        data[key] = pos_1_mapped_end

                elif strand == "-" and read.is_read1 == True and read.is_reverse == True:
                    if read.reference_end != pos_1 or read.cigartuples[-1][0] == 4:
                        continue
                    else:
                        pos_1_mapped_end = pos_1 + 1
                        for cigar in  read.cigartuples[::-1]:
                            cigar_id = cigar[0]
                            length = cigar[1]
                            if cigar_id == 0:
                                pos_1_mapped_end -= length
                            elif cigar_id == 2:
                                pos_1_mapped_end -= length
                            elif cigar_id == 3:
                                break
                            else:
                                continue

                    if value > pos_1_mapped_end:
                        data[key] = pos_1_mapped_end
    
    with open(options.output_file, "w") as output_file:
        for key, values in data.items():
            print("{chr}\t{pos_0}\t{pos_1}\t{TSS}\t.\t{strand}".format(chr=key[0], TSS = key[1], strand = key[2], pos_0 = values - 1, pos_1 = values),
                   file=output_file)
