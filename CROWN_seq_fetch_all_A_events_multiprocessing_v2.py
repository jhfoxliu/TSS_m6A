from Bio import SeqIO
import argparse
import pysam
from Bio.Seq import reverse_complement
import sys, os
import multiprocessing
from multiprocessing import Process,Pool
import signal
import time
from time import gmtime, strftime

def call_all_bases_from_a_read(chr, BAM, main_pid):
    output_temp = {}  # structure, (chr, pos, strand, base) for first base -> (chr, pos, strand) for A -> [A, T, C, G]

    with pysam.AlignmentFile(BAM, "rb") as input_BAM:
        chr_names = {}
        for SQ in input_BAM.header["SQ"]:
            chr_names[SQ["SN"]] = 1

        # chr_v2 = chr
        # if "M" in chr_names:
        #     if chr_v2 == "MT":
        #         chr_v2 = "M"
        # elif "MT" in chr_names:
        #     if chr_v2 == "M":
        #         chr_v2 = "MT"

        for read in input_BAM.fetch(chr, multiple_iterators=True):
            if read.is_secondary == False and read.is_unmapped == False:
                if read.is_read1 and read.is_reverse == False:
                    align_pairs = read.get_aligned_pairs(matches_only=False)               
                    TSN_detected = False
                    TSN_key = None

                    N = -1
                    EJC = 0

                    cigars = [[i[0], i[1]] for i in read.cigartuples][::-1]
                    cigar = cigars.pop()
                    last_cigar_status = cigar[0]
                    for item in align_pairs:
                        current_pos = item[1]

                        if cigars and cigar[1] == 0:
                            last_cigar_status = cigar[0]
                            cigar = cigars.pop()
                            if last_cigar_status == 3 and cigar[0] == 0:
                                EJC += 1
                        last_cigar_status = cigar[0]
                        status = cigar[0]
                        cigar[1] -= 1
                        if status == 0:
                            N += 1
                        if N > read.query_length - options.ignore_read_end:
                            continue
                        
                        if item[1] is not None:
                            current_base = reference_genome[chr][current_pos]
                            key = (read.reference_name, current_pos, "+", current_base)
                            if TSN_detected == False:
                                TSN_detected = True
                                TSN_key = (read.reference_name, current_pos, "+", current_base) # 0-based 
                                if TSN_key not in output_temp:
                                    output_temp[TSN_key] = {}
                            
                            if current_base == "A":
                                if key not in output_temp[TSN_key] and item[0] is not None:
                                    output_temp[TSN_key][key] = {"A":0, "T":0, "C":0, "G": 0, "N": 0, "dist_to_5p": N, "EJC": EJC}
                                if item[0] is not None:
                                    output_temp[TSN_key][key][read.query_sequence[item[0]]] += 1
                            elif key == TSN_key:
                                if key not in output_temp[TSN_key] and item[0] is not None:
                                    output_temp[TSN_key][key] = {"A":0, "T":0, "C":0, "G": 0, "N": 0, "dist_to_5p": N, "EJC": EJC}
                                if item[0] is not None:
                                    output_temp[TSN_key][key][read.query_sequence[item[0]]] += 1
                        else:
                            continue

                elif read.is_read1 and read.is_reverse == True:
                    align_pairs = read.get_aligned_pairs(matches_only=False)[::-1]
                    TSN_detected = False

                    N = -1
                    EJC = 0
                    cigars = [[i[0], i[1]] for i in read.cigartuples] # [::-1]
                    cigar = cigars.pop()
                    last_cigar_status = cigar[0]
                    for item in align_pairs:
                        current_pos = item[1]

                        if cigars and cigar[1] == 0:
                            last_cigar_status = cigar[0]
                            cigar = cigars.pop()
                        if last_cigar_status == 3 and cigar[0] == 0:
                            EJC += 1
                        last_cigar_status = cigar[0]
                        status = cigar[0]
                        cigar[1] -= 1
                        if status == 0:
                            N += 1
                        if N > read.query_length - options.ignore_read_end:
                            continue
                        if item[1] is not None:
                            current_base = reverse_complement(reference_genome[chr][current_pos])
                            key = (read.reference_name, current_pos, "-", current_base)
                            if TSN_detected == False:
                                TSN_detected = True
                                TSN_key = (read.reference_name, current_pos, "-", current_base) # 0-based 
                                if TSN_key not in output_temp:
                                    output_temp[TSN_key] = {}
                            if current_base == "A":
                                if key not in output_temp[TSN_key] and item[0] is not None:
                                    output_temp[TSN_key][key] = {"A":0, "T":0, "C":0, "G": 0, "N": 0, "dist_to_5p": N, "EJC": EJC}
                                if item[0] is not None:
                                    output_temp[TSN_key][key][reverse_complement(read.query_sequence[item[0]])] += 1
                            elif key == TSN_key:
                                if key not in output_temp[TSN_key] and item[0] is not None:
                                    output_temp[TSN_key][key] = {"A":0, "T":0, "C":0, "G": 0, "N": 0, "dist_to_5p": N, "EJC": EJC}
                                if item[0] is not None:
                                    output_temp[TSN_key][key][reverse_complement(read.query_sequence[item[0]])] += 1
                        else:
                            continue

        temp_fn = "{}_{}_temp.all_bases.txt".format(chr, main_pid)
        with open(temp_fn, "w") as output:
            for TSS_key, subdict in output_temp.items():
                for A_key, bases_counts in subdict.items():
                    output.write("{chr}\t{pos}\t{strand}\t{ref_base}\t{ref_cov}\t{target_pos}\t{A}\t{T}\t{C}\t{G}\t{N}\t{dist}\t{EJC}\n".format(
                        chr = TSS_key[0], pos = TSS_key[1] + 1, strand = TSS_key[2], ref_base = TSS_key[3], 
                        ref_cov = output_temp[TSS_key][TSS_key]["A"] + output_temp[TSS_key][TSS_key]["T"] + output_temp[TSS_key][TSS_key]["C"] + output_temp[TSS_key][TSS_key]["G"],
                        target_pos = A_key[1] + 1,
                        A = bases_counts["A"], T = bases_counts["T"], C = bases_counts["C"], G = bases_counts["G"], N = bases_counts["N"], dist = bases_counts["dist_to_5p"], EJC = bases_counts["EJC"]
                    ))
        
        return temp_fn

def signal_handler(sig,frame):
	pool.terminate()
	sys.exit()

if __name__ == "__main__":
    description = """"""
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-r","--ref",dest="reference",required=True, help="Reference fasta")
    group_required.add_argument("-b","--input_bam",dest="input_bam",required=True, help="Input BAM file")
    group_required.add_argument("-o","--output",dest="output",required=False, default="TSS_non_conversion.txt", help="Output txt file name")
    group_required.add_argument("-p",dest="process",default=4, type=int,help="Number of processors, default=4")
    group_required.add_argument("--ignore",dest="ignore_read_end",default=5, type=int,help="Number of reads at read end to be ignored, default=5")
    options = parser.parse_args()

    reference = options.reference
    BAM = options.input_bam

    reference_genome = {}
    for seq in SeqIO.parse(reference,"fasta"):
        reference_genome[seq.id] = str(seq.seq).upper()

    with pysam.AlignmentFile(BAM, "rb") as input_BAM:
        chr_names = {}
        for SQ in input_BAM.header["SQ"]:
            chr_names[SQ["SN"]] = 1

    main_pid = os.getpid()
    all_temp_files = "{}_{}_temp.all_bases.txt".format("*", main_pid)

    sys.stderr.write("[%s] Running...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    signal.signal(signal.SIGINT,signal_handler)
    # call_all_bases_from_a_read("X", BAM, main_pid)
    pool = multiprocessing.Pool(options.process)
    try:
        for chr in chr_names.keys():
            pool.apply_async(call_all_bases_from_a_read,args=(chr, BAM, main_pid,))
        pool.close()
        pool.join()

        # print "Merging tmp files..."
        sys.stderr.write("[%s] Merging TEMPs...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        os.system("cat %s > %s" % (all_temp_files, options.output))
        # print "Removing pileup tmp files..."
        os.system('rm %s' % all_temp_files)
    finally:
        pool.terminate()

    sys.stderr.write("[%s] Genome pileup finished.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))