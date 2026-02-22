import pandas as pd
import sys
import argparse
import numpy as np

if __name__ == "__main__":
    description = """"""
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i","--input_csv",dest="input_csv",help="input CSV file")
    group_required.add_argument("-b","--input_bed",dest="input_bed",help="input BED file")
    group_required.add_argument("-o","--output_file",dest="output_file",help="output file")
    options = parser.parse_args()

    df = pd.read_csv(options.input_csv, index_col=0, header=0, low_memory=False)

    df2 = pd.read_csv(options.input_bed, index_col=None, sep="\t", names=["chr", "end_0", "end_1", "TSS", "x", "strand", "xx", "xxx", "xxxx", "xxxxx", "ENST", "Gene", "dist"], low_memory=False)
    # dist, positive = annot in downstream, negative = upstream
    df2 = df2.drop(["x", "xx", "xxx", "xxxx", "xxxxx"], axis=1)
    # print(df2)
    data = {}

    def find_best(x):
        if x["dist"].min() == 0:
            return 0
        else:
            x_over_zero = x[x["dist"]>0]
            if x_over_zero.shape[0] > 0:
                return x_over_zero["dist"].min()
            else:
                return x["dist"].max()

    df3 = df2.groupby(by=["chr", "end_0", "end_1", "TSS", "strand"])[["dist"]].apply(lambda x: find_best(x)).reset_index()
    df3.columns = ["chr", "end_0", "end_1", "TSS", "strand", "dist"]

    def fix_dist(x):
        if x["dist"] <= 2000 and x["dist"] >= -2000:
            if x["strand"] == "+":
                return x["end_1"], x["dist"] + x["end_1"]
            elif x["strand"] == "-":
                return x["end_1"], x["end_1"] - x["dist"]
        else:
            return x["end_1"], x["end_1"]
    df3["end_1_original"], df3["end_1_fix"]  = zip(*df3.apply(lambda x: fix_dist(x), axis=1))
    df3.to_csv(options.output_file + ".fixed.pairs.csv")
    
    data_original = {}
    data_fixed = {}
    
    for idx, row in df3.iterrows():
        key = (row["chr"], row["TSS"], row["strand"])
        data_original[key] = row["end_1_original"]
        data_fixed[key] = row["end_1_fix"]
    
    df["first_exon_end_original"] = df.apply(lambda x: data_original.get((x["chr"], x["TSS"], x["strand"])), axis=1)
    df["first_exon_end_fixed"] = df.apply(lambda x: data_fixed.get((x["chr"], x["TSS"], x["strand"])), axis=1)
    
    df["m6A_to_exon_end_original"] = np.abs(df["Target"] - df["first_exon_end_original"])
    df["m6A_to_exon_end_fixed"] = np.abs(df["Target"] - df["first_exon_end_fixed"])
    
    df.to_csv(options.output_file)



