# coding: utf-8
import pandas as pd
import numpy as np
from scipy.stats import binom
import argparse




def test_row(row):
    af1 = row['af1']
    depth1 = row['depth1']
    depth2 = row['depth2']
    p = af1
    if p == 1.0:
        p = 0.99
    pval1 = binom.cdf(0, depth1, p)
    pval2 = binom.cdf(0, depth2, p)
    return pval1 < 0.05 and pval2 < 0.05

def parse_info(AF_line):
    split = AF_line.split(";")
    af1_text = [x for x in split if "AF1" in x][0]
    af1 = float(af1_text.split("=")[1])
    return af1

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_bed")
    parser.add_argument("--depth_template")
    parser.add_argument("--clone")
    parser.add_argument("--output_bed")
    args = parser.parse_args()

    readdepth_clones = list(set(["A", "B", "C"]) - set(args.clone))
    depth_1_file = args.depth_template.replace("CLONE", readdepth_clones[0])
    depth_2_file = args.depth_template.replace("CLONE", readdepth_clones[1])

    df_1 = pd.read_csv(depth_1_file, sep = "\t", header = None, names = ["chrom", "chromEnd", "depth" + "1"])
    df_2 = pd.read_csv(depth_2_file, sep = "\t", header = None, names = ["chrom", "chromEnd", "depth" + "2"])

    # Parse input
    names = ['chrom', 'chromStart', 'chromEnd', 'strand', 'score', 'ref', 'alt', 'alt_strand', 'allele_info', 'info1']
    input_bed_df = pd.read_csv(args.input_bed, sep = "\t", header = None, names = names, index_col = False)
    input_bed_df['af1'] = input_bed_df.allele_info.apply(parse_info)

    df_merged = pd.merge(input_bed_df, df_1, how = 'left', on = ['chrom', 'chromEnd'])
    df_merged = pd.merge(df_merged, df_2, how = 'left', on = ['chrom', 'chromEnd'])

    df_merged.depth1.fillna(0, inplace = True)
    df_merged.depth2.fillna(0, inplace = True)

    keep_variant = df_merged.apply(test_row, axis = 1)

    df_merged = df_merged.loc[keep_variant]

    df_merged.to_csv(args.output_bed, sep = "\t", columns = names, index = False, header = False)


    


    




