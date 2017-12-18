# coding: utf-8

import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv")
    parser.add_argument("--sample_name")
    parser.add_argument("--clone")
    args = parser.parse_args()

    clone_assign_df = pd.read_csv(args.input_csv)

    bcodes = clone_assign_df.loc[clone_assign_df.clone == args.clone].barcode

    barcodes = list(bcodes)

    paths = ["data/{}/rna/bam/{}_cell_{}.bam".format(args.sample_name, args.sample_name, barcode) for barcode in barcodes]

    print("\n".join(paths))
