# coding: utf-8
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import binom
import vcf
import gzip
import argparse

from joblib import Parallel, delayed



def get_depth(chr, snv_pos, file):
    arg = "{}:{}-{}".format(chr, snv_pos, snv_pos)

    result = subprocess.check_output(['samtools', 'depth', '-r', arg, file])
    result = result.decode("utf-8")

    depth = -1
    if len(result) == 0:
        depth = '0'
    else:
        depth = result.split("\t")[2].strip()

    return depth

def wrap_depth_function(series, clone):
    path = 'path' + clone
    return get_depth(series['chr'], series['pos'], series[path])

def filter_record(bam_1, bam_2, record):
    chr = record[0]
    pos = record[1] 
    af1 = record[2]
    af1 = round(af1 / 0.5) * 0.5 # Round to nearest .5
    depth1 = float(get_depth(chr, pos, bam_1))
    depth2 = float(get_depth(chr, pos, bam_2))

    p = af1
    if p == 1.0:
        p = 0.99

    pval1 = binom.cdf(0, depth1, p)
    pval2 = binom.cdf(0, depth2, p)

    if pval1 < 0.05 and pval2 < 0.05:
        return "{}_{}".format(chr, pos)
    else:
        return None



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf")
    parser.add_argument("--bam_template")
    parser.add_argument("--clone")
    parser.add_argument("--output_vcf")
    args = parser.parse_args()

    readdepth_clones = list(set(["A", "B", "C"]) - set(args.clone))
    bam_1 = args.bam_template.replace("CLONE", readdepth_clones[0])
    bam_2 = args.bam_template.replace("CLONE", readdepth_clones[1])

    vcf_reader = vcf.Reader(open(args.input_vcf, 'rb'))
    vcf_writer = vcf.Writer(open(args.output_vcf, 'w'), vcf_reader)

    records = [[record.CHROM, record.POS, float(record.INFO['AF1'])] for record in vcf_reader]
    keys = ["{}_{}".format(record.CHROM, record.POS) for record in vcf_reader]

    record_dict = dict(zip(keys, [record for record in vcf_reader]))


    passed_record_keys = Parallel(n_jobs = 20)(delayed(filter_record)(bam_1, bam_2, record) for record in records)

    passed_record_keys = [x for x in passed_record_keys if x is not None]

    for key in passed_record_keys:
        vcf_writer.write_record(record_dict[key])
    


    




