#!/usr/bin/env python3
# Only keep regions covered in at least a coverage of 10X in gnomAD

import pandas as pd
import argparse
import hail as hl

descr = "Filter_gnomAD_coverage_file.py \n"
descr += "usage: python src/Filter_gnomAD_coverage_file.py -d 100g -e 64g -c 50 -t /data-tmp/antoine_data/hailTMP \n"
parser = argparse.ArgumentParser(description=descr)
parser.add_argument('-d', '--driver_memory', required=True, help='memory for Spark driver (e.g. 100g)')
parser.add_argument('-e', '--executor_memory', required=False, default="64g", help='memory for Spark executors (e.g. 64g)')
parser.add_argument('-c', '--cores', required=True, help='number of cores to use (e.g. 12)')
parser.add_argument('-t', '--tmp_hail', required=True, help='path to your hail temporary folder')

args = parser.parse_args()

hl.init(
    master=f"local[{args.cores}]",
    spark_conf={
        "spark.local.dir": args.tmp_hail,
        "spark.driver.extraJavaOptions": f"-Djava.io.tmpdir={args.tmp_hail}",
        "spark.driver.memory": args.driver_memory,
        "spark.executor.memory": args.executor_memory
    }
)

gnomad = hl.import_table('data/gnomad.exomes.coverage.summary.tsv.bgz')
gnomad = gnomad.drop('median',
                     'over_1',
                     'over_5',
                     'over_15',
                     'over_20',
                     'over_25',
                     'over_30',
                     'over_50',
                     'over_100',
                     'mean')
gnomad = gnomad.filter((gnomad.chrom == 'X') | (gnomad.chrom == 'Y') | (gnomad.chrom == 'MT'), keep = False)
gnomad = gnomad.annotate(chrom = hl.int(gnomad.chrom))
gnomad = gnomad.annotate(over_10 = hl.float(gnomad.over_10))
gnomad_flt = gnomad.filter(gnomad.over_10 >= 0.9, keep = True)
# Convert to Pandas
gnomad_pd = gnomad_flt.to_pandas()
gnomad_pd.pos = pd.to_numeric(gnomad_pd.pos)
gnomad_pd.chrom = pd.to_numeric(gnomad_pd.chrom)
# Sort values by chromosome and start site
gnomad_pd = gnomad_pd.sort_values(["chrom", "pos"]).reset_index(drop=True)

def merge_vectorized(sub):
    # Create a vector of "previous" end value which is initialized by 0
    prev_pos = sub["pos"].shift(fill_value=0)
    # Breaks if start is higher than prev_end + 1 (gets true)
    break_mask = sub["pos"] > (prev_pos + 1)
    # group_id is being incremented at each break (true =+1; false =+0)
    group_id = break_mask.cumsum()
    # Aggregate by group: start = min(start), end = max(end)
    merged = sub.groupby(group_id, as_index=False).agg(
        chrom=("chrom","first"),
        start=("pos","min"),
        end=("pos","max")
    )
    return merged[["chrom","start","end"]]

pseudoBed = (
    gnomad_pd.groupby("chrom", group_keys=False)
      .apply(merge_vectorized)
      .sort_values(["chrom", "start"])
)

pseudoBed.to_csv('intermediate_files/filtered_coverage_gnomad.tsv', header=True, sep="\t", index=False)
