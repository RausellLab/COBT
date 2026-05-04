#!/usr/bin/env python3
# Keep only regions with at least 10X of median coverage and reformat chromosomes

import os, sys
import pandas as pd

argvL = sys.argv
cov_file = str(argvL[1])
nbr_samples = int(argvL[2])

# Import the coverage file
patientsCov = pd.read_csv(os.getcwd()+'/'+cov_file, sep="\t", header=None,
    names=["chrom", "start", "end", "cov10X"],
    dtype={"chrom": "str", "start": "int32", "end": "int32", "cov10X": "int32"},
    compression="infer"
)
# Filter out chromosomes X, Y and MT
patientsCov = patientsCov[~patientsCov["chrom"].isin(["X", "Y", "MT"])].copy()
# Only keep regions with more than 90% of samples sequenced at 10X
patientsCov = patientsCov[(patientsCov["cov10X"]/nbr_samples) >= 0.9].copy()
# Convert chromosomes to int to get them in the right order
patientsCov["chrom"] = patientsCov["chrom"].astype("int")
# Sort values by chromosome and start site
patientsCov = patientsCov.sort_values(["chrom", "start"]).reset_index(drop=True)

def merge_vectorized(sub):
    # Create a vector of "previous" end value which is initialized by 0
    prev_end = sub["end"].shift(fill_value=0)
    # Breaks if start is higher than prev_end + 1 (gets true)
    break_mask = sub["start"] > (prev_end + 1)
    # group_id is being incremented at each break (true =+1; false =+0)
    group_id = break_mask.cumsum()
    # Aggregate by group: start = min(start), end = max(end)
    merged = sub.groupby(group_id, as_index=False).agg(
        chrom=("chrom", "first"),
        start=("start", "min"),
        end=("end", "max")
    )
    return merged[["chrom", "start", "end"]]

pseudoBed = (
    patientsCov.groupby("chrom", group_keys=False)
      .apply(merge_vectorized)
      .sort_values(["chrom", "start"])
)

# Export pseudo-bed file
pseudoBed.to_csv(os.getcwd()+'/intermediate_files/filtered_coverage_patients.tsv',
    sep="\t",
    header=True,
    index=False
)

