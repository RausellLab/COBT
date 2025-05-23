# -*- coding: utf-8 -*-
# antoine.favier@institutimagine.org
# Only keep regions covered in at least a median of 10X in gnomAD

import sys, os, pandas as pd
import hail as hl
from joblib import Parallel, delayed

argvL = sys.argv
core_nbr = int(argvL[1])
max_mem = str(argvL[2])
hail_tmp = str(argvL[3])

os.environ['PYSPARK_SUBMIT_ARGS'] = "--driver-memory "+max_mem+" pyspark-shell"
os.environ['TMPDIR'] = hail_tmp
hl.init(spark_conf={"spark.local.<200b>​dir":hail_tmp,
                    "spark.driver.extraJavaOptions":"-Djava.io.tmpdir="+hail_tmp})


gnomad = hl.import_table('data/gnomad.exomes.coverage.summary.tsv.bgz', )
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
gnomad = gnomad.filter(gnomad.over_10 >= 0.9, keep = True)
gnomad_flt = gnomad.to_pandas()
gnomad_flt.pos = pd.to_numeric(gnomad_flt.pos)
gnomad_flt.chrom = pd.to_numeric(gnomad_flt.chrom)

chromosome = range(1,23)

def gnomad_filter(chro):
    chrom = []
    start = []
    end = []
    chromToKeep = (gnomad_flt['chrom'] == chro)
    gnomad_chr = gnomad_flt[chromToKeep]
    for row in gnomad_chr.index:
         gnomad_chr.loc[row, "pos"]
         if start == []:
             chrom.append(chro)
             start.append(gnomad_chr.loc[row, "pos"])
             posi = gnomad_chr.loc[row, "pos"]
         elif gnomad_chr.loc[row, "pos"] == posi+1:
             posi = gnomad_chr.loc[row, "pos"]
         elif gnomad_chr.loc[row, "pos"] > posi+1:
             end.append(posi)
             chrom.append(chro)
             start.append(gnomad_chr.loc[row, "pos"])
             posi = gnomad_chr.loc[row, "pos"]
    if len(start) > len(end):
        end.append(posi)
    df = pd.DataFrame({'chr' : chrom,
                       'start' : start,
                       'end' : end})
    df = df[['chr','start','end']]
    df2 = df.sort_values(by=['chr', 'start'])
    df2.to_csv(os.getcwd()+'/intermediate_files/TEMP/filtered_coverage_gnomad_chr{}.tsv'.format(chro),
               header=True, sep="\t", index=False)

Parallel(n_jobs=core_nbr)(delayed(gnomad_filter)(chro) for chro in chromosome) # loop for each chromosome with parallelization

# Concatenate all files
dff = pd.read_csv(os.getcwd()+'/intermediate_files/TEMP/filtered_coverage_gnomad_chr1.tsv', sep='\t')
for chro in chromosome[1:len(chromosome)]:
    dff2 = pd.read_csv(os.getcwd()+'/intermediate_files/TEMP/filtered_coverage_gnomad_chr{}.tsv'.format(chro),
                        sep='\t', dtype = {"Chromosome" : "str"})
    dff = pd.concat([dff,dff2])
dff.to_csv('intermediate_files/filtered_coverage_gnomad.tsv', header=True, sep="\t", index=False)
