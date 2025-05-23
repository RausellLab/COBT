# -*- coding: utf-8 -*-
# antoine.favier@institutimagine.org
# Keep only regions with at least 10X of median coverage and reformat chromosomes

import os, sys
import hail as hl
import pandas as pd
from joblib import Parallel, delayed

argvL = sys.argv
core_nbr = int(argvL[1])
max_mem = str(argvL[2])
hail_tmp = str(argvL[3])
cov_file = str(argvL[4])

os.environ['PYSPARK_SUBMIT_ARGS'] = "--driver-memory "+max_mem+" pyspark-shell"
os.environ['TMPDIR'] = hail_tmp
hl.init(spark_conf={"spark.local.<200b>​dir":hail_tmp,
                    "spark.driver.extraJavaOptions":"-Djava.io.tmpdir="+hail_tmp})

# Import the coverage file
if (cov_file.endswith('.gz') | cov_file.endswith('.bgz')):
    patientsCov = hl.import_table(cov_file,
                                  types={'f0': hl.tstr, 'f1': hl.tint32,
                                         'f2':hl.tint32, 'f3':hl.tint32},
                                  no_header=True, force_bgz=True)
else:
    patientsCov = hl.import_table(cov_file,
                              types={'f0': hl.tstr, 'f1': hl.tint32,
                                     'f2':hl.tint32, 'f3':hl.tint32},
                              no_header=True)

patientsCov = patientsCov.rename({'f0':'chrom', 'f1':'start', 'f2':'end', 'f3':'median'})
# Filter X, Y and MT chromosomes
patientsCov = patientsCov.filter((patientsCov.chrom == 'X') | (patientsCov.chrom == 'Y') |
                                 (patientsCov.chrom == 'MT'), keep = False)
patientsCov = patientsCov.filter((patientsCov.median/497) >= 0.9, keep = True)
# Export as a pandas dataframe
patientsCov10X = patientsCov.to_pandas()
patientsCov10X.chrom = pd.to_numeric(patientsCov10X.chrom)

chromosome = range(1,23)

def table_filter(chro):
    start = []
    end = []
    chrom = []
    chromToKeep = (patientsCov10X['chrom'] == chro)
    patientsCov10X_chr = patientsCov10X[chromToKeep]
    for row in patientsCov10X_chr.index:
        if start == []:
            starti = patientsCov10X_chr.loc[row, "start"]
            endi = patientsCov10X_chr.loc[row, "end"]
            chrom.append(chro)
            start.append(starti)
        elif patientsCov10X_chr.loc[row, "start"] == endi+1:
            endi = patientsCov10X_chr.loc[row, "end"]
        elif patientsCov10X_chr.loc[row, "start"] > endi+1:
            end.append(endi)
            starti=patientsCov10X_chr.loc[row, "start"]
            endi=patientsCov10X_chr.loc[row, "end"]
            chrom.append(chro)
            start.append(starti)
    if len(start) > len(end):
        end.append(patientsCov10X_chr.loc[row, "end"])
    df = pd.DataFrame({'chr' : chrom,
                       'start' : start,
                       'end' : end})
    df = df[['chr','start','end']]
    df2 = df.sort_values(by=['chr', 'start'])
    df2.to_csv(os.getcwd()+'/intermediate_files/TEMP/filtered_coverage_patients_chr{}.tsv'.format(chro),
               header=True, sep="\t", index=False)

Parallel(n_jobs=core_nbr)(delayed(table_filter)(chro) for chro in chromosome) # loop for each chromosome with parallelization

# Concatenate all files
dff = pd.read_csv(os.getcwd()+'/intermediate_files/TEMP/filtered_coverage_patients_chr1.tsv',sep='\t')
for chro in chromosome[1:len(chromosome)]:
    dff2 = pd.read_csv(os.getcwd()+'/intermediate_files/TEMP/filtered_coverage_patients_chr{}.tsv'.format(chro),
                        sep='\t', dtype = {"Chromosome" : "str"})
    dff = pd.concat([dff,dff2])
dff.to_csv('intermediate_files/filtered_coverage_patients.tsv',header=True,sep="\t",index=False)
