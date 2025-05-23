#!/data/antoine_data/anaconda2/bin/python2
# -*- coding: utf-8 -*-
# antoine.favier@institutimagine.org
# script to calculate one probability per gene

import os
import pandas as pd
import numpy as np

predictfile = os.getcwd()+'/data/gnomad_all_chr_predicted.tsv'
data = pd.read_csv(predictfile, sep="\t")
genelist = []
for gene in data.geneID:
	if gene not in genelist:
		genelist.append(gene)

pergenefile = os.getcwd()+'/data/gnomAD_count_table.tsv'
with open(pergenefile, 'w') as o:
    o.write("%s\n" % '\t'.join(['chromosome',
                                'transcript',
                                'geneID',
                                'length',
                                'obs_syn',
                                'obs_mis',
                                'obs_ptv'
                                ]))
    for gene in genelist:
        chrom = str(list(data.loc[data.geneID == gene].chr)[0])
        trans = str(list(data.loc[data.geneID == gene].transcript)[0])
        length = str(sum(data.loc[data.geneID == gene].Region_Length))
        tot_syn = str(np.nansum(data.loc[data.geneID == gene].Region_syn_obs))
        tot_mis = str(np.nansum(data.loc[data.geneID == gene].Region_mis_obs))
        tot_ptv = str(np.nansum(data.loc[data.geneID == gene].Region_ptv_obs))
        o.write("%s\n" % '\t'.join([chrom,
                                    trans,
                                    gene,
                                    length,
                                    tot_syn,
                                    tot_mis,
                                    tot_ptv]))
