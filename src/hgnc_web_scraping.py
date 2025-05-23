# -*- coding: utf-8 -*-
# antoine.favier@institutimagine.org
# Script to download HGNC gene symbols & other info

import urllib, os
import pandas as pd
data = urllib.request.urlopen('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt')
tab = pd.read_csv(data, sep='\t', header=0, low_memory=False)
tab.to_csv(path_or_buf=os.getcwd()+'/data/hgnc_protein_coding.txt', sep='\t', na_rep='NA', index=False)
