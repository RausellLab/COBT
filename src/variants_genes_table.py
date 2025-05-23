# -*- coding: utf-8 -*-
# antoine.favier@institutimagine.org
# Export variant tables for all variants kept in the after case-only pipeline

import os, argparse
import hail as hl
import pandas as pd
import numpy as np

descr = "variants_genes_table.py \n"
descr += "usage: python src/variants_genes_table.py -m 100g -t /data-tmp/antoine_data/hailTMP \n"
parser = argparse.ArgumentParser(description=descr)
parser.add_argument('-m', '--max_memory', required=True, help='maximum memory to allocate')
parser.add_argument('-t', '--tmp_hail', required=True, help='path to your hail temporary folder')

args = parser.parse_args()
max_mem = args.max_memory
hail_tmp = args.tmp_hail

os.environ['PYSPARK_SUBMIT_ARGS'] = "--driver-memory "+max_mem+" pyspark-shell"
os.environ['TMPDIR'] = hail_tmp
hl.init(spark_conf={"spark.local.<200b>​dir":hail_tmp,
                    "spark.driver.extraJavaOptions":"-Djava.io.tmpdir="+hail_tmp})

# Load the matrix table
mat = hl.read_matrix_table(os.getcwd()+'/data/hail_matrix_patients.mt')

# Extract transcript information
mat = mat.annotate_rows(info = mat.info.annotate(transcript = mat.info.CSQ.split('\\|')[6]))

prot_truncating = ["stop_gained","splice_acceptor_variant",
                   "splice_donor_variant", "frameshift_variant"]
prot_truncating = hl.literal(prot_truncating)

# Load hits table
syn_hits = pd.read_csv(os.getcwd()+'/results/synonymous_hits.tsv', sep='\t')
mis_hits = pd.read_csv(os.getcwd()+'/results/missense_hits.tsv', sep='\t')
ptv_hits = pd.read_csv(os.getcwd()+'/results/ptv_hits.tsv', sep='\t')
alt_hits = pd.read_csv(os.getcwd()+'/results/protein_altering_hits.tsv', sep='\t')

# Lists of hit genes for each variant consequence
synGenes = hl.literal(syn_hits['ENSG'].tolist())
misGenes = hl.literal(mis_hits['ENSG'].tolist())
ptvGenes = hl.literal(ptv_hits['ENSG'].tolist())
altGenes = hl.literal(alt_hits['ENSG'].tolist())

# filter for missense variants
matsyn = mat.filter_rows(synGenes.contains(mat.info.ensg), keep=True)
matsyn = matsyn.filter_rows(matsyn.info.conseq == "synonymous_variant", keep=True)
matmis = mat.filter_rows(misGenes.contains(mat.info.ensg), keep=True)
matmis = matmis.filter_rows(matmis.info.conseq == "missense_variant", keep=True)
matptv = mat.filter_rows(ptvGenes.contains(mat.info.ensg), keep=True)
matptv = matptv.filter_rows(prot_truncating.contains(matptv.info.conseq), keep=True)
matalt = mat.filter_rows(altGenes.contains(mat.info.ensg), keep=True)
matalt = matalt.filter_rows(((matalt.info.conseq == "missense_variant") |
                             prot_truncating.contains(matalt.info.conseq)), keep=True)

# Function to export the tables
def export_variant_table(matcons):
    mat_het = matcons.group_rows_by(variant=matcons.locus,
                                    rsid=matcons.rsid,
                                    alleles=matcons.alleles,
                                    consequence=matcons.info.conseq,
                                    geneSymbol=matcons.info.genesymbol,
                                    transcript=matcons.info.transcript).aggregate(
                                        n=hl.agg.count_where(matcons.GT.is_het()))
    mat_hom = matcons.group_rows_by(variant=matcons.locus).aggregate(
        n=hl.agg.count_where(matcons.GT.is_hom_var()))
    t_het = mat_het.localize_entries('entry_structs', 'columns')
    t_het = t_het.select(entries = t_het.entry_structs.map(lambda entry: entry.n))
    het_data = np.array(t_het.entries.collect())
    t_hom = mat_hom.localize_entries('entry_structs', 'columns')
    t_hom = t_hom.select(entries = t_hom.entry_structs.map(lambda entry: entry.n))
    hom_data = np.array(t_hom.entries.collect())
    samples = mat_het.cols().key_by().to_pandas()
    variants = mat_het.rows().key_by().to_pandas()
    variants = variants.rename(columns={"variant.contig":"chrom", "variant.position":"pos"})
    het = pd.DataFrame(data=het_data, columns=samples.s.to_list())
    het = het.apply(lambda row: row[row != 0].index, axis=1)
    het.name = "heterozygous_mutated_patients"
    hom = pd.DataFrame(data=hom_data, columns=samples.s.to_list())
    hom = hom.apply(lambda row: row[row != 0].index, axis=1)
    hom.name = "homozygous_mutated_patients"
    returnTable = pd.concat([variants, het, hom], axis=1)
    return(returnTable)

exportSYN = export_variant_table(matsyn)
exportSYN['alleles'] = [','.join(map(str, l)) for l in exportSYN['alleles']]
exportSYN['heterozygous_mutated_patients'] = [','.join(map(str, l)) for l in exportSYN['heterozygous_mutated_patients']]
exportSYN['heterozygous_mutated_patients'] = exportSYN['heterozygous_mutated_patients'].replace('', np.nan)
exportSYN['homozygous_mutated_patients'] = [','.join(map(str, l)) if not type(l) == float else np.nan for l in exportSYN['homozygous_mutated_patients']]
exportSYN['homozygous_mutated_patients'] = exportSYN['homozygous_mutated_patients'].replace('', np.nan)
exportSYN.to_csv('results/synonymous_variants_case_only.tsv', sep='\t', na_rep='NA', index=False)

exportMIS = export_variant_table(matmis)
exportMIS['alleles'] = [','.join(map(str, l)) for l in exportMIS['alleles']]
exportMIS['heterozygous_mutated_patients'] = [','.join(map(str, l)) for l in exportMIS['heterozygous_mutated_patients']]
exportMIS['heterozygous_mutated_patients'] = exportMIS['heterozygous_mutated_patients'].replace('', np.nan)
exportMIS['homozygous_mutated_patients'] = [','.join(map(str, l)) if not type(l) == float else np.nan for l in exportMIS['homozygous_mutated_patients']]
exportMIS['homozygous_mutated_patients'] = exportMIS['homozygous_mutated_patients'].replace('', np.nan)
exportMIS.to_csv('results/missense_variants_case_only.tsv', sep='\t', na_rep='NA', index=False)

exportPTV = export_variant_table(matptv)
exportPTV['alleles'] = [','.join(map(str, l)) for l in exportPTV['alleles']]
exportPTV['heterozygous_mutated_patients'] = [','.join(map(str, l)) for l in exportPTV['heterozygous_mutated_patients']]
exportPTV['heterozygous_mutated_patients'] = exportPTV['heterozygous_mutated_patients'].replace('', np.nan)
exportPTV['homozygous_mutated_patients'] = [','.join(map(str, l)) if not type(l) == float else np.nan for l in exportPTV['homozygous_mutated_patients']]
exportPTV['homozygous_mutated_patients'] = exportPTV['homozygous_mutated_patients'].replace('', np.nan)
exportPTV.to_csv('results/protein_truncating_variants_case_only.tsv', sep='\t', na_rep='NA', index=False)

exportALT = export_variant_table(matalt)
exportALT['alleles'] = [','.join(map(str, l)) for l in exportALT['alleles']]
exportALT['heterozygous_mutated_patients'] = [','.join(map(str, l)) for l in exportALT['heterozygous_mutated_patients']]
exportALT['heterozygous_mutated_patients'] = exportALT['heterozygous_mutated_patients'].replace('', np.nan)
exportALT['homozygous_mutated_patients'] = [','.join(map(str, l)) if not type(l) == float else np.nan for l in exportALT['homozygous_mutated_patients']]
exportALT['homozygous_mutated_patients'] = exportALT['homozygous_mutated_patients'].replace('', np.nan)
exportALT.to_csv('results/protein_altering_variants_case_only.tsv', sep='\t', na_rep='NA', index=False)
