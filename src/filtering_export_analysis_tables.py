# -*- coding: utf-8 -*-
# antoine.favier@institutimagine.org
# Keep only 1 transcript (APPRIS PRINCIPAL) per variant and only synonymous and missense
# Keep only variants at low frequency both in the patients cohort and in gnomAD
# And export all the needed tables for analysis

import os, sys
import hail as hl
import numpy as np
import pandas as pd

argvL = sys.argv
max_mem = str(argvL[1])
hail_tmp = str(argvL[2])

os.environ['PYSPARK_SUBMIT_ARGS'] = "--driver-memory "+max_mem+" pyspark-shell"
os.environ['TMPDIR'] = hail_tmp
hl.init(spark_conf={"spark.local.<200b>​dir":hail_tmp,
                    "spark.driver.extraJavaOptions":"-Djava.io.tmpdir="+hail_tmp})

# import the list of APPRIS PRINCIPAL transcripts
transcriptlist = os.getcwd()+'/data/transcript_list.tsv'

# Create a list of all APPRIS principal transcripts
appris = []
with open(transcriptlist) as tab:
    for t in tab:
        appris.append(t.lstrip().split()[1])
appris = hl.literal(appris)

# Create a dictionary to change the chromosome name while importing the vcf
dict_chr=dict()
for x in [str(x) for x in range(1,23)] + ['X', 'Y']:
    dict_chr['chr'+x]=x

# Import the 1st row of the vcf
vcf = pd.read_csv(os.getcwd()+'/data/patients_annotated.vcf.bgz',
                  compression='gzip', comment='#', chunksize=1,
                  delim_whitespace=True, header=None)
# Check if the CHR column has only numbers or a 'chr' character and import the
# annotated vcf as a hail matrix table
chrTest = pd.DataFrame(vcf.get_chunk(1)).iloc[0,0]
if chrTest[0:3] == 'chr':
    patients = hl.import_vcf(os.getcwd()+'/data/patients_annotated.vcf.bgz',
                              contig_recoding=dict_chr,
                              reference_genome='GRCh37', force_bgz=True)
else:
    patients = hl.import_vcf(os.getcwd()+'/data/patients_annotated.vcf.bgz',
                          reference_genome='GRCh37', force_bgz=True)

# Keep only "trustworthy well-covered regions" from patients and gnomAD
bed = hl.import_bed(os.getcwd()+'/data/gnomad_patients_all_chr.bed', reference_genome='GRCh37')
regions = [row['interval'] for row in bed.collect()]
patients_flt_reg = hl.filter_intervals(patients, regions)

# Update AC, AF, AN, NS after removing non-european samples
patients_flt_reg = hl.variant_qc(patients_flt_reg)
patients_flt_reg = patients_flt_reg.annotate_rows(
    info = patients_flt_reg.info.annotate(AC=patients_flt_reg.variant_qc.AC[1],
                                          AF=patients_flt_reg.variant_qc.AF[1],
                                          AN=patients_flt_reg.variant_qc.AN))
# Remove AC=0 variants (that may appear when removing non European samples)
patients_cor = patients_flt_reg.filter_rows(patients_flt_reg.info.AC > 0, keep=True)

# Genotype filtering
patients_cor = patients_cor.filter_rows(hl.agg.mean(patients_cor.GQ) > 20.0)
patients_cor = patients_cor.filter_rows(hl.agg.mean(patients_cor.DP) > 8)
 
# filter out samples < 95% call rate
patients_cor = hl.sample_qc(patients_cor, name='sample_qc')
patients_prun = patients_cor.filter_cols(patients_cor.sample_qc.call_rate < 0.95, keep=False)

# Re-update AC, AF, AN, NS after pruning samples
patients_prun = hl.variant_qc(patients_prun)
patients_prun = patients_prun.annotate_rows(
    info = patients_prun.info.annotate(AC=patients_prun.variant_qc.AC[1],
                                       AF=patients_prun.variant_qc.AF[1],
                                       AN=patients_prun.variant_qc.AN))
# Remove AC=0 variants (that may appear when pruning samples)
patients_prun = patients_prun.filter_rows(patients_prun.info.AC > 0, keep=True)

##### KEEP ONLY ONE TRANSCRIPT PER VARIANT AND FILTER #####

# Only keep the consequences for which the transcript is APPRIS PRINCIPAL
patients_prun = patients_prun.annotate_rows(info = patients_prun.info.annotate(CSQ = hl.filter(
    lambda x: appris.contains((x.split('\\|')[6])), patients_prun.info.CSQ)))

# Add a column with the consequences of the transcripts kept
patients_prun = patients_prun.annotate_rows(info = patients_prun.info.annotate(conseq = hl.map(
    lambda x: x.split('\\|')[1], patients_prun.info.CSQ)))

# Function to prioritize coding consequences (PTV > missense > synonymous)
prot_truncating = ["stop_gained","splice_acceptor_variant",
                   "splice_donor_variant", "frameshift_variant"]
prot_truncating = hl.literal(prot_truncating)
def Prioritize_CSQ(x):
    expr = (hl.case()
            .when(x == 'synonymous_variant', 1)
            .when(x == 'missense_variant', 2)
            .when(prot_truncating.contains(x), 3)
            .or_missing())
    return expr

# Create a new column with the score of prioritization of the consequences
patients_prun = patients_prun.annotate_rows(info = patients_prun.info.annotate(prio = hl.map(
    lambda x: Prioritize_CSQ(x), patients_prun.info.conseq)))
# Create a new column with index of the highest prioritization score
patients_prun = patients_prun.annotate_rows(info = patients_prun.info.annotate(
    prio_idx = patients_prun.info.prio.index(hl.nanmax(patients_prun.info.prio))))
# Remove variants without any synonymous, missense or PTV coding consequence
patients_prun = patients_prun.filter_rows(hl.is_nan(patients_prun.info.prio_idx), keep = False)
# Keep only the transcript with the highest prioritization score
patients_prun = patients_prun.annotate_rows(info = patients_prun.info.annotate(
    CSQ = patients_prun.info.CSQ[patients_prun.info.prio_idx]))
patients_prun = patients_prun.annotate_rows(info = patients_prun.info.drop('prio', 'prio_idx'))
# Replace conseq with the only remaining conseq
patients_prun = patients_prun.annotate_rows(info = patients_prun.info.annotate(
    conseq = patients_prun.info.CSQ.split('\\|')[1]))

# Add the gnomAD global Allele Frequency
gnomad = hl.read_table('/data-cbl/gnomad/gnomad_hailtables/gnomad.exomes.r2.1.1.sites.ht')
gnomad = gnomad.key_by('locus')
gnomad_index = gnomad.freq_index_dict['gnomad'].collect()[0]
patients_gnomad = patients_prun.annotate_rows(freq = gnomad[patients_prun.locus].freq)
# Only keep variants with gnomAD AF <= 0.1% (keep also missing data from gnomAD)
patients_gnomad = patients_gnomad.annotate_rows(info = patients_gnomad.info.annotate(
    AF_gnomAD = patients_gnomad.freq.AF[gnomad_index]))
patients_gnomad = patients_gnomad.filter_rows((patients_gnomad.info.AF_gnomAD <= 0.001) |
                                              (hl.is_missing(patients_gnomad.info.AF_gnomAD)), keep = True)

# Keep only "PASS" filters from the VCF and AF <= 1% in the patients cohort
patients_gnomad = patients_gnomad.filter_rows((patients_gnomad.info.AF <= 0.01), keep = True)

# Extract ENSG and gene symbol
patients_gnomad = patients_gnomad.annotate_rows(info = patients_gnomad.info.annotate(
    ensg = patients_gnomad.info.CSQ.split('\\|')[4]))
patients_gnomad = patients_gnomad.annotate_rows(info = patients_gnomad.info.annotate(
    genesymbol = patients_gnomad.info.CSQ.split('\\|')[3]))

# Export the matrix table
patients_gnomad.write(os.getcwd()+'/data/hail_matrix_patients.mt')

##### EXPORT THE TABLES #####

# Create a dictionary with the number of variants for each gene
var_nbr = patients_gnomad.aggregate_rows(
    hl.agg.group_by(hl.struct(chrom = patients_gnomad.locus.contig,
                              ensg=patients_gnomad.info.ensg,
                              genesymbol=patients_gnomad.info.genesymbol),
                    hl.struct(
                        total = hl.agg.count(),
                        num_syn = hl.agg.count_where(patients_gnomad.info.conseq == 'synonymous_variant'),
                        num_mis = hl.agg.count_where(patients_gnomad.info.conseq == 'missense_variant'),
                        num_PTV = hl.agg.count_where(prot_truncating.contains(patients_gnomad.info.conseq)))))
# Export the dictionary in a table
with open(os.getcwd()+'/data/variant_number.tsv', 'w') as out:
    print('Writing Variant number data table')
    out.write('Chromosome\tGene_ID\tGene_symbol\tVariant_number\tSynonymous_variants\tMissense_variants\tPTV\n')
    for key, value in var_nbr.items():
        out.write('%s\n' % '\t'.join(map(str, [
            key.chrom,
            key.ensg,
            key.genesymbol,
            value.total,
            value.num_syn,
            value.num_mis,
            value.num_PTV])))

# Create a table for synonymous variants, for missense variants and PTV
data_syn = patients_gnomad.filter_rows(patients_gnomad.info.conseq == 'synonymous_variant', keep = True)
data_mis = patients_gnomad.filter_rows(patients_gnomad.info.conseq == 'missense_variant', keep = True)
data_ptv = patients_gnomad.filter_rows(prot_truncating.contains(patients_gnomad.info.conseq), keep = True)

# Create a table with numbers of homo/heterozygote synonymous/missense variants for each sample/gene couple
# First let's create 4 matrix tables with the counts of variants per sample/gene couple
hetero_syn = data_syn.group_rows_by(ensg=data_syn.info.ensg).aggregate(
        n=hl.agg.count_where(data_syn.GT.is_het()))
homo_syn = data_syn.group_rows_by(ensg=data_syn.info.ensg).aggregate(
    n=hl.agg.count_where(data_syn.GT.is_hom_var()))
hetero_mis = data_mis.group_rows_by(ensg=data_mis.info.ensg).aggregate(
        n=hl.agg.count_where(data_mis.GT.is_het()))
homo_mis = data_mis.group_rows_by(ensg=data_mis.info.ensg).aggregate(
    n=hl.agg.count_where(data_mis.GT.is_hom_var()))
hetero_ptv = data_ptv.group_rows_by(ensg=data_ptv.info.ensg).aggregate(
        n=hl.agg.count_where(data_ptv.GT.is_het()))
homo_ptv = data_ptv.group_rows_by(ensg=data_ptv.info.ensg).aggregate(
    n=hl.agg.count_where(data_ptv.GT.is_hom_var()))
# Function to create the four the 4 tables to export
def tableExport(matTable, tableName):
    t = matTable.localize_entries('entry_structs', 'columns')
    t = t.select(entries = t.entry_structs.map(lambda entry: entry.n))
    numpy_data = np.array(t.entries.collect())
    ensg = matTable.rows().key_by().to_pandas()
    samples = matTable.cols().key_by().to_pandas()
    df = pd.DataFrame(data=numpy_data, index=ensg.ensg.to_list(), columns=samples.s.to_list())
    df.to_csv(os.getcwd()+'/data/'+tableName+'.tsv', sep='\t', header=True, index=True)
# Export the tables
tableExport(hetero_syn, 'Hetero_syn_var')
tableExport(homo_syn, 'Homo_syn_var')
tableExport(hetero_mis, 'Hetero_mis_var')
tableExport(homo_mis, 'Homo_mis_var')
tableExport(hetero_ptv, 'Hetero_ptv_var')
tableExport(homo_ptv, 'Homo_ptv_var')

# Create a dictionary with counts for each sample the number of synonymous/missense/PTV hetero/homozygote variants
samp_var_nbr = patients_gnomad.aggregate_entries(
    hl.agg.group_by(hl.struct(samples = patients_gnomad.s),
    hl.struct(
        syn_hom = hl.agg.count_where((patients_gnomad.info.conseq == 'synonymous_variant') &
                                     (patients_gnomad.GT.is_hom_var())),
        syn_het = hl.agg.count_where((patients_gnomad.info.conseq == 'synonymous_variant') &
                                     (patients_gnomad.GT.is_het())),
        mis_hom = hl.agg.count_where((patients_gnomad.info.conseq == 'missense_variant') &
                                     (patients_gnomad.GT.is_hom_var())),
        mis_het = hl.agg.count_where((patients_gnomad.info.conseq == 'missense_variant') &
                                     (patients_gnomad.GT.is_het())),
        ptv_hom = hl.agg.count_where(prot_truncating.contains(patients_gnomad.info.conseq) &
                                     (patients_gnomad.GT.is_hom_var())),
        ptv_het = hl.agg.count_where(prot_truncating.contains(patients_gnomad.info.conseq) &
                                     (patients_gnomad.GT.is_het())))))
# Export the dictionary in a table
with open(os.getcwd()+'/data/variant_per_individual.tsv', 'w') as out:
    print('Writing Variant number per sample data table')
    out.write('\tsyn_homo\tsyn_hetero\tmis_homo\tmis_hetero\tptv_homo\tptv_hetero\n')
    for key, value in samp_var_nbr.items():
        out.write('%s\n' % '\t'.join(map(str, [
            key.samples,
            value.syn_hom,
            value.syn_het,
            value.mis_hom,
            value.mis_het,
            value.ptv_hom,
            value.ptv_het])))
