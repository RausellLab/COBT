# -*- coding: utf-8 -*-
# antoine.favier@institutimagine.org
# Script to compute in gnomAD for all samples
# the counts of synonymous missense and PTVs variants in the pre-defined regions

import hail as hl
import argparse, os

descr = "gnomAD_count_proba.py \n"
descr += "usage: python src/gnomAD_count_proba.py -m 100g -t /data-tmp/antoine_data/hailTMP -i data/gnomad_all_chr.bed -o data/gnomad_all_chr_predicted.bed \n"
parser = argparse.ArgumentParser(description=descr)
parser.add_argument('-m', '--max_memory', required=True, help='maximum memory to allocate')
parser.add_argument('-t', '--tmp_hail', required=True, help='path to your hail temporary folder')
parser.add_argument('-i', '--subRegionF', required=True, help='input subRegion file')
parser.add_argument('-o', '--output', required=True, help='output expected probability tsv file')

args = parser.parse_args()
max_mem = args.max_memory
hail_tmp = args.tmp_hail
subRegionF = args.subRegionF
output = args.output

os.environ['PYSPARK_SUBMIT_ARGS'] = "--driver-memory "+max_mem+" pyspark-shell"
os.environ['TMPDIR'] = hail_tmp
hl.init(spark_conf={"spark.local.<200b>​dir":hail_tmp,
                    "spark.driver.extraJavaOptions":"-Djava.io.tmpdir="+hail_tmp})

max_AN = 125748*2 # Number of gnomAD samples in exomes overall

# Import gnomad table
gnomad = hl.read_table('data/gnomad.exomes.r2.1.1.sites.ht')
# Keep only vep annotation and freq
gnomad = gnomad.drop('allele_info',
                     'qual',
                     'age_hist_het',
                     'age_hist_hom',
                     'popmax',
                     'faf',
                     'rsid',
                     'lcr',
                     'decoy',
                     'segdup', 
                     'nonpar', 
                     'variant_type', 
                     'allele_type', 
                     'n_alt_alleles', 
                     'was_mixed', 
                     'has_star', 
                     'qd', 
                     'pab_max', 
                     'info_MQRankSum', 
                     'info_SOR', 
                     'info_InbreedingCoeff', 
                     'info_ReadPosRankSum', 
                     'info_FS', 
                     'info_QD', 
                     'info_MQ', 
                     'info_DP', 
                     'transmitted_singleton', 
                     'fail_hard_filters', 
                     'info_POSITIVE_TRAIN_SITE', 
                     'info_NEGATIVE_TRAIN_SITE', 
                     'omni', 
                     'mills', 
                     'n_nonref', 
                     'tp', 
                     'rf_train', 
                     'rf_label', 
                     'rf_probability', 
                     'singleton', 
                     'was_split', 
                     'score', 
                     'rank', 
                     'singleton_rank', 
                     'biallelic_rank', 
                     'adj_biallelic_singleton_rank', 
                     'adj_rank', 
                     'adj_biallelic_rank', 
                     'adj_singleton_rank', 
                     'biallelic_singleton_rank', 
                     'filters', 
                     'gq_hist_alt',
                     'gq_hist_all',
                     'dp_hist_alt',
                     'dp_hist_all',
                     'ab_hist_alt')
gnomad = gnomad.key_by('locus')

prot_truncating = ["stop_gained","splice_acceptor_variant",
                   "splice_donor_variant", "frameshift_variant"]
prot_truncating = hl.literal(prot_truncating)

# Import bed region file
reg = hl.import_bed(subRegionF, reference_genome='GRCh37')
colToSplit = reg.target.split(':')
reg = reg.annotate(transcript = colToSplit[0],
                   exon = hl.int(colToSplit[1]),
                   subRegionID = colToSplit[3],
                   subRegionReg = ([colToSplit[3]+'|'+
                                    colToSplit[4].split('-')[0]+'|'+
                                    colToSplit[4].split('-')[1]]),
                   regionLength = reg.interval.end.position-reg.interval.start.position,
                   geneID = colToSplit[5],
                   geneSymb = colToSplit[6])
reg = reg.drop('target')
# Transform intervals in loci to join with gnomad table
reg = reg.annotate(start = reg.interval.start.position)
reg = reg.annotate(end = reg.interval.end.position)
reg = reg.annotate(chrom = reg.interval.end.contig)
reg = reg.annotate(pos = hl.range(reg.start,reg.end))
reg = reg.explode(reg.pos)
reg = reg.annotate(locus = hl.locus(reg.chrom, reg.pos, reference_genome='GRCh37'))
reg = reg.drop('pos', 'chrom', 'start', 'end')
reg = reg.key_by('locus')

# Joining region table and gnomad table
reg_gnomad = reg.annotate(vep = gnomad[reg.locus].vep, freq = gnomad[reg.locus].freq)

gnomad_index = gnomad.freq_index_dict['gnomad'].collect()[0]
# Filter poorly sequenced variants
reg_gnomad_flt = reg_gnomad.filter((reg_gnomad.freq.AC[gnomad_index] > 0) &
                                   (reg_gnomad.freq.AN[gnomad_index] > (0.8 * max_AN)), keep = True)
transcriptIndex = reg_gnomad_flt.vep.transcript_consequences.transcript_id.index(reg_gnomad_flt.transcript)
        
# Synonymous
reg_gnomad_syn = reg_gnomad_flt.filter(
    reg_gnomad_flt.vep.transcript_consequences.consequence_terms[transcriptIndex][0] == "synonymous_variant", keep = True)
reg_gnomad_syn = reg_gnomad_syn.drop('vep', 'freq')
reg_gnomad_syn_count = reg_gnomad_syn.group_by(reg_gnomad_syn.interval).aggregate(n_syn=hl.agg.count())
reg_gnomad_syn = reg_gnomad_syn.key_by('interval')
reg_gnomad_syn = reg_gnomad_syn.annotate(n_syn = reg_gnomad_syn_count[reg_gnomad_syn.interval].n_syn)
reg_gnomad_syn = reg_gnomad_syn.distinct()
# Missense
reg_gnomad_mis = reg_gnomad_flt.filter(
    reg_gnomad_flt.vep.transcript_consequences.consequence_terms[transcriptIndex][0] == "missense_variant", keep = True)
reg_gnomad_mis = reg_gnomad_mis.drop('vep', 'freq')
reg_gnomad_mis_count = reg_gnomad_mis.group_by(reg_gnomad_mis.interval).aggregate(n_mis=hl.agg.count())
reg_gnomad_mis = reg_gnomad_mis.key_by('interval')
reg_gnomad_mis = reg_gnomad_mis.annotate(n_mis = reg_gnomad_mis_count[reg_gnomad_mis.interval].n_mis)
reg_gnomad_mis = reg_gnomad_mis.distinct()
# Protein Truncating Variants
reg_gnomad_ptv = reg_gnomad_flt.filter(
    prot_truncating.contains(reg_gnomad_flt.vep.transcript_consequences.consequence_terms[transcriptIndex][0]), keep = True)
reg_gnomad_ptv = reg_gnomad_ptv.drop('vep', 'freq')
reg_gnomad_ptv_count = reg_gnomad_ptv.group_by(reg_gnomad_ptv.interval).aggregate(n_ptv=hl.agg.count())
reg_gnomad_ptv = reg_gnomad_ptv.key_by('interval')
reg_gnomad_ptv = reg_gnomad_ptv.annotate(n_ptv = reg_gnomad_ptv_count[reg_gnomad_ptv.interval].n_ptv)
reg_gnomad_ptv = reg_gnomad_ptv.distinct()
        
reg_gnomad_exportTMP = reg_gnomad_syn.annotate(n_mis = reg_gnomad_mis[reg_gnomad_syn.interval].n_mis)
reg_gnomad_export = reg_gnomad_exportTMP.annotate(n_ptv = reg_gnomad_ptv[reg_gnomad_exportTMP.interval].n_ptv)
reg_gnomad_export = reg_gnomad_export.drop('locus')
reg_gnomad_export = reg_gnomad_export.key_by()
reg_gnomad_export = reg_gnomad_export.annotate(interval = hl.str(reg_gnomad_export.interval))
subregfile = reg_gnomad_export.to_pandas()

transD = {} # transcript -> subRegion -> exon -> chr, start, end, obs, proba, etc
# For each requested regions
with open(output, 'w') as out:
    for line in subregfile.iterrows():
        idx = line[1]
        chrom = int(idx['interval'].split(':')[0].split('[')[1])
        gStart = int(idx['interval'].split(':')[1].split('-')[0])
        gEnd = int(idx['interval'].split(':')[2].split(')')[0])
        transcript = idx['transcript']
        exon = idx['exon']
        subRegionReg = idx['subRegionReg'][0]
        regionLength = idx['regionLength']
        geneID = idx['geneID']
        geneSymb = idx['geneSymb']
        obsN_syn = idx['n_syn']
        obsN_mis = idx['n_mis']
        obsN_ptv = idx['n_ptv']
        # Prediction on the given Region
        if transcript not in transD:
            transD[transcript] = {subRegionReg :{
                'geneID': geneID,
                'geneSymb': geneSymb,
                'reg':str(gStart)+'-'+str(gEnd),
                'chr':chrom,
                'exon': str(exon),
                'regionLength':regionLength,
                'obsN_syn':obsN_syn,
                'obsN_mis':obsN_mis,
                'obsN_ptv':obsN_ptv}}
        elif subRegionReg not in transD[transcript]:
            transD[transcript][subRegionReg] = {
                'geneID': geneID,
                'geneSymb': geneSymb,
                'reg':str(gStart)+'-'+str(gEnd),
                'chr':chrom,
                'exon': str(exon),
                'regionLength':regionLength,
                'obsN_syn':obsN_syn,
                'obsN_mis':obsN_mis,
                'obsN_ptv':obsN_ptv}
        else:
            transD[transcript][subRegionReg]['reg'] += '|' + str(gStart)+'-'+str(gEnd)
            transD[transcript][subRegionReg]['exon'] += '|' + str(exon)
            transD[transcript][subRegionReg]['regionLength'] += regionLength
            transD[transcript][subRegionReg]['obsN_syn'] += obsN_syn
            transD[transcript][subRegionReg]['obsN_mis'] += obsN_mis
            transD[transcript][subRegionReg]['obsN_ptv'] += obsN_ptv
    # write the output
    out.write("transcript\tgeneID\tgeneSymbol\tchr\tRegionID\tInputOrgRegions\tMappedGencodeGenomicRegions\tRegion_Length\tRegion_MappedExons\tRegion_syn_obs\tRegion_mis_obs\tRegion_ptv_obs\n")
    print("Writing output...")
    for local_transcript in transD:
        # total prediction of this transcript
        for subR in transD[local_transcript]:
            subRegion = transD[local_transcript][subR]
            RegID = subR.split('|')[0]
            RegStart = subR.split('|')[1]
            RegEnd = subR.split('|')[2]
            out.write("%s\n" % '\t'.join(map(str, [
                local_transcript,
                subRegion['geneID'],
                subRegion['geneSymb'],
                subRegion['chr'],
                RegID,
                RegStart + '-' + RegEnd,
                subRegion['reg'],
                subRegion['regionLength'],
                subRegion['exon'],
                subRegion['obsN_syn'],
                subRegion['obsN_mis'],
                subRegion['obsN_ptv']])))
    
