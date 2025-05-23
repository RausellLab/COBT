#!/bin/bash
# antoine.favier@institutimagine.org

help()
{
  echo "Usage:
    Path to conda installation (path the the conda.sh file) [ -p | --conda_path ]
    Number of cores you want to use for parallelization [ -c | --cores ]
    Maximum memory to allocate per core [ -m | --max_memory ]
    Path to the bgzipped vcf of the case samples [ -v | --variants_samples ]
    Path to the coverage file of the case samples [ -f | --cov_samples ]
    Path to a folder to store hail temporary files [ -t | --tmp_hail ]
    Path to the VEP cache directory [ -d | --vep_cache_dir ]"
  exit 2
}

SHORT=p:,c:,m:,v:,f:,t:,,d:,h
LONG=conda_path:,cores:,max_memory:,variants_samples:,cov_samples:,tmp_hail:,vep_cache_dir:,help
OPTS=$(getopt -a --options $SHORT --longoptions $LONG -- "$@")

if [ "$#" -eq 0 ]; then
        help
fi

eval set -- "$OPTS"

while :
do
  case "$1" in
    -p | --conda_path)
      conda_path="$2"
      shift 2
      ;;
    -c | --cores )
      cores="$2"
      shift 2
      ;;
    -m | --max_memory )
      max_memory="$2"
      shift 2
      ;;
    -v | --variants_samples )
      variants_samples="$2"
      shift 2
      ;;
    -f | --cov_samples )
      cov_samples="$2"
      shift 2
      ;;
    -t | --tmp_hail )
      tmp_hail="$2"
      shift 2
      ;;
    -d | --vep_cache_dir)
      vep_cache_dir="$2"
      shift 2
      ;;
    -h | --help )
      help
      ;;
    -- )
      shift;
      break
      ;;
    * )
      echo "Unexpected option: $1"
      help
      ;;
  esac
done

source $conda_path
conda activate COBT

mkdir intermediate_files
mkdir intermediate_files/TEMP

##### CREATE REGION FILES #####

# Create 1-based tsv file (pseudo-bed file) with all regions covered in at least 10X in 90% of patients
PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell" python src/Filter_patients_coverage_file.py $cores $max_memory $tmp_hail $cov_samples
# Filter gnomAD coverage file to create a 1-based tsv file (pseudo-bed file) of regions covered at 10X in at least 90% of the samples
PYSPARK_SUBMIT_ARGS="--driver-memory 120g pyspark-shell" python src/Filter_gnomAD_coverage_file.py $cores $max_memory $tmp_hail
# Intersect 1000G 1-based tsv mask with gnomAD
intersectBed -a intermediate_files/filtered_coverage_gnomad.tsv -b intermediate_files/filtered_coverage_patients.tsv > intermediate_files/gnomAD_patients_1_based.tsv
# Add some columns to the pseudo-bedfile and create more tmp files
python src/reformat.py intermediate_files/gnomAD_patients_1_based.tsv intermediate_files/TEMP/gnomAD_patients_1_based_reg.tsv
awk '{print $0"\t"$2"\t"$3}' intermediate_files/TEMP/gnomAD_patients_1_based_reg.tsv > intermediate_files/TEMP/gnomAD_patients_1_based_reg_more_cols.tsv
# Create APPRIS principal transcripts file
grep 'PRINCIPAL' data/appris_data.principal.txt > intermediate_files/appris_principal.txt
gunzip -c /data-cbl/gencode/gencode.v19.annotation.gtf.gz > intermediate_files/gencode.gtf
# Create the tsv file (1-based bed file) with all APPRIS Principal transcripts
python src/restriction.py intermediate_files/appris_principal.txt intermediate_files/gencode.gtf data/gencode.CDS.appris_all.tsv all
# Create the gencode file with only ONE APPRIS Principal transcript per gene
python src/restriction.py intermediate_files/appris_principal.txt intermediate_files/gencode.gtf intermediate_files/gencode.CDS.appris.tsv one
sort -V -k 1 -k 2 intermediate_files/gencode.CDS.appris.tsv > data/gencode.CDS.appris.sort.tsv
# Intersect gnomAD well-covered regions with patients well-covered regions
intersectBed -a intermediate_files/TEMP/gnomAD_patients_1_based_reg_more_cols.tsv -b data/gencode.CDS.appris.sort.tsv -wb | awk '{print $1"\t"$2"\t"$3"\t"$11"\t"$13"\t"$4"\t"$6"\t"$7"\t"$12}' | sort -V -k 1 -k 2 > intermediate_files/filtered_coverage_gnomad_patients_1_based_gencode.tsv
# Create a 0-based bedfile with the same regions
awk -v n=1 '{print $1, $2-n, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' intermediate_files/filtered_coverage_gnomad_patients_1_based_gencode.tsv > intermediate_files/filtered_coverage_gnomad_patients_0_based_gencode.bed
# Remove duplicated regions in the bedfile
awk '!a[$1$2$3]++' intermediate_files/filtered_coverage_gnomad_patients_0_based_gencode.bed > intermediate_files/filtered_coverage_gnomad_patients_0_based_gencode_unduplicated.bed
# Remove segmental duplications and create the bed file to extract regions from VCFs
subtractBed -a intermediate_files/filtered_coverage_gnomad_patients_0_based_gencode_unduplicated.bed -b <(zcat data/hg19-blacklist.v2.bed.gz | cut -f1,2,3 | sed 's/chr//') > intermediate_files/gnomad_patients_all_chr.bed
awk '{print $1"\t"$2"\t"$3"\t"$4":"$5":"$6":"$7"-"$8":"$9}' intermediate_files/gnomad_patients_all_chr.bed > data/gnomad_patients_all_chr.bed

##### ANNOTATE VCFs #####

# First convert from multi to bi-allelic
bcftools norm -Oz -m -any $variants_samples -o data/patients_biallelic.vcf.bgz

conda activate vep_env
vep -i data/patients_biallelic.vcf.bgz -o data/patients_annotated.vcf.bgz -e --species homo_sapiens -a GRCh37 --format vcf --no_stats --fork $cores --cache --vcf --dir_cache $vep_cache_dir --offline --compress_output bgzip > intermediate_files/TEMP/annotation_error.log

##### RUN THE PREDICTIONS OF NUMBER OF VARIANTS IN GNOMAD #####

conda activate COBT
python src/gnomAD_count_proba.py -m $max_memory -t $tmp_hail -i data/gnomad_patients_all_chr.bed -o data/gnomad_all_chr_predicted.tsv
# Create a file with one count of variants per gene in gnomAD
python src/per_gene_count_gnomAD.py

##### FILTERING AND CREATION OF ANALYSIS TABLES #####

# Create a list of the transcripts to keep
cut -f2 data/gnomAD_count_table.tsv | sort | uniq -c | grep -v "transcript" > data/transcript_list.tsv
# Keep only one APPRIS PRINCIPAL transcript per variant, filter by allele frequencies of the cohort and gnomAD and export all needed tables
python src/filtering_export_analysis_tables.py $max_memory $tmp_hail

##### ANALYZING THE RESULTS #####

# Download HGNC table and other info
python src/hgnc_web_scraping.py

mkdir plots
mkdir results

