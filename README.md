# Case-only Burden Test (COBT)

COBT is a rare variant collapsing method that does not require controls. The method is a parametric test that uses [gnomAD](https://gnomad.broadinstitute.org/) exome data v2.1.1.\
The method can be applied on WES and panel sequencing data aligned on GrCH37/hg19.

Requirements:
* conda (anaconda or miniconda, tested on miniconda v4.9.2)
* R (tested on v4.1.)
You can use your own version of R, but we recommand using the conda Renv provided with this software.

First you need to clone this repo:
```
git clone https://github.com/cbl-imagine/Case-Only-Burden-Test.git
```
And then you will need to download and process several files, then you will be able to run the software on your workstation (server with multiple cores is recommended).

## DOWNLOADS AND PREPARATION
First some files need to be downloaded and processed to run the tool.

### Create the conda environments
First, change directory to the folder where you installed COBT:
```
cd COBT
```
Then, create the main environement, the specific environment for VEP annotation and the R environment:
```
conda env create -f case_only_burden_test.yml
conda env create -f vep_env.yml
conda env create -f Renv.yml
```
And create the specific environment for VEP annotation:
```
conda env create -f vep_env.yml
```
You also need to set up the VEP cache file.
**Replace `.cache_dir` with the path to your cache directory (for example `/home/my_work/.vep_cache/`)**:
```
vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh37 --CACHEDIR .path_to_cache_dir
```

### Downloads
Download gnomAD exomes files (**this is very long**):
```
wget -P data/ https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz
conda activate COBT
gsutil cp -r gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht data/
conda deactivate
```
Download Gencode v19 file:
```
wget -P data/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
```
Download APPRIS database for GRCh37 assembly and Gencode v19 (and creation of a file containing only APPRIS PRINCIPAL transcripts):
```
wget -P data/ http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh37/appris_data.principal.txt
```
Download ENCODE blacklist from [Boyle lab](https://doi.org/10.1038/s41598-019-45839-z):
```
wget -P data/ https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
```
### Additional explanations for some files

The exclusion table `data/Table_S7_gene_exclusion_list_final.txt` comes from [Fuentes Fajardo et al. (2012)](https://doi.org/10.1002/humu.22033)'s paper.
The protein coding information table `data/proteincodinggene_db_biomart_hg19_grch37.tsv` comes from [BioMart](http://www.ensembl.org/info/data/biomart/index.html).

### Prepare a coverage file
In order to account for regions supposedly sufficiently covered in depth (we recommend a coverage depth of 10X in at least 90% of the samples), you need to provide a 1-based bed-like file (tabulated separated columns file with 4 columns: chromosome\t start\t end\t number of samples with at least 10X of coverage depth), the file should preferably to be bgzipped.\
**Chromosome in the first column should only be a digit, no `chr` character before**

## RUN THE PROGRAM

Run the script that will prepare all the files you need to perform the burden testing. You have to specify several parameters:
* the path to your conda (anaconda or miniconda) installation, more precisevly to your conda.sh file in order to run conda in the shell script `-p`
* the number of cores you want to use for parallelization `-c`
* the maximum memory to allocate for hail scripts in gigabytes `-m`
* the path to your variant file (bgzipped vcf) `-v`
* the path to your coverage file (tsv) `-f`
* the path to your hail temporary folder `-t`
* the path to your VEP cache folder `-d`

Here is an example:
```
bash src/PIPELINE.sh -p /usr/local/miniconda3/etc/profile.d/conda.sh -c 16 -m 100g -v /data-cbl/Ciliome/ciliomes_Nov2020.vcf.gz -f /data-cbl/Ciliome/10.count.497ciliomes.tsv.gz -t /data-tmp/antoine_data/hailTMP -d /home/my_work/.vep_cache/
```
Once you have all the needed files, **you can run the Case-only Burden test**, choosing if you want to run it *only on genes or on pathways* as well. You can choose between:
* [Hallmark](https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#H)
* [Reactome](https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C2)
* [KEGG](https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C2)

Pathways must be comma-separated. Here are a few examples:
```
conda activate Renv
Rscript src/case_only_burden_test.R
Rscript src/case_only_burden_test.R -p kegg
Rscript src/case_only_burden_test.R -p kegg,hallmark,reactome
```
Histograms of the distributions of p-values and manhattan plots of the p-values will be stored in `plots/` folder and lists of genes significantly enriched in synonymous, missense, protein truncating variants and protein altering variants will be stored in the `results/` folder.

If you want to have access to the variants that were used to detect the hit genes, run this additional script still specifying:
* the maximum memory to allocate for hail scripts in gigabytes `-m`
* the path to your hail temporary folder `-t`
As in the following example:
```
conda activate COBT
python src/variants_genes_table.py -m 100g -t /data-tmp/antoine_data/hailTMP
```
The tables gathering the variants found in the burden genes identfied by the test will be stored in the `results/` folder.
