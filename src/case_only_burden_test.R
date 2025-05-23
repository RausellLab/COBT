#!/usr/bin/env Rscript
# antoine.favier@institutimagine.org
# Test if some genes are over-mutated for a specific individual

rm(list=ls())

library(optparse)
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(gridExtra)
library(grid)
library(tibble)
library(ggrepel)
library(purrr)

option_list <- list(make_option(c("-p", "--pathways"), default="no",
                                help="comma separated list of pathways for pathway analysis between (without spaces)
                                hallmark
                                reactome
                                kegg
                                no if no pathway analysis is recquired
                                (default %default)"))

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

pathwaylist <- unlist(strsplit(opt$pathways, ","))

if(length(pathwaylist) > 3){
  stop("you cannot have more than 3 pathway databases")
}else if(!all(pathwaylist %in% c("hallmark", "reactome", "kegg"))){
  stop("you probably made a typing error in the pathways")
}

theme_hist <- theme(panel.background=element_rect(fill="white"),
                    panel.grid=element_line(colour="gray"),
                    plot.title=element_text(hjust=0.5, size=18),
                    panel.grid.minor=element_blank(),
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=14))

# Load blacklist of overmutated genes
blacklist_gene_symb <- as.character(scan(file="data/Table_S7_gene_exclusion_list_final.txt", what="", sep="\n"))

# Theoretical probabilities of mutation
prob_theor <- read_tsv("data/gnomAD_count_table.tsv") %>%
  dplyr::mutate(chromosome=as.character(chromosome)) %>%
  dplyr::mutate(transcript=as.character(transcript)) %>%
  dplyr::mutate(geneID=as.character(geneID)) %>%
  dplyr::mutate(obs_alt=obs_mis+obs_ptv)

# HGNC protein-coding table
hgnc <- read_delim("data/hgnc_protein_coding.txt", col_names=T, delim="\t") %>%
  dplyr::select(symbol, ensembl_gene_id) %>% dplyr::rename(geneSymbol=symbol, ensg=ensembl_gene_id)

# Table with ENSG and gene symbols correspondance and other stuff
genedb <- read_tsv("data/proteincodinggene_db_biomart_hg19_grch37.tsv") %>% dplyr::select(-c("Gene start (bp)", "Gene end (bp)")) %>%
  dplyr::rename(ENSG="Gene stable ID", transcript="Transcript stable ID", chrom="Chromosome/scaffold name",
                TransStart="Transcript start (bp)", TransEnd="Transcript end (bp)", TSS="Transcription start site (TSS)",
                TranscriptLength="Transcript length (including UTRs and CDS)", geneSymbol="Gene name", GC_perc="Gene % GC content") %>%
  filter(chrom %in% prob_theor$chromosome) %>% filter(transcript %in% prob_theor$transcript) %>%
  dplyr::mutate(usedtrans=map_dbl(ENSG, function(x) prob_theor$length[which(prob_theor$geneID == x)])) %>%
  dplyr::mutate(ratiotrans=(usedtrans/TranscriptLength)*100)
# Poorly covered genes
poorlyCovGenes <- subset(genedb, ratiotrans <= 25)$ENSG
# Remove too poorly genes
genedb <- subset(genedb, ratiotrans > 25)

# Remove blacklisted and poorly covered genes from theoretical probabilities
blacklist <- hgnc %>% filter(geneSymbol %in% blacklist_gene_symb) %>% dplyr::select(ensg) %>% pull
prob_theor <- prob_theor %>% filter(!(geneID %in% blacklist)) %>% filter(!(geneID %in% poorlyCovGenes))


# MAIN TABLES -------------------------------------------------------------


# Summary of the variant numbers of the datasets
# Number of variants per gene
genes_patients <- read_tsv("data/variant_number.tsv") %>%
  dplyr::mutate(Gene_ID=as.character(Gene_ID)) %>%
  dplyr::mutate(Chromosome=as.character(Chromosome)) %>%
  dplyr::mutate(Protein_altering_variants=Missense_variants+PTV)
genes_patients <- subset(genes_patients, !(Gene_ID %in% poorlyCovGenes))

no_proba <- subset(genes_patients, !Gene_ID %in% prob_theor$geneID)$Gene_ID
# Add genes for which a probability was computed with 0 variants to the genes_patients table
zero_mut <- subset(prob_theor, !geneID %in% genes_patients$Gene_ID)$geneID
zero_mutTab <- data.frame(matrix(data=0, nrow=length(zero_mut), ncol=length(colnames(genes_patients))))
colnames(zero_mutTab) <- colnames(genes_patients)
zero_mutTab$Gene_ID <- zero_mut
zero_mutTab$Chromosome <- prob_theor$chromosome[which(prob_theor$geneID %in% zero_mut)]
genedbTMP <- genedb[match(unique(genedb$ENSG), genedb$ENSG),]
zero_mutTab$Gene_symbol <- as.character(sapply(zero_mut, function(x)
  ifelse(x %in% genedbTMP$ENSG,
         as.character(genedbTMP$geneSymbol[which(genedbTMP$ENSG == x)]),
         NA)))
rm(genedbTMP)
genes_patients <- rbind(genes_patients, zero_mutTab)

genedb_patients <- genedb %>% full_join(genes_patients, by=c("geneSymbol"="Gene_symbol")) %>%
  dplyr::select(chrom, geneSymbol, ENSG, transcript, TranscriptLength, TransStart, ratiotrans, usedtrans,
                Synonymous_variants, Missense_variants, PTV, Protein_altering_variants) %>%
  rename("transcriptLength"="TranscriptLength", "TSS"="TransStart", "ratioTrans"="ratiotrans",
         "usedTrans"="usedtrans", "syn_var"="Synonymous_variants", "mis_var"="Missense_variants",
         "prot_alt_var"="Protein_altering_variants")

# Synonymous
tot_syn <- 2*as.matrix(read.table("data/Homo_syn_var.tsv", header=T, sep="\t", check.names=F, row.names=1)) +
  as.matrix(read.table("data/Hetero_syn_var.tsv", header=T, sep="\t", check.names=F, row.names=1))
tot_syn <- subset(tot_syn, !(rownames(tot_syn) %in% no_proba)) # remove genes with no probability computed
tot_syn <- subset(tot_syn, !(rownames(tot_syn) %in% blacklist)) # remove blacklisted genes
tot_syn <- subset(tot_syn, !(rownames(tot_syn) %in% poorlyCovGenes)) # remove poorly covered genes
ENSG_non_mut <- subset(prob_theor, !geneID %in% rownames(tot_syn))$geneID
mat_non_mut <- matrix(data=0, nrow=length(ENSG_non_mut), ncol=length(colnames(tot_syn)))
rownames(mat_non_mut) <- ENSG_non_mut
tot_syn <- rbind(tot_syn, mat_non_mut)
# Missense
tot_mis <- 2*as.matrix(read.table("data/Homo_mis_var.tsv", header=T, sep="\t", check.names=F, row.names=1)) +
  as.matrix(read.table("data/Hetero_mis_var.tsv", header=T, sep="\t", check.names=F, row.names=1))
tot_mis <- subset(tot_mis, !(rownames(tot_mis) %in% no_proba)) # remove genes with no probability computed
tot_mis <- subset(tot_mis, !(rownames(tot_mis) %in% blacklist)) # remove blacklisted genes
tot_mis <- subset(tot_mis, !(rownames(tot_mis) %in% poorlyCovGenes)) # remove poorly covered genes
ENSG_non_mut <- subset(prob_theor, !geneID %in% rownames(tot_mis))$geneID
mat_non_mut <- matrix(data=0, nrow=length(ENSG_non_mut), ncol=length(colnames(tot_mis)))
rownames(mat_non_mut) <- ENSG_non_mut
tot_mis <- rbind(tot_mis, mat_non_mut)
# Protein Truncating Variants
tot_ptv <- 2*as.matrix(read.table("data/Homo_ptv_var.tsv", header=T, sep="\t", check.names=F, row.names=1)) +
  as.matrix(read.table("data/Hetero_ptv_var.tsv", header=T, sep="\t", check.names=F, row.names=1))
tot_ptv <- subset(tot_ptv, !(rownames(tot_ptv) %in% no_proba)) # remove genes with no probability computed
tot_ptv <- subset(tot_ptv, !(rownames(tot_ptv) %in% blacklist)) # remove blacklisted genes
tot_ptv <- subset(tot_ptv, !(rownames(tot_ptv) %in% poorlyCovGenes)) # remove poorly covered genes
ENSG_non_mut <- subset(prob_theor, !geneID %in% rownames(tot_ptv))$geneID
mat_non_mut <- matrix(data=0, nrow=length(ENSG_non_mut), ncol=length(colnames(tot_ptv)))
rownames(mat_non_mut) <- ENSG_non_mut
tot_ptv <- rbind(tot_ptv, mat_non_mut)
# Protein Truncating Variants added to Missense (Protein Altering)
tot_alt <- tot_ptv+tot_mis[rownames(tot_ptv),]
  
all_mut <- list("tot_syn"=tot_syn, "tot_mis"=tot_mis, "tot_ptv"=tot_ptv, "tot_alt"=tot_alt, "genedb_patients"=genedb_patients)


# CASE-ONLY BURDEN TEST ---------------------------------------------------


case_only_mut <- function(tot_mut, mut){
  nind <- dim(tot_mut)[2]
  med_mut <- median(colSums(tot_mut))
  var_ind_mut <- sapply(colnames(tot_mut), function(x) sum(tot_mut[,x]))
  tot_mut_TEMP <- tot_mut
  tot_mut <- t(sapply(rownames(tot_mut_TEMP), function(x) (tot_mut_TEMP[x,]*med_mut)/var_ind_mut))
  tot_mut[is.na(tot_mut)] <- 0 # some NAs appear because of non-muted patients
  tot_mut_sum <- round(rowSums(tot_mut))
  
  m_i_mut=nind*med_mut # n times the median of the number of mutations per individual
  p_values_mut <- rep(0, length(tot_mut_sum))
  names(p_values_mut) <- names(tot_mut_sum)
  gnomad_tot_mut <- sum(prob_theor[[paste0("obs_", mut, sep="")]])
  p_values_mut <- sapply(names(tot_mut_sum),
                         function(x) ifelse(((prob_theor[[paste0("obs_", mut, sep="")]][which(prob_theor$geneID == x)] == 0) && (tot_mut_sum[x] == 0)),
                                            NA,
                                            ppois(tot_mut_sum[x],
                                                  m_i_mut*(prob_theor[[paste0("obs_", mut, sep="")]][which(prob_theor$geneID == x)]/gnomad_tot_mut),
                                                  lower.tail=F)))
  p_values_mut_adj <- p.adjust(p_values_mut, method="bonferroni")
  chisq <- qchisq(p_values_mut, 1, lower.tail=FALSE)
  lambda <- median(na.omit(chisq))/qchisq(0.5, 1)
  newchisq <- chisq/lambda
  p_lambda_corrected <- pchisq(newchisq, df=1, lower.tail=FALSE)
  p_lambda_corrected_adj <- p.adjust(p_lambda_corrected, method="bonferroni")
  write.table(x=tot_mut_TEMP, file=paste("data/", mut, "_table.txt", sep=""),
              quote=F, sep="\t", row.names=T, col.names=T)
  returnT <- tibble("ensg"=names(p_values_mut), "pval"=p_values_mut) %>%
    inner_join(tibble("ensg"=names(p_values_mut_adj), "pval_adj"=p_values_mut_adj)) %>%
    inner_join(tibble("ensg"=names(p_lambda_corrected), "pval_lambda_cor"=p_lambda_corrected)) %>%
    inner_join(tibble("ensg"=names(p_lambda_corrected_adj), "pval_lambda_cor_adj"=p_lambda_corrected_adj))
  return(returnT)
}
  
# Synonymous
synT <- case_only_mut(all_mut$tot_syn, "syn") %>% rename_with(.fn=~gsub("pval", "pvalsyn", .), .cols=contains("pval"))

# Missense
misT <- case_only_mut(all_mut$tot_mis, "mis") %>% rename_with(.fn=~gsub("pval", "pvalmis", .), .cols=contains("pval"))

# PTVs
ptvT <- case_only_mut(all_mut$tot_ptv, "ptv") %>% rename_with(.fn=~gsub("pval", "pvalptv", .), .cols=contains("pval"))

# Protein Altering
altT <- case_only_mut(all_mut$tot_alt, "alt") %>% rename_with(.fn=~gsub("pval", "pvalalt", .), .cols=contains("pval"))

mutT <- synT %>% inner_join(misT, by="ensg") %>% inner_join(ptvT, by="ensg") %>% inner_join(altT, by="ensg")

# For the whole dataset
data <- list(tot_syn=all_mut$tot_syn, tot_mis=all_mut$tot_mis, tot_ptv=all_mut$tot_ptv,
             tot_alt=all_mut$tot_alt, genedb_patients=all_mut$genedb_patients, case_only=mutT)


# PLOT CASE-ONLY BURDEN TEST DISTRIBUTIONS --------------------------------


synHist <- ggplot(data$case_only, aes_string(x="pvalsyn_lambda_cor")) + geom_histogram(color="black", fill="white", bins=20) +
  theme_hist + ggtitle("Synonymous case-only burden test") + xlab("p-values")
misHist <- ggplot(data$case_only, aes_string(x="pvalmis_lambda_cor")) + geom_histogram(color="black", fill="white", bins=20) +
  theme_hist + ggtitle("Missense case-only burden test") + xlab("p-values")
ptvHist <- ggplot(data$case_only, aes_string(x="pvalptv_lambda_cor")) + geom_histogram(color="black", fill="white", bins=20) +
  theme_hist + ggtitle("PTV case-only burden test") + xlab("p-values")
altHist <- ggplot(data$case_only, aes_string(x="pvalalt_lambda_cor")) + geom_histogram(color="black", fill="white", bins=20) +
  theme_hist + ggtitle("Protein Altering variants case-only burden test") + xlab("p-values")

pdf("plots/Histograms_case_only_burden_test.pdf", width=14, height=8)
grid.arrange(synHist, misHist, ptvHist, altHist, ncol=2, nrow=2)
dev.off()

manPlotWholeCohort <- function(pvalmut, main){
  # Create the table data table for the plot
  dat <- data$genedb_patients %>% dplyr::inner_join(tibble("ENSG"=names(data$case_only[[paste(pvalmut)]]),
                                                              "P"=data$case_only[[paste(pvalmut)]],
                                                              "P_adj"=data$case_only[[paste0(pvalmut, "_adj")]]),
                                                            by="ENSG") %>%
    dplyr::select(geneSymbol, chrom, TSS, P, P_adj) %>% dplyr::rename(CHR=chrom, SNP=geneSymbol, BP=TSS) %>%
    dplyr::mutate(CHR=as.integer(CHR))
  # Compute the cumulative distribution for the plot
  datPlot <- dat %>%
    # Compute chromosome size
    dplyr::group_by(CHR) %>% dplyr::summarise(chr_len=max(BP)) %>%
    # Calculate cumulative position of each chromosome
    dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    dplyr::left_join(dat, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    dplyr::arrange(CHR, BP) %>% dplyr::mutate(BPcum=BP+tot) %>%
    # Add annotation information (is or is not to annotate)
    dplyr::mutate(is_annotate=ifelse(P_adj <= 0.05, "yes", "no"))
  # Add chromosome names to the cumulative positions
  xAxisDat <- datPlot %>% dplyr::group_by(CHR) %>%
    dplyr::summarize(center=(max(BPcum) + min(BPcum))/2)
  # Create the plot itself
  ggplot(datPlot, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.6, size=1.3) + scale_color_manual(values=rep(c("grey", "darkblue"), 22)) +
    # Customize X axis
    scale_x_continuous(label=xAxisDat$CHR, breaks=xAxisDat$center) +
    # Remove space between plot area and x axis
    scale_y_continuous(expand=c(0, 0)) +
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=subset(datPlot, is_annotate=="yes"), aes(label=SNP), size=2) +
    # Add the line at Bonferroni correction level
    geom_hline(yintercept=-log10(0.05/nrow(data$case_only)), color="darkred") +
    # Custom the theme and add axis names and the title
    labs(title=main, x="Chromosome", y=expression(-log[10](italic(p-value)))) +
    theme_bw() +
    theme(legend.position="none", panel.border=element_blank(), plot.title=element_text(hjust=0.5),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
}

pdf("plots/Manhattan_plots_whole_cohort.pdf", width=14, height=8)
manPlotWholeCohort("pvalsyn_lambda_cor", "Synonymous")
manPlotWholeCohort("pvalmis_lambda_cor", "Missense")
manPlotWholeCohort("pvalptv_lambda_cor", "Protein Truncating Variants")
manPlotWholeCohort("pvalalt_lambda_cor", "Protein Altering Variants")
dev.off()


# NUMBER OF HITS FOR WHOLE COHORT -----------------------------------------


data$exportTable <- data$genedb_patients %>% inner_join(data$case_only, by=c("ENSG"="ensg")) %>%
  left_join(prob_theor, by=c("ENSG"="geneID")) %>% dplyr::select(-any_of(c("transcript.y", "chromosome"))) %>%
  left_join(tibble("ENSG"=rownames(data$tot_syn),
                   "n_ind_mut_syn"=sapply(rownames(data$tot_syn),
                                          function(x) length(which(data$tot_syn[x,] != 0)))), by="ENSG") %>%
  left_join(tibble("ENSG"=rownames(data$tot_syn), "sum_accros_ind_syn"=rowSums(data$tot_syn)), by="ENSG") %>%
  left_join(tibble("ENSG"=rownames(data$tot_mis),
                   "n_ind_mut_mis"=sapply(rownames(data$tot_mis),
                                          function(x) length(which(data$tot_mis[x,] != 0)))), by="ENSG") %>%
  left_join(tibble("ENSG"=rownames(data$tot_mis), "sum_accros_ind_mis"=rowSums(data$tot_mis)), by="ENSG") %>%
  left_join(tibble("ENSG"=rownames(data$tot_ptv),
                   "n_ind_mut_PTV"=sapply(rownames(data$tot_ptv),
                                          function(x) length(which(data$tot_ptv[x,] != 0)))), by="ENSG") %>%
  left_join(tibble("ENSG"=rownames(data$tot_ptv), "sum_accros_ind_PTV"=rowSums(data$tot_ptv)), by="ENSG") %>%
  left_join(tibble("ENSG"=rownames(data$tot_alt),
                   "n_ind_mut_alt"=sapply(rownames(data$tot_alt),
                                          function(x) length(which(data$tot_alt[x,] != 0)))), by="ENSG") %>%
  left_join(tibble("ENSG"=rownames(data$tot_alt), "sum_accros_ind_alt"=rowSums(data$tot_alt)), by="ENSG") %>%
  dplyr::rename(transcript=transcript.x, gnomAD_transcript_length=length,
                NbrSynVar=syn_var, NbrMisVar=mis_var, NbrPTV=PTV, NbrProtAltVar=prot_alt_var,
                gnomAD_syn=obs_syn, gnomAD_mis=obs_mis, gnomAD_PTV=obs_ptv, gnomAD_alt=obs_alt) %>%
  dplyr::select(chrom, geneSymbol, ENSG, transcript, transcriptLength, usedTrans, ratioTrans, NbrSynVar, NbrMisVar,
                NbrPTV, NbrProtAltVar, n_ind_mut_syn, sum_accros_ind_syn, n_ind_mut_mis, sum_accros_ind_mis,
                n_ind_mut_PTV, sum_accros_ind_PTV, n_ind_mut_alt, sum_accros_ind_alt, gnomAD_transcript_length,
                gnomAD_syn, gnomAD_mis, gnomAD_PTV, gnomAD_alt, pvalsyn_lambda_cor, pvalmis_lambda_cor, pvalptv_lambda_cor,
                pvalalt_lambda_cor, pvalsyn_lambda_cor_adj, pvalmis_lambda_cor_adj, pvalptv_lambda_cor_adj,
                pvalalt_lambda_cor_adj)

# Synonymous

synHits <- data$exportTable %>% filter(pvalsyn_lambda_cor_adj <= 0.05) %>%
  dplyr::select(chrom, geneSymbol, ENSG, transcript, transcriptLength, usedTrans, ratioTrans, NbrSynVar,
                n_ind_mut_syn, sum_accros_ind_syn, gnomAD_syn, pvalsyn_lambda_cor, pvalsyn_lambda_cor_adj) %>%
  dplyr::mutate(patients_with_variant=map_chr(ENSG, function(x) paste(unlist(names(which(data$tot_syn[x,] != 0))), collapse=";"))) %>%
  dplyr::rename(nb_var=NbrSynVar, n_mut_ind=n_ind_mut_syn, sum_accros_ind=sum_accros_ind_syn, gnomAD_nb_var=gnomAD_syn,
                pval=pvalsyn_lambda_cor, pval_adj=pvalsyn_lambda_cor_adj)

write.table(synHits, file="results/synonymous_hits.tsv", sep="\t", row.names=F, quote=F)

# Missense

misHits <- data$exportTable %>% filter(pvalmis_lambda_cor_adj <= 0.05) %>%
  dplyr::select(chrom, geneSymbol, ENSG, transcript, transcriptLength, usedTrans, ratioTrans, NbrMisVar,
                n_ind_mut_mis, sum_accros_ind_mis, gnomAD_mis, pvalmis_lambda_cor, pvalmis_lambda_cor_adj) %>%
  dplyr::mutate(patients_with_variant=map_chr(ENSG, function(x) paste(unlist(names(which(data$tot_mis[x,] != 0))), collapse=";"))) %>%
  dplyr::rename(nb_var=NbrMisVar, n_mut_ind=n_ind_mut_mis, sum_accros_ind=sum_accros_ind_mis, gnomAD_nb_var=gnomAD_mis,
                pval=pvalmis_lambda_cor, pval_adj=pvalmis_lambda_cor_adj)

write.table(misHits, file="results/missense_hits.tsv", sep="\t", row.names=F, quote=F)

# PTVs

ptvHits <- data$exportTable %>% filter(pvalptv_lambda_cor_adj <= 0.05) %>%
  dplyr::select(chrom, geneSymbol, ENSG, transcript, transcriptLength, usedTrans, ratioTrans, NbrPTV,
                n_ind_mut_PTV, sum_accros_ind_PTV, gnomAD_PTV, pvalptv_lambda_cor, pvalptv_lambda_cor_adj) %>%
  dplyr::mutate(patients_with_variant=map_chr(ENSG, function(x) paste(unlist(names(which(data$tot_ptv[x,] != 0))), collapse=";"))) %>%
  dplyr::rename(nb_var=NbrPTV, n_mut_ind=n_ind_mut_PTV, sum_accros_ind=sum_accros_ind_PTV, gnomAD_nb_var=gnomAD_PTV,
                pval=pvalptv_lambda_cor, pval_adj=pvalptv_lambda_cor_adj)

write.table(misHits, file="results/ptv_hits.tsv", sep="\t", row.names=F, quote=F)

# Protein Altering

altHits <- data$exportTable %>% filter(pvalalt_lambda_cor_adj <= 0.05) %>%
  dplyr::select(chrom, geneSymbol, ENSG, transcript, transcriptLength, usedTrans, ratioTrans, NbrProtAltVar,
                n_ind_mut_alt, sum_accros_ind_alt, gnomAD_alt, pvalalt_lambda_cor, pvalalt_lambda_cor_adj) %>%
  dplyr::mutate(patients_with_variant=map_chr(ENSG, function(x) paste(unlist(names(which(data$tot_alt[x,] != 0))), collapse=";"))) %>%
  dplyr::rename(nb_var=NbrProtAltVar, n_mut_ind=n_ind_mut_alt, sum_accros_ind=sum_accros_ind_alt, gnomAD_nb_var=gnomAD_alt,
                pval=pvalalt_lambda_cor_adj, pval_adj=pvalalt_lambda_cor_adj)

write.table(altHits, file="results/protein_altering_hits.tsv", sep="\t", row.names=F, quote=F)


# TESTING BY PATHWAY ------------------------------------------------------

if("no" %in% pathwaylist){
  print("Done!")
}else{
  # Import the pathways
  library(msigdbr)
  
  pathways <- NULL
  if("hallmark" %in% pathwaylist){
    pathways$hallmark <- msigdbr(species="Homo sapiens", category="H")
  }
  if("reactome" %in% pathwaylist){
    pathways$reactome <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
  }
  if("kegg" %in% pathwaylist){
    pathways$kegg <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:KEGG")
  }
  
  # Create a coordination between pathways and genes
  path_db <- lapply(pathways, function(x) genedb %>% inner_join(x, by=c("ENSG"="human_ensembl_gene")) %>%
                      dplyr::select(chrom, ENSG, geneSymbol, gs_name, TranscriptLength, usedtrans, ratiotrans))
  
  # gnomAD counts of mutations per pathway
  prob_theor_path <- lapply(path_db, function(x) prob_theor %>% inner_join(x, by=c("geneID"="ENSG")) %>%
                              dplyr::select(gs_name, obs_syn, obs_mis, obs_ptv, obs_alt) %>%
                              group_by(gs_name) %>% dplyr::summarize(obs_syn=sum(obs_syn), obs_mis=sum(obs_mis),
                                                                     obs_ptv=sum(obs_ptv), obs_alt=sum(obs_alt)))
  
  # Function to capitalize the first letter of a string
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  # Function to create the count matrices by pathway using the count matrices by genes
  mat_path <- function(dat, path){
    returnList <- NULL
    for(mut in c("syn", "mis", "ptv", "alt")){
      returnList[[paste0("tot_", mut)]] <- path_db[[path]] %>% dplyr::select(ENSG, gs_name) %>%
        inner_join(rownames_to_column(as.data.frame(dat[[paste0("tot_", mut)]])), by=c("ENSG"="rowname")) %>%
        dplyr::select(-one_of("ENSG")) %>% group_by(gs_name) %>% summarise_all(sum) %>%
        column_to_rownames("gs_name") %>% as.matrix
    }
    returnList$path_db_patients <- dat$genedb_patients %>% right_join(path_db[[path]], by="ENSG") %>%
      dplyr::select(gs_name, syn_var, mis_var, PTV, prot_alt_var) %>% group_by(gs_name) %>%
      dplyr::summarize(syn_var=sum(syn_var), mis_var=sum(mis_var), PTV=sum(PTV), prot_alt_var=sum(prot_alt_var))
    return(returnList)
  }
  
  # Function to compute case-only for all pathways and all mutation types
  case_only_path <- function(dat){
    returnList <- NULL
    for(path in names(dat)){
      returnList[[path]] <- NULL
      for(mut in c("syn", "mis", "ptv", "alt")){
        tot_mut <- dat[[path]][[paste0("tot_", mut)]]
        nind <- dim(tot_mut)[2]
        med_mut <- median(colSums(tot_mut))
        var_ind_mut <- sapply(colnames(tot_mut), function(x) sum(tot_mut[,x]))
        tot_mut_TEMP <- tot_mut
        tot_mut <- t(sapply(rownames(tot_mut_TEMP), function(x) (tot_mut_TEMP[x,]*med_mut)/var_ind_mut))
        tot_mut[is.na(tot_mut)] <- 0 # some NAs appear because of non-muted patients
        tot_mut_sum <- round(rowSums(tot_mut))
        if(sum(tot_mut_sum) == 0){
          returnList[[path]][[mut]] <- setNames(list(structure(rep(NA, length(tot_mut_sum)), names=names(tot_mut_sum), class="numeric"),
                                                     structure(rep(NA, length(tot_mut_sum)), names=names(tot_mut_sum), class="numeric"),
                                                     tot_mut_TEMP),
                                                c(paste0("pval", mut), paste0("pval", mut, "_adj"), paste0("tab", firstup(mut))))
        }else{
          m_i_mut=nind*med_mut # n times the median of the number of mutations per individual
          p_values_mut <- rep(0, length(tot_mut_sum))
          names(p_values_mut) <- names(tot_mut_sum)
          gnomad_tot_mut <- sum(prob_theor_path[[path]][paste0("obs_", mut)])
          p_values_mut <- sapply(names(tot_mut_sum),
                                 function(x) ifelse(((prob_theor_path[[path]][[paste0("obs_", mut)]][which(prob_theor_path[[path]]$gs_name == x)] == 0) && (tot_mut_sum[x] == 0)),
                                                    NA,
                                                    ppois(tot_mut_sum[x],
                                                          m_i_mut*(prob_theor_path[[path]][[paste0("obs_", mut)]][which(prob_theor_path[[path]]$gs_name == x)]/gnomad_tot_mut),
                                                          lower.tail=F)))
          p_values_mut_adj <- p.adjust(p_values_mut, method="bonferroni")
          chisq <- qchisq(p_values_mut, 1, lower.tail=FALSE)
          lambda <- median(na.omit(chisq))/qchisq(0.5, 1)
          newchisq <- chisq/lambda
          p_lambda_corrected <- pchisq(newchisq, df=1, lower.tail=FALSE)
          p_lambda_corrected_adj <- p.adjust(p_lambda_corrected, method="bonferroni")
          returnList[[path]][[mut]] <- setNames(list(p_lambda_corrected, p_lambda_corrected_adj, tot_mut_TEMP),
                                                c(paste0("pval", mut), paste0("pval", mut, "_adj"), paste0("tab", firstup(mut))))
        }
      }
    }
    return(returnList)
  }
  
  # Function to group initial matrix and case-only results by pathway
  group_mat_res <- function(dat, datTMP){
    returnList <- NULL
    for(path in names(dat)){
      returnList[[path]] <- list("tot_syn"=dat[[path]]$tot_syn, "tot_mis"=dat[[path]]$tot_mis,
                                 "tot_ptv"=dat[[path]]$tot_ptv, "tot_alt"=dat[[path]]$tot_alt,
                                 "path_db_patients"=dat[[path]]$path_db_patients,
                                 "pvalsyn"=datTMP[[path]]$syn$pvalsyn, "pvalmis"=datTMP[[path]]$mis$pvalmis,
                                 "pvalptv"=datTMP[[path]]$ptv$pvalptv, "pvalalt"=datTMP[[path]]$alt$pvalalt,
                                 "pvalsyn_adj"=datTMP[[path]]$syn$pvalsyn_adj, "pvalmis_adj"=datTMP[[path]]$mis$pvalmis_adj,
                                 "pvalptv_adj"=datTMP[[path]]$ptv$pvalptv_adj, "pvalalt_adj"=datTMP[[path]]$alt$pvalalt_adj,
                                 "tabSyn"=datTMP[[path]]$syn$tabSyn, "tabMis"=datTMP[[path]]$mis$tabMis,
                                 "tabPtv"=datTMP[[path]]$ptv$tabPtv, "tabAlt"=datTMP[[path]]$alt$tabAlt)
    }
    return(returnList)
  }
  
  # Function to get a summary table with the main information from each pathway database
  exportTablePath <- function(dat, path){
    returnTable <- dat[[path]]$path_db_patients %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalsyn),
                        "pvalsyn"=dat[[path]]$pvalsyn), by="gs_name") %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalmis),
                        "pvalmis"=dat[[path]]$pvalmis), by="gs_name") %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalptv),
                        "pvalptv"=dat[[path]]$pvalptv), by="gs_name") %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalalt),
                        "pvalalt"=dat[[path]]$pvalalt), by="gs_name") %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalsyn_adj),
                        "pvalsyn_adj"=dat[[path]]$pvalsyn_adj), by="gs_name") %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalmis_adj),
                        "pvalmis_adj"=dat[[path]]$pvalmis_adj), by="gs_name") %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalptv_adj),
                        "pvalptv_adj"=dat[[path]]$pvalptv_adj), by="gs_name") %>%
      inner_join(tibble("gs_name"=names(dat[[path]]$pvalalt_adj),
                        "pvalalt_adj"=dat[[path]]$pvalalt_adj), by="gs_name") %>%
      left_join(prob_theor_path[[path]], by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabSyn),
                       "n_ind_mut_syn"=sapply(rownames(dat[[path]]$tabSyn),
                                              function(x) length(which(dat[[path]]$tabSyn[x,] != 0)))), by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabSyn),
                       "sum_accros_ind_syn"=rowSums(dat[[path]]$tabSyn)), by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabMis),
                       "n_ind_mut_mis"=sapply(rownames(dat[[path]]$tabMis),
                                              function(x) length(which(dat[[path]]$tabMis[x,] != 0)))), by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabMis),
                       "sum_accros_ind_mis"=rowSums(dat[[path]]$tabMis)), by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabPtv),
                       "n_ind_mut_PTV"=sapply(rownames(dat[[path]]$tabPtv),
                                              function(x) length(which(dat[[path]]$tabPTV[x,] != 0)))), by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabPtv),
                       "sum_accros_ind_PTV"=rowSums(dat[[path]]$tabPtv)), by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabAlt),
                       "n_ind_mut_alt"=sapply(rownames(dat[[path]]$tabAlt),
                                              function(x) length(which(dat[[path]]$tabAlt[x,] != 0)))), by="gs_name") %>%
      left_join(tibble("gs_name"=rownames(dat[[path]]$tabAlt),
                       "sum_accros_ind_alt"=rowSums(dat[[path]]$tabAlt)), by="gs_name") %>%
      dplyr::rename(NbrSynVar=syn_var, NbrMisVar=mis_var, NbrPTV=PTV, NbrProtAltVar=prot_alt_var,
                    gnomAD_syn=obs_syn, gnomAD_mis=obs_mis, gnomAD_PTV=obs_ptv, gnomAD_alt=obs_alt) %>%
      dplyr::select(gs_name, NbrSynVar, NbrMisVar, NbrPTV, NbrProtAltVar, n_ind_mut_syn, sum_accros_ind_syn,
                    n_ind_mut_mis, sum_accros_ind_mis, n_ind_mut_PTV, sum_accros_ind_PTV, n_ind_mut_alt, sum_accros_ind_alt,
                    gnomAD_syn, gnomAD_mis, gnomAD_PTV, gnomAD_alt, pvalsyn, pvalmis, pvalptv, pvalalt, pvalsyn_adj,
                    pvalmis_adj, pvalptv_adj, pvalalt_adj)
    return(returnTable)
  }
  
  # Function to get everything in a single list
  case_only_path_list <- function(dat){
    ReturnList <- lapply(names(pathways), function(x) mat_path(dat, x))
    names(ReturnList) <- names(pathways)
    # Case-only burden test
    ReturnList_tmp <- case_only_path(ReturnList)
    ReturnList <- group_mat_res(ReturnList, ReturnList_tmp)
    rm(ReturnList_tmp)
    # Create the export tables
    ReturnList$hallmark$export <- exportTablePath(ReturnList, "hallmark")
    ReturnList$reactome$export <- exportTablePath(ReturnList, "reactome")
    ReturnList$kegg$export <- exportTablePath(ReturnList, "kegg")
    return(ReturnList)
  }
  
  # Function to plot histograms
  plotHistoPath <- function(dat, multipleDat=F){
    mutationText <- list("syn"="synonymous", "mis"="missense", "ptv"="protein truncating", "alt"="protein altering")
    for(path in names(dat)){
      par(mfrow=c(2,2), oma=c(0,0,2,0))
      for(mut in names(mutationText)){
        if(!all(is.na(dat[[path]][[paste0("pval", mut)]]))){
          hist(data_path[[path]][[paste0("pval", mut)]], breaks=20, xlab="p-values",
               main=paste("Histogram of", mutationText[[mut]], "variants case-only burden test p-values"))
        }else{
          plot.new()
        }
      }
      if(multipleDat == F){
        mtext(paste(toupper(path)), outer=TRUE, cex=1.5)
      }else{
        mtext(paste(toupper(path),toupper(gsub("_path", "", gsub("_", " ", as.character(substitute(dat)))))),
              outer=TRUE, cex=1.5)
      }
    }
  }
  
  # Create the matrices
  data_path <- case_only_path_list(data)
  
  # Plot the histograms of the p-values distributions
  pdf("plots/Histograms_case_only_burden_test_whole_cohort_pathways.pdf", width=14, height=8)
  plotHistoPath(data_path)
  dev.off()
  
  # Analyzing the results
  
  # Function to create the tables to export
  pathTabExport <- function(path, mut){
    if(mut == "syn"){
      returnTab <- pathways[[path]] %>% rename("ensg"="ensembl_gene") %>%
        left_join(tibble(sums=rowSums(data$tot_syn), ensg=names(rowSums(data$tot_syn))), by="ensg") %>%
        dplyr::mutate(mutated_gene=if_else(sums == 0, "not_mutated", "mutated")) %>%
        left_join(tibble(ensg=names(data$case_only$pvalsyn_lambda_cor_adj),
                         pval=data$case_only$pvalsyn_lambda_cor_adj), by="ensg") %>%
        group_by(gs_name) %>% distinct(human_gene_symbol, .keep_all=T) %>%
        summarise(kegg_nb_genes=n(),
                  dataset_nb_genes=sum(!is.na(mutated_gene)),
                  dataset_nb_mutated_genes=sum(!is.na(mutated_gene) & mutated_gene == "mutated"),
                  dataset_nb_signif_genes=sum(!is.na(pval) & (pval <= 0.05))) %>%
        full_join(data_path$kegg$export, by="gs_name") %>% filter(pvalsyn_adj <= 0.05) %>%
        dplyr::rename(NbrVar=NbrSynVar, n_ind_mut=n_ind_mut_syn, sum_accros_ind=sum_accros_ind_syn, gnomAD_Var=gnomAD_syn,
                      pval=pvalsyn, pval_adj=pvalsyn_adj) %>%
        dplyr::select(gs_name, kegg_nb_genes, dataset_nb_genes, dataset_nb_mutated_genes, dataset_nb_signif_genes,
                      NbrVar, n_ind_mut, sum_accros_ind, gnomAD_Var, pval, pval_adj)
    }else if(mut == "mis"){
      returnTab <- pathways[[path]] %>% rename("ensg"="ensembl_gene") %>%
        left_join(tibble(sums=rowSums(data$tot_mis), ensg=names(rowSums(data$tot_mis))), by="ensg") %>%
        dplyr::mutate(mutated_gene=if_else(sums == 0, "not_mutated", "mutated")) %>%
        left_join(tibble(ensg=names(data$case_only$pvalmis_lambda_cor_adj),
                         pval=data$case_only$pvalmis_lambda_cor_adj), by="ensg") %>%
        group_by(gs_name) %>% distinct(human_gene_symbol, .keep_all=T) %>%
        summarise(kegg_nb_genes=n(),
                  dataset_nb_genes=sum(!is.na(mutated_gene)),
                  dataset_nb_mutated_genes=sum(!is.na(mutated_gene) & mutated_gene == "mutated"),
                  dataset_nb_signif_genes=sum(!is.na(pval) & (pval <= 0.05))) %>%
        full_join(data_path$kegg$export, by="gs_name") %>% filter(pvalmis_adj <= 0.05) %>%
        dplyr::rename(NbrVar=NbrMisVar, n_ind_mut=n_ind_mut_mis, sum_accros_ind=sum_accros_ind_mis, gnomAD_Var=gnomAD_mis,
                      pval=pvalmis, pval_adj=pvalmis_adj) %>%
        dplyr::select(gs_name, kegg_nb_genes, dataset_nb_genes, dataset_nb_mutated_genes, dataset_nb_signif_genes,
                      NbrVar, n_ind_mut, sum_accros_ind, gnomAD_Var, pval, pval_adj)
    }else if(mut == "ptv"){
      returnTab <- pathways[[path]] %>% rename("ensg"="ensembl_gene") %>%
        left_join(tibble(sums=rowSums(data$tot_ptv), ensg=names(rowSums(data$tot_ptv))), by="ensg") %>%
        dplyr::mutate(mutated_gene=if_else(sums == 0, "not_mutated", "mutated")) %>%
        left_join(tibble(ensg=names(data$case_only$pvalptv_lambda_cor_adj),
                         pval=data$case_only$pvalptv_lambda_cor_adj), by="ensg") %>%
        group_by(gs_name) %>% distinct(human_gene_symbol, .keep_all=T) %>%
        summarise(kegg_nb_genes=n(),
                  dataset_nb_genes=sum(!is.na(mutated_gene)),
                  dataset_nb_mutated_genes=sum(!is.na(mutated_gene) & mutated_gene == "mutated"),
                  dataset_nb_signif_genes=sum(!is.na(pval) & (pval <= 0.05))) %>%
        full_join(data_path$kegg$export, by="gs_name") %>% filter(pvalptv_adj <= 0.05) %>%
        dplyr::rename(NbrVar=NbrPTV, n_ind_mut=n_ind_mut_PTV, sum_accros_ind=sum_accros_ind_PTV, gnomAD_Var=gnomAD_PTV,
                      pval=pvalptv, pval_adj=pvalptv_adj) %>%
        dplyr::select(gs_name, kegg_nb_genes, dataset_nb_genes, dataset_nb_mutated_genes, dataset_nb_signif_genes,
                      NbrVar, n_ind_mut, sum_accros_ind, gnomAD_Var, pval, pval_adj)
    }else if(mut == "alt"){
      returnTab <- pathways[[path]] %>% rename("ensg"="ensembl_gene") %>%
        left_join(tibble(sums=rowSums(data$tot_alt), ensg=names(rowSums(data$tot_alt))), by="ensg") %>%
        dplyr::mutate(mutated_gene=if_else(sums == 0, "not_mutated", "mutated")) %>%
        left_join(tibble(ensg=names(data$case_only$pvalalt_lambda_cor_adj),
                         pval=data$case_only$pvalalt_lambda_cor_adj), by="ensg") %>%
        group_by(gs_name) %>% distinct(human_gene_symbol, .keep_all=T) %>%
        summarise(kegg_nb_genes=n(),
                  dataset_nb_genes=sum(!is.na(mutated_gene)),
                  dataset_nb_mutated_genes=sum(!is.na(mutated_gene) & mutated_gene == "mutated"),
                  dataset_nb_signif_genes=sum(!is.na(pval) & (pval <= 0.05))) %>%
        full_join(data_path$kegg$export, by="gs_name") %>% filter(pvalalt_adj <= 0.05) %>%
        dplyr::rename(NbrVar=NbrProtAltVar, n_ind_mut=n_ind_mut_alt, sum_accros_ind=sum_accros_ind_alt, gnomAD_Var=gnomAD_alt,
                      pval=pvalalt, pval_adj=pvalalt_adj) %>%
        dplyr::select(gs_name, kegg_nb_genes, dataset_nb_genes, dataset_nb_mutated_genes, dataset_nb_signif_genes,
                      NbrVar, n_ind_mut, sum_accros_ind, gnomAD_Var, pval, pval_adj)
    }
    return(returnTab)
  }
  
  if("hallmark" %in% pathwaylist){
    export_syn_hallmark <- pathTabExport("hallmark", "syn")
    export_mis_hallmark <- pathTabExport("hallmark", "mis")
    export_ptv_hallmark <- pathTabExport("hallmark", "ptv")
    export_alt_hallmark <- pathTabExport("hallmark", "alt")
    write.table(export_syn_hallmark, "results/hallmark_syn_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_mis_hallmark, "results/hallmark_mis_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_ptv_hallmark, "results/hallmark_ptv_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_alt_hallmark, "results/hallmark_alt_hits.tsv", sep="\t", quote=F, row.names=F)
  }
  
  if("reactome" %in% pathwaylist){
    export_syn_reactome <- pathTabExport("reactome", "syn")
    export_mis_reactome <- pathTabExport("reactome", "mis")
    export_ptv_reactome <- pathTabExport("reactome", "ptv")
    export_alt_reactome <- pathTabExport("reactome", "alt")
    write.table(export_syn_reactome, "results/reactome_syn_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_mis_reactome, "results/reactome_mis_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_ptv_reactome, "results/reactome_ptv_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_alt_reactome, "results/reactome_alt_hits.tsv", sep="\t", quote=F, row.names=F)
  }
  
  if("reactome" %in% pathwaylist){
    export_syn_KEGG <- pathTabExport("kegg", "syn")
    export_mis_KEGG <- pathTabExport("kegg", "mis")
    export_ptv_KEGG <- pathTabExport("kegg", "ptv")
    export_alt_KEGG <- pathTabExport("kegg", "alt")
    write.table(export_syn_KEGG, "results/KEGG_syn_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_mis_KEGG, "results/KEGG_mis_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_ptv_KEGG, "results/KEGG_ptv_hits.tsv", sep="\t", quote=F, row.names=F)
    write.table(export_alt_KEGG, "results/KEGG_alt_hits.tsv", sep="\t", quote=F, row.names=F)
  }
  
  print("Done!")
}
