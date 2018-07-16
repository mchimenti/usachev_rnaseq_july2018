## Analysis of Usachev/zhihong RNA-Seq mouse MCU knockout in hippocampus and cortex
## Date: 7.13.2018
## Author: Michael Chimenti
## Organism: mm10 / mouse 
## Aligners: hisat2 / salmon
## Design: balanced case/control in two tissues: hip and cortex
## Reps: 3
## total samples: 12

##########
## Imports
##########

#source("http://bioconductor.org/biocLite.R")

#source("http://bioconductor.org/biocLite.R")
#biocLite("COMBINE-lab/wasabi")       #install wasabi the first time

devtools::install_github("pachterlab/sleuth")  #update sleuth 

#pseudoalignment analysis

library(sleuth)
library(cowplot)

#annotation
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library(dplyr)

setwd("~/iihg/RNA_seq/usachev/project_usachev_zhihong_july2018/analysis/")

#########################
## kallisto > sleuth
#########################

base_dir <- "~/iihg/RNA_seq/usachev/project_usachev_zhihong_july2018"

kal_dirs <- file.path(base_dir, 
                      c("cox1_salmon", "cox2_salmon", "cox3_salmon",
                        "cox4_salmon", "cox5_salmon", "cox6_salmon",
                        "hipp1_salmon", "hipp2_salmon","hipp3_salmon",
                        "hipp4_salmon", "hipp5_salmon","hipp6_salmon"))

## create a R matrix containing sample names and conditions from a text file

## 'awt' = wt = wild type
## sleuth puts the factor levels in alphabetical order and chooses the first as the reference
## to keep the "WT" as the ref level we must put "a" in front of it 

s2c <- read.table(file.path(base_dir, "for_sleuth_meta.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- arrange(s2c, sname)

## add a column called "kal_dirs" containing paths to the data
s2c <- mutate(s2c, path = kal_dirs)
colnames(s2c) <- c("sample","geno","tissue","path")
## Get common gene names for transcripts

## this section queries Ensemble online database for gene names associated with transcript IDs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
t2g <- getBM(attributes = c("ensembl_transcript_id", "transcript_version",
                            "ensembl_gene_id", "external_gene_name", "description",
                            "transcript_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

#####################
## KO vs. WT 
#####################

#drop the hipp1 outlier 
#s2c <- filter(s2c, sample != "hipp1")

so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE)

## QC plots
plot_pca(so, color_by = 'geno', pc_y = 2, text_labels = FALSE)
plot_pca(so, color_by = 'tissue', pc_y = 2, text_labels = TRUE)
plot_sample_heatmap(so)

## fit models 
so <- sleuth_fit(so, ~geno + tissue, 'full')
so <- sleuth_fit(so, ~tissue, 'no_tissue')

so <- sleuth_lrt(so, 'no_tissue', 'full')
so <- sleuth_wt(so, 'genoko')
so <- sleuth_wt(so, 'tissuehip')
tests(so)

## EDA 
sleuth_live(so)

## get results 
so_tab_wt_geno <- sleuth_results(so, 'genoko')
so_tab_lrt_geno <- sleuth_results(so, 'no_tissue:full', 'lrt', show_all=TRUE)

## plotting by hand b/c the plotting function is broken ... arggh 

transcripts <- head(so_tab_wt_geno, 10)$target_id

tabd_df <- so$obs_norm[so$obs_norm$target_id %in% transcripts,]
tabd_df <- dplyr::select(tabd_df, target_id, sample, tpm)
tabd_df <- reshape2::dcast(tabd_df, target_id ~sample, value.var = 'tpm')

rownames(tabd_df) <- tabd_df$target_id
tabd_df$target_id <- NULL
trans_mat <- as.matrix(log(tabd_df + 1))

s2c <- so$sample_to_covariates
rownames(s2c) <- s2c$sample
annotation_cols = setdiff(colnames(so$sample_to_covariates), 'sample')
s2c <- s2c[, annotation_cols, drop = FALSE]

color_high = '#581845'
color_mid = '#FFC300'
color_low = '#DAF7A6'
colors <- colorRampPalette(c(color_low, color_mid, color_high))(100)

pheatmap::pheatmap(trans_mat, annotation_col = s2c, color = colors,
                        cluster_cols = TRUE)
