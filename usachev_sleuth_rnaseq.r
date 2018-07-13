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

#pseudoalignment analysis

library(sleuth)
library(cowplot)

#annotation
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library(dplyr)

#pathway
library(pathview)
library(gage)
library(gageData)

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
s2c <- read.table(file.path(base_dir, "for_sleuth_meta.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- arrange(s2c, sname)
colnames(s2c) <- c("sample","geno","tissue","path")
## add a column called "kal_dirs" containing paths to the data
s2c <- mutate(s2c, path = kal_dirs)

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

so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE, aggregation_column = "ens_gene")
plot_pca(so, color_by = 'geno', pc_y = 3, text_labels = FALSE)
plot_pca(so, color_by = 'tissue', pc_y = 3, text_labels = TRUE)

sleuth_live(so)

plot_loadings(so, pc_input=2)
#plot_bootstrap(so, 'ENSMUST00000082402.1', color_by = 'type')

plot_sample_heatmap(so)

so <- sleuth_fit(so, ~geno + tissue, 'full')
so <- sleuth_fit(so, ~tissue, 'no_tissue')

so <- sleuth_lrt(so, 'no_tissue', 'full')
so <- sleuth_wt(so, 'genowt')
so <- sleuth_wt(so, 'tissuehip')
tests(so)

