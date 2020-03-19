#' ---
#' title: Btrecase on psoriasis data
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule Btrecase_analysis from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Snakefile

## Load R packages and functions:
library(data.table)
library(ggplot2)

source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')


#' ## Summary of associations

#' Get summary output for the genes run by Btrecase and format_btrecase_output
sum <- fread(snakemake@input[['stan']])

#' Number of genes and distribution of the number of snps per gene run by Btrecase
length(unique(sum$Gene_id))
summary(sum[, .N, Gene_id][,N])

#' Number of genes with non-null associations and distribution of the number of snps with non-null associations

length(unique(sum[log2_aFC_null== "no",Gene_id]))
summary(sum[log2_aFC_null== "no", .N, Gene_id][,N])

#' ## QC

#' * Histograms of log2 fold change for all SNPs tested
ggplot(sum, aes(x=log2_aFC_mean)) +  geom_histogram(color="black", fill="white") + theme_bw()

#' * Correlation between dispersion of the log2 FC per gene and number of fSNPs per gene
sd.FC <- sum[, .(sd.log2FC=sd(abs(log2_aFC_mean))), Gene_id]
sd.FC <- merge(sd.FC, unique(sum[, .(Gene_id, n.fsnps)]), by="Gene_id")
ggplot(sd.FC, aes(y=sd.log2FC, x=n.fsnps)) + geom_point(shape=1) +
    geom_smooth(method=lm,  se=FALSE) +
    ggtitle("SD for log2FC per gene vs number of fSNPs") +
    theme_bw()

#' * Correlation between dispersion of the log2 FC per gene and EAF of tested SNPs

sd.EAF <- sum[, .(sd.EAF=sd(tag.EAF)), Gene_id]
sd.FC <- merge(sd.FC, sd.EAF, by="Gene_id")
ggplot(sd.FC, aes(y=sd.log2FC, x=sd.EAF)) + geom_point(shape=1) +
    geom_smooth(method=lm,  se=FALSE) +
    ggtitle("SD for log2FC per gene vs SD EAF") +
    theme_bw()


#' Plot associations by gene
gene_coord <- snakemake@input[['gene_coord']]
geneStEnd <- fread(gene_coord)
plots <- lapply(unique(sum$Gene_id), function(i) gene_plot(sum,gene=i,geneStEnd, gwas="yes"))

## some interesting genes to follow up from GWAS side:
## 1- IKBKE, NFKB1Z and REL ( this one a bit unlucky with r2, sig) but the 3 are components of the NFkB pathway implicated in psoriasis.
## 2- ERAP1, process peptides for HLA classI pathway
## 3-CAST has been implicated in skin diseases
## 4-RASIP1 ras pathway always important for signalling, but not super convincing
## 5-FUT2 more convincing but unclear link to psoriasis


#' ## Plots per gene
#' * Blue dotted lines indicate gene boundaries
#' * min.p the minimun p value for the Fisher test of proportion of het fSNPs in sample and reference panel

plots


