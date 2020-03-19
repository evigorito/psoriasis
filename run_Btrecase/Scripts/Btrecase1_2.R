#' ---
#' title: Btrecase on psoriasis data
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---


library(data.table)
library(ggplot2)

source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

## repeated brecase run
dir <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/output"
files <- list.files(path=dir, pattern = "ENSG[0-9]*\\.exfSNPp4.noGT.stan.summary.txt", recursive=TRUE, full.names=TRUE)

sum2 <- rbindlist(lapply(files, fread))
## add gene names
sum <- fread( "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/results/full.stan.summary.txt")
sum2  <- merge(sum2, sum[,.(Gene_id,Gene_name,tag,gene_dist, rs_id,FC,r2TagGwas)], by=c("Gene_id","tag"),
               all.x=TRUE)

gene_coord <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt"
geneStEnd <- fread(gene_coord)
plots2 <- lapply(unique(sum2$Gene_id), function(i) gene_plot(sum2,gene=i,geneStEnd, gwas="yes"))
plots2

## look at the correlation of counts for ERAP1 and CAST

counts <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts/Psoriasis_skin.txt")
counts <- counts[gene_id %in% unique(sum2[Gene_name %in% c("ERAP1", "CAST"), Gene_id]),]

counts.n <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts/normal_skin.txt")
counts.n <- counts.n[gene_id %in% unique(sum[Gene_name %in% c("ERAP1", "CAST"), Gene_id]),]

cor(t(counts.n[,-1]))
