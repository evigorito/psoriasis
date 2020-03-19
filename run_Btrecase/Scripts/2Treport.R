#' ---
#' title: Comparing eQTL in normal vs psoriatic skin run in one or independent models
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule btrecase2T_report from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Snakefile

## Load R packages and functions:
library(data.table)
library(ggplot2)
library(cowplot)
library(parallel)


source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")

##' # Process outout from joint model

##' A joint model was formulated allowing different slope and
##' intercept for psoriasis and normal tissue and re-parametrised
##' defining ba = (bn+bp)/2 ; average and bd=(bp-bd)/2 mean
##' difference. The priors used were ba ~ N(0, 0.04) (same as single
##' tissue) and bd ~ N(0, 0.01), which corresponds to sd(ba)/2. We
##' expect bn and bp consistent in direction, mostly null, sometimes
##' same non-null effect but occasionally null in one tissue and
##' significant in the other.

sum <- fread(snakemake@input[['stan_sum2T']])
## sum <- fread('/mrc-bsu/scratch/ev250/psoriasis/Btrecase/results/normal_pso_JointModel_summary.txt')

## make significance columns based on coefficient estimates
sum <- add.sig(sum, type="Signif", cols=c("Signif.bp", "Signif.bn"), lab=c("bp", "bn"), newCol="Signif.bp.bn")
sum <- add.sig(sum, type="Signif", cols=c("Signif.ba", "Signif.bd"), lab=c("ba", "bd"), newCol="Signif.ba.bd")

## Compare output when tissues were run independently
sum.in <- fread(snakemake@input[['stan_sum']])
#'/mrc-bsu/scratch/ev250/psoriasis/Btrecase/results/normal_pso_ci99_summary.R'
## select same genes and tags

sum.in <- sum.in[Gene_id %in% sum$Gene_id & tag %in% sum$tag,]
## sum.in long to wide
skin=unique(sum.in$skin)

sum.inwide <- Reduce(function(x,y) merge(x,y, by=c("Gene_id","tag", "Gene_name","CHROM","POS","Alleles","rs_id", "gene_dist","table","RankP","CaseMedian", "FC","tag.pos","r2TagGwas","chrom"), suffixes=c(".bp", ".bn")), lapply(skin, function(i) sum.in[skin==i,]))

sum.inwide <- add.sig(sum.inwide, type="null", cols=c("log2_aFC_null.bp", "log2_aFC_null.bn"), lab=c("bp", "bn"), newCol="Signif.bp.bn")

##' # Compare Joint model with independent models

joint.babd <- ggplot(sum, aes(x=log2_aFC_mean.bp, y=log2_aFC_mean.bn, color=Signif.ba.bd)) + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("Joint Model") + geom_hline(yintercept=0, linetype="dashed", color= "gray") + geom_vline(xintercept=0, linetype="dashed", color= "gray")

joint.bpbn <- ggplot(sum, aes(x=log2_aFC_mean.bp, y=log2_aFC_mean.bn, color=Signif.bp.bn)) + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("Joint Model")+ geom_hline(yintercept=0, linetype="dashed", color= "gray") + geom_vline(xintercept=0, linetype="dashed", color= "gray")

indep.bpbn <- ggplot(sum.inwide, aes(x=log2_aFC_mean.bp, y=log2_aFC_mean.bn, color=Signif.bp.bn)) + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("Two models")+ geom_hline(yintercept=0, linetype="dashed", color= "gray") + geom_vline(xintercept=0, linetype="dashed", color= "gray")


plot_grid(joint.bpbn, indep.bpbn, joint.babd, ncol=2)


## For each gene get the SNP with min p-value in gtex

setkey(sum, Gene_id, pval_nominal)

sum[!is.na(pval_nominal) , .SD[1], .(Gene_id)][, c("Gene_name", "tag", grep("Signif\\.[a-z]*\\.[a-z]*", names(sum), value=T), "tss_distance", "maf", "pval_nominal"), with=F]


## look for the most significant SNP associated with each gene in previous psoriasis eQTL data
eQTL <- snakemake@input[['eqtl_pso']]
names.eqtl <- c("pso"=grep("PP53", eQTL) ,"normal"=grep("NN57", eQTL))
eqtl_pso <- lapply(eQTL[names.eqtl], fread)
names(eqtl_pso) <- names(names.eqtl)
eqtl <- rbindlist(eqtl_pso, idcol="skin")
eqtl <- eqtl[Gene.Symbol %in% unique(sum$Gene_name),]
setkey(eqtl, Gene.Symbol, PVALUE)
eqtl[, .SD[1], Gene.Symbol]


## look at plots per gene
gene_coord <- snakemake@input[['gene_coord']]
geneStEnd <- fread(gene_coord)

plots1M <- lapply(unique(sum$Gene_id), function(i) gene_plot_2T(sum,gene=i,geneStEnd,coef=c("bn", "bp"),
                                                                ci=c(0.5,99.5),
                                                                sig="Signif.ba.bd") )

names(plots1M) <- unique(sum$Gene_id)

## Plots for independent models
plots <- lapply(unique(sum.in$Gene_id), function(i) gene_plot_skin(sum.in,gene=i,geneStEnd,
                                                                  gwas="yes", var="skin",ci=c(0.5,99.5)))

names(plots) <- unique(sum.in$Gene_id)

## select genes run in joint model

plots <- plots[unique(sum$Gene_id)]

lapply(names(plots1M), function(i) plot_grid(plots1M[[i]], plots[[i]], nrow=2) )
