log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Load R packages and functions:
library(data.table)
library(ggplot2)
library(xlsx)
library(parallel)
mc.cores = snakemake@threads

source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

#' Format Btrecase output to include GWAS data and DEG
#' 
#' @param dir directory with stan outputs to merge and format
#' @param pat pattern to match stan output, defaults to "ENSG[0-9]*\\.noGT.stan.summary.txt"
#' @param gwas full path to table with GWAS data, see snakefile 
#' @param drg full name to drg excel file, see snakefile
#' @param colclass character vector with col classes to use when openning excel drg file
#' @param legend character vector with full path to legend files for all autosomes
#' @param haps character vector with full path to hap files for all autosomes
#' @param out full name to output file to save formatted output
#' @export
#' @return text file with formatted output
#' format_btrecase()

format_btrecase <- function(dir, pat="ENSG[0-9]*\\.noGT.stan.summary.txt", gwas, drg, colclass, legend, haps, out){
    files <- list.files(path=dir, pattern = pat, recursive=TRUE, full.names=TRUE)

    sum <- rbindlist(lapply(files, fread))

    ## add gene_name (biomart)
    sum <- add.name(sum)
    
    ## open gwas file
    gwas <- fread(gwas)

    ## open drg file
    print(colclass)
    drg=data.table(read.xlsx2(drg,
                              sheetIndex=1 ,
                              colClasses=colclass))
    
    ## merge sum with gwas by gene_id (unique id, better than gene_name)
    sum <- merge(sum, gwas[,.(gene_id,CHROM, POS, Alleles, rs_id,gene_dist,table)], by.x="Gene_id", by.y="gene_id", sort=FALSE, all.x=T)

    ## merge sum with drg
    sum <- merge(sum, drg[, .(Gene.Symbol,RankP,CaseMedian, FC)], by.x="Gene_name",by.y="Gene.Symbol", all.x=T, sort=F)
    

    ## Add r2 between tag SNP and GWAS hit when applicable

    ## For each gene get CHROM, min and max position for the range to extract
    sum[, tag.pos:=as.numeric(gsub(":.*", "", tag))]
    setkey(sum, Gene_id, tag.pos)
    minPos <- sum[, min(c(tag.pos, POS)), Gene_id]
    maxPos <- sum[, max(c(tag.pos, POS)), Gene_id]
    minmax <- merge(minPos,maxPos, by="Gene_id", suffixes=c(".min", ".max"))

    ## get GWAS info from sum
    minmax <- merge(minmax, unique(sum[, .(Gene_id,CHROM, POS, Alleles, rs_id)]), by="Gene_id", all.x=T)

    ## remove NA (no GWAS hit)
    minmax <- minmax[!is.na(POS),]


    ## get r2 for tag SNPs and GWAS hit 

    tagr2 <-mclapply(1:nrow(minmax), function(i){
        file1 <- grep(paste0("chr",minmax[i,CHROM],".legend.gz"), legend,value=T)
        file2 <- grep(paste0("chr",minmax[i,CHROM],".hap.gz"), haps,value=T)
        cw <- as.numeric(unlist(minmax[i, .(V1.min, V1.max)]))
        mat <- haps.range(file1=file1, file2=file2, cw=cw, maf=0)
        rsnps <- sum[Gene_id == minmax[i,Gene_id], tag] 
        gwas <- grep(minmax[i,POS], rownames(mat), value=T)
        mat <- mat[c(rsnps,gwas),]
        r2 <- (cor2(t(mat)))^2
        r2 <- r2[gwas,rsnps,drop=FALSE]   
        return(r2) 
    },
    mc.cores=mc.cores)
    
    names(tagr2) <- minmax$Gene_id

    ## add tagr2 to sum
    for (i in seq_along(tagr2)){
        sum[Gene_id==names(tagr2)[i], r2TagGwas:=tagr2[[i]][1,] ]
    }

    write.table(sum, out, row.names=FALSE)
}

dir=snakemake@params[['direc']]
gwas=snakemake@input[['gwas']]
drg=snakemake@input[['drg']]
colclass=snakemake@params[['colclass']]
legend <- unlist(snakemake@input[['leRef']])
haps <- unlist(snakemake@input[['hapRef']])
out <- snakemake@output[[1]]

pat <- snakemake@params[['pattern']]
if(!exists("pat")) pat <- "ENSG[0-9]*\\.noGT.stan.summary.txt"

format_btrecase(dir, pat, gwas, drg, colclass, legend, haps, out)
