log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Load R packages and functions:
library(data.table)
library(ggplot2)
library(parallel)
mc.cores = snakemake@threads

source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R", verbose=FALSE)
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R', verbose=FALSE)


#' Prepare summary of associations combining Normal and psoriatic skin and compute r2 between rsnps per gene
#'
#' @param pso95 summary already made when running psoriasis with CI95, to get gene names, GWAS snp, r2 GWAS snp and tag snps
#' @param gene_coord text file with cols gene_id and chrom, other cols allowed
#' @param lefiles character vector with full names to reference panel legend files per chromosome
#' @param hapfiles character vector with full names to reference panel hap files per chromosome
#' @param maf value to filter snps, defaults to 0.05 which was used to select rsnps to run
#' @param stan.dir path to directory with input stan summaries and tags
#' @param r2 full name to file to save list of matrices with r2 between rsnps per gene
#' @param stan_sum full name to file to save stan summary for all genes run both in normal and psoriatic skin
#' @return save 2 files, list of matrices with r2 between rsnps per gene and stan summary combining normal and psoriatic skin
#' @export
#' skin_sumr2()

skin_sumr2 <- function(pso95, gene_coord,  lefiles, hapfiles, maf=0.05, stan.dir, r2, stan_sum){

    skin=c("Psoriasis_skin", "normal_skin")
    out.files <- lapply(skin, function(i) list.files(path=stan.dir, pattern=paste0("ENSG[0-9]*\\.", i,".noGT.stan.summary.txt"), recursive=TRUE, full.names=TRUE))
    names(out.files) <- skin

    all.tags <- lapply(skin, function(i) list.files(path=stan.dir, pattern=paste0("ENSG[0-9]*\\.", i,".noGT.tags.lookup.txt"), recursive=TRUE, full.names=TRUE))
    names(all.tags) <- skin
    
    fsnps.fish <- lapply(skin, function(i) list.files(path=stan.dir, pattern=paste0("ENSG[0-9]*\\.", i,".fsnps.het.fisher.test.txt"), recursive=TRUE, full.names=TRUE))
    names(fsnps.fish) <- skin

    

    ## Summary of association per skin type
    sum99 <- rbindlist(lapply(out.files,function(i) rbindlist(lapply(i,fread))), idcol="skin")

    ## Get summary output for Btrecase in psoriatic skin with CI95
    pso95 <- fread(pso955)
    
    ## Extract gene name, GWAS hit and r2 of tag to gwas hit plus expression FC from pso95
    sum99 <- merge(sum99, pso95[,c("Gene_id","Gene_name","tag","CHROM","POS", "Alleles","rs_id",
                               "gene_dist","table","RankP","CaseMedian","FC", "tag.pos","r2TagGwas")],
               by=c("Gene_id","tag"),
               all.x=T,
               sort=F)

    ## Add chrom to stan99
    gene.coord <- fread(gene_coord)
    sum99 <- merge(sum99, gene.coord[,.(gene_id,chrom)], by.x="Gene_id", by.y="gene_id")

###### Look for fsnps as rsnps in sum99
    ## get tags
    tags <- rbindlist(lapply(all.tags,function(i) rbindlist(lapply(i,fread))), idcol="skin")
    tags.run <- merge(tags, unique(sum99[,.(chrom, Gene_id,tag, skin)]), by= c("Gene_id","tag", "skin"))

    ## get fSNPs fisher output and add EAF
    fsnps <- rbindlist(lapply(fsnps.fish,function(i) rbindlist(lapply(i,fread))), idcol="skin")
    fsnps <- merge(fsnps, gene.coord[,.(gene_id,chrom)], by="gene_id")
    setkey(fsnps, chrom)
    eaf <- mclapply(unique(fsnps$chrom), function(i) snp.eaf(file1=grep(paste0("chr",i,".legend.gz"), lefiles,value=T),
                                                             snps=unique(fsnps[chrom==i,fsnp])),
                    mc.cores=mc.cores)
    names(eaf) <- unique(fsnps$chrom)
    fsnps <- merge(fsnps, rbindlist(eaf,idcol="chrom"), by.x=c("chrom", "fsnp"), by.y=c("chrom", "snp"))
    setnames(fsnps, c("fsnp","eaf"), c("fSNP.id","fSNP.EAF"))
        

    ## merge fsnps with tags.run to get tagging snp for fsnp
    fsnps.run <- merge(fsnps,tags.run,  by.y=c("Gene_id", "SNP","skin", "chrom"), by.x=c("gene_id","fSNP.id","skin", "chrom"), all.x=T)
    fsnps.run <- fsnps.run[!is.na(tag),][, c("OR", "pvalue"):=NULL]

    ## add fSNPs to sum99

    sum99 <- merge(sum99, fsnps.run, by.x=c("Gene_id", "skin","chrom", "tag"), by.y=c("gene_id", "skin","chrom", "tag"), all.x=T)
    sum99[!is.na(fSNP.id) , tag.fsnp.op:="no"][(fSNP.EAF < 0.5 & tag.EAF > 0.5) | (fSNP.EAF > 0.5 & tag.EAF < 0.5), tag.fsnp.op:="yes"]
    
    ## save sum99
    write.table(sum99, stan_sum, row.names=F)

    ## get r2, add tag.pos to sum99
    sum99[, tag.pos:=as.numeric(gsub(":.*", "", tag))]
    setkey(sum99, chrom, Gene_id, tag.pos)

    tagr2 <- mclapply(unique(sum99$Gene_id), function(i){
        dt <- sum99[Gene_id==i,]
        file1 <- grep(paste0("chr",unique(dt[,chrom]),".legend.gz"), lefiles,value=T)
        file2 <- grep(paste0("chr",unique(dt[,chrom]),".hap.gz"), hapfiles,value=T)
        cw <- c(min(dt$tag.pos), max(dt$tag.pos))
        mat <- haps.range(file1=file1,file=file2,cw=cw, maf=maf)
        mat <- mat[unique(dt$tag),]
        r2 <- (cor2(t(mat)))^2
        return(r2)
    },
    mc.cores=mc.cores)
    names(tagr2) <- unique(sum99$Gene_id)

    ## save tagr2
    saveRDS(tagr2, r2)    

}

pso95 <- snakemake@input[['pso95']]
gene_coord <- snakemake@input[['gene_coord']]
lefiles <- unlist(snakemake@input[['leRef']])
hapfiles <- unlist(snakemake@input[['hapRef']])
stan.dir <- snakemake@params[['out_dir']]
r2 <- snakemake@output[['r2']]
stan_sum <- snakemake@output[['stan_sum']]

skin_sumr2(pso95, gene_coord, lefiles, hapfiles, maf=0.05, stan.dir, r2, stan_sum)

## gene_coord="/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt"
