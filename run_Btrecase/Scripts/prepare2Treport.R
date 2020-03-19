
## Load R packages and functions:
library(data.table)


source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")

#' Merges stan output from 2 Tissues in one model for all genes run. 
#'
#' @param sum character vector with full names to stan summaries by gene
#' @param tag character vector with full names to stan tagging SNPs by gene
#' @param gtex character vector with full names to gtex files to use
#' @param out full name to save output
#' @keywords stan summary eQTL gtex
#' @export
#' @return saves file with merged and formated stan outputs
#' prep.sum()

prep.sum <- function(sum, tag,gtex, out){
    sum <- rbindlist(lapply(sum, fread))
    ## add gene_name (biomart)
    sum <- add.name(sum)
    ## add variable to code for joint significance in ba bd / bp bn
    sum <- add.sig(sum, type="Signif", cols=c("Signif.ba", "Signif.bd"), lab=c("ba", "bd"), newCol="Signif.ba.bd")
    sum <- add.sig(sum, type="Signif", cols=c("Signif.bp", "Signif.bn"), lab=c("bp", "bn"), newCol="Signif.bp.bn")

    ## check if the SNPs in sum have a reported association in GTEx skin sun or not sun exposed.
    ## also check if an association was reported in the psoriasis eQTL paper
    gtex.skin <- list.files(snakemake@params[['gtex_dir']], full.names=T)

    ## extract eQTL for the genes of interest in each dataset, then subset by tag, otherwise too slow

    gtexList <- lapply(gtex.skin, function(i) {
        dt=fread(cmd=paste0("zcat ", i, " | grep -E '", paste(unique(sum$Gene_id), collapse="|"), "' "),header=F)
        names(dt)  <- names(fread(cmd=paste0("zcat ", i, " | head -1 " )))
        return(dt)
    })
    names(gtexList) <- c("no_sun", "sun")

    ## format to match sum Gene_id and taf col
    gtex <- rbindlist(lapply(gtexList, function(i) {
        i[, Gene_id:=gsub("\\..*","", gene_id)]
        i[, Chrom:=gsub( "_.*", "", variant_id)]
        i[, SNP:=sub("^[0-9]*_", "", variant_id)]
        i[, SNP:=sub("_b37", "", SNP)]
        i[, SNP:=gsub("_", ":", SNP)]
    }), idcol="Sun_exposure")

    ## for each gene-variant pair in gtex select the row with min p value
    setkey(gtex,Gene_id, SNP, pval_nominal)

    gtex <- gtex[, .SD[1], .(Gene_id, SNP)]

    sum2 <- merge(sum, gtex, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), all.x=T)
    ## some tags are not in gtex data, look if any of the tagging snp can replace the missing value

    mis <- sum2[is.na(slope), .(Gene_id, tag)]

    tags <- rbindlist(lapply(tag, fread))

    ## add tag col in gtex (same as SNP) so I can use function tag.comp
    gtex[, tag:=SNP]

    ## get alternative SNP in gtex that is in "tags" for a gene and tag pair
    alt.tag <- rbindlist(tag.comp(DT1=mis, DT2=tags, DT3=gtex, s=".2T"))

    ## replace nas in sum2 with alternative SNP from gtex and add SNP col with SNP id
    sum2[, SNP:=tag]
    for (i in nrow(alt.tag)){
        sum2[Gene_id == alt.tag$Gene_id[i] & tag == alt.tag$tag.2T[i],
             names(gtex)[c(3:(ncol(gtex)-1),2)] := gtex[Gene_id ==alt.tag$Gene_id[i] & tag == alt.tag$tag[i],names(gtex)[c(3:(ncol(gtex)-1),2)], with=F] ]

    }

    
    

   write.table(sum2, out, row.names=FALSE)
    
}

    
prep.sum(sum= snakemake@input[['stan_sum2T']],
         tag = snakemake@input[['stan_tags2T']],
         gtex = snakemake@params[['gtex_dir']],
         out = snakemake@output[['summary']]
         )
