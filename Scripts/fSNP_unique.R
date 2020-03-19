library(data.table)

#' Select unique fSNPs to genes when no strand information is available
#'
#' @param fsnps file with fSNP coordinates mapped to gene_id
#' @param out file with unique fSNP per gene
#' @return saves file with unique fSNPs to genes
#' @export
#' u.fsnps()

u.fsnps <- function(fsnps, out){
    comb <- fread(fsnps)
    ## remove duplicated fSNPs
    comb <- comb[ !ID %in% comb[duplicated(ID),ID],]
    write.table(comb, file=out, row.names=F)
    
}



u.fsnps( fsnp=snakemake@input[['fsnps']],
        out=snakemake@output[['ufsnps']])
