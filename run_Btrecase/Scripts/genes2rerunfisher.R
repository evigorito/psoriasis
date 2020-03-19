library(data.table)

#' Select genes with non-null associations but Fisher test below cut-off to re-run with Btrecase
#' @param sum full name to file with full summary from stan, extra cols allowed
#' @param min.p minimun p-value for Fisher test of fSNPs to re run test
#' @param out full name to file to write results
#' @keywords fisher genes
#' @export
#' @return saves text file with gene_id to re-run
#' genesFisher()

genesFisher <- function(sum, min.p,out){
    sum <- fread(snakemake@input[['summary']])
    u <- unique(sum[min.p.fsnp<= min.p & log2_aFC_null == "no", .(Gene_id)])
    write.table(u, out, row.names=F)
}

genesFisher(snakemake@input[['summary']], snakemake@params[['p']], snakemake@output[['fisher']])

