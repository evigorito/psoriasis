library(data.table)


## Aux function for report.R


## aux function to look at variables in gwas by gene_dist
f <- function(dt, var1, var2=NULL, colname, d){
    var1.d <- t(sapply(d, function(i) summary(dt[gene_dist <=i,.N,var1][,N])))
    var1.Number <- sapply(d, function(i) length(unique(dt[gene_dist <=i,get(var1)])))
    
    if(!is.null(var2)){     
        var2.Number <- sapply(d, function(i) length(unique(dt[gene_dist <=i,get(var2)])))
        var1.d <- cbind(var1.d, var2.Number)
    }
    
    var1.d <- cbind(var1.d, var1.Number)  
    colnames(var1.d)[(ncol(var1.d)-length(colname) + 1):ncol(var1.d)] <- colname
    rownames(var1.d) <- paste0(d/1000,"_KB")
    return(var1.d)
}

## aux function to look at variables in drg by CaseMedian and FC
f2 <- function(dt, var1, var2, range1){
    var1.2 <-t(sapply(2:length(range1),
                       function(i) summary(dt[get(var1)> range1[(i-1)] & get(var1) <=range1[i],get(var2)])))
        
    Number.Genes <- sapply(2:length(range1),
                           function(i) nrow(dt[get(var1)> range1[(i-1)] & get(var1) <=range1[i],]))
    
    var1.2 <- cbind(var1.2,Number.Genes)
    rownames(var1.2) <- c("Low", "Medium", "High")
    return(var1.2)
}

## aux function to look at variables in drg by CaseMedian and distance to gene
f3 <- function(dt, gwas, var1, var2, d){
    exp <- t(sapply(d, function(i) summary(dt[get(var1) %in% unlist(gwas[gene_dist<=i, 'gene_name']) ,get(var2)])))
    
    Number.Genes=sapply(d, function(i) length(unique(dt[get(var1) %in% unlist(gwas[gene_dist<=i, 'gene_name']), get(var1)])))

    exp <- cbind(exp,  Number.Genes)
    rownames(exp) <- paste0(d/1000,"_KB")
    return(exp)

}

