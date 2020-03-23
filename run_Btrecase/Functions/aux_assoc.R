library(data.table)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(mygene)
library(ggrepel)
library(cowplot)
library(biovizBase)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
library(biovizBase)
library(GenomicAlignments)
library(GenomicFeatures)
library(ggbio)
library(annotate)
library(parallel)
library(biomaRt)

source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")

#' Plot to visualize gene-snps associations
#'
#' @param sum data table with stan output 
#' @param gene character vector with gene ID
#' @param geneStEnd data table  with gene_id, gene start and end columns
#' @param gwas whether to highlight the SNP in sum with the highest r2 with the gwas hit, when there is one
#' @keywords stan plot gene-snps compare
#' @export
#' @return ggplot object
#' gene_plot()
gene_plot <- function(sum,gene, geneStEnd, gwas=NULL){
    dt <- sum[Gene_id==gene,]
    stend <- unlist(geneStEnd[gene_id==gene, .(start,end)])
    ## create col for position, use tag.ngt because it has more entries
    dt[, Position:=as.numeric(gsub(":.*","", tag))]
    dt[, log2_aFC_null:=factor(log2_aFC_null, levels=c("yes", "no"))]

    ## modify cred int when necessary

    dt[,abs.ci.low:=ifelse(`log2_aFC_mean`> 0, `log2_aFC_2.5%`, -`log2_aFC_2.5%`)]
    dt[,abs.ci.h:=ifelse(`log2_aFC_mean`> 0, `log2_aFC_97.5%`, -`log2_aFC_97.5%`)]
    ## set colors for consistency across plots
    colors <- c("#D55E00","#009E73")
    names(colors) <- c("yes","no")
    man.col <- colors[names(colors) %in% dt$log2_aFC_null]

    ## get gene name for ggtitle
    geneName <- unique(dt[["Gene_name"]])

    ## get source of gene
    sourceGene <- ifelse(is.na(unique(dt$gene_dist)), paste("DRG/FC",round(unique(dt$FC),2)), paste("GWAS",unique(dt$rs_id) ))
    
    p <- ggplot(dt, aes(x=Position, y=abs(log2_aFC_mean))) +
        geom_point(aes(colour=`log2_aFC_null`)) +
        scale_colour_manual(values=man.col) +
        geom_hline(yintercept=0) +
        geom_errorbar(data=dt, aes(ymin=abs.ci.low, ymax=abs.ci.h), color="grey",linetype="dashed",width=0.1) +
        ggtitle(paste(geneName,"from", sourceGene, "\nmedian info=", round(median(dt$info),2),", min.p=",formatC(min(dt$min.p.fsnp),0))) +
        ylab(label="|log2(aFC)|") +
        geom_vline(xintercept=stend, linetype="dotted", color="blue")+
        theme_bw()

    if(!is.null(gwas)){
        ## gwas hit
        if(any(!is.na(dt$r2TagGwas))){
            dt[, r2:=""][r2TagGwas==max(r2TagGwas), r2:= paste0("r2= ",as.character(round(r2TagGwas,2)))]
            p <- p + geom_text_repel(aes(label=r2), size=4, data=dt) 
        }
    }
    return(p)
    
}
    

#' Add gene name by gene id
#'
#' @param sum data table with stan output 
#' @keywords gene name
#' @export
#' @return data table with gene_name column
#' add.name()
add.name <- function(sum){
    ##mart <-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
    ##mart <-  useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ## g.name <- data.table(getBM(attributes=c("ensembl_gene_id","external_gene_name"),
    ##                             filters="ensembl_gene_id",
    ##                             values=unique(sum[["Gene_id"]]),
    ##                            mart=mart))
    g.name <- as.data.table(getGenes(unique(sum[['Gene_id']]), fields='symbol'))[, .(query,symbol)]
    setnames(g.name, names(g.name), c("Gene_id", "Gene_name"))
    DT <- merge(sum, g.name, by="Gene_id", sort=F, all.x=T)
    return(DT)
}

#' Add rsid 
#'
#' @param sum data table with stan output, requires tag and chrom column
#' @param col with variable POS:REF:ALT, defaults to tag
#' @keywords rsid
#' @export
#' @return data table with rsid column
#' add.rsid()
add.rsid <- function(sum, pos="tag"){
    snpmart = useEnsembl(biomart = "snp", dataset="hsapiens_snp", GRCh=37)
    ## Need to make small intervals to call variants from biomart. I will groups snps per gene.
    dt1 <- unique(sum[,.(Gene_id, tag, chrom)])
    dt1[, Position:=as.numeric(gsub(":.*","", get(pos)))]
    dt1[, allele:=gsub("[0-9]+:","",get(pos))][, allele:=gsub(":", "/", allele)]

    ## Prepare list for getting rsid
    temp <- unique(dt1[,.(chr=chrom, start=min(Position), end=max(Position)),Gene_id])

    rs.id <- mclapply(1:nrow(temp), function(i){
        dt <- data.table(getBM(attributes = c('refsnp_id','allele', "chrom_start"), 
                         filters = c('chr_name','start','end'), 
                         values =list(temp[i, chr], temp[i,start], temp[i, end]),
                         mart=snpmart))                      
        tmp <- merge(dt1[Gene_id==temp[i,Gene_id],], dt, by.x=c("Position", "allele"), by.y=c("chrom_start", "allele"))
        return(tmp)
    }, mc.cores=parallel::detectCores())
    ## exclude failures
    tmp <- rbindlist(rs.id[which(sapply(rs.id, class)!="try-error")])

    sum <- merge(sum, unique(tmp[,.(get(pos), refsnp_id)]), by.x=pos, by.y="V1",all.x=T)
    return(sum)
       
}
    
   

#' Plot to compare gene-snps associations by a categorical variable
#'
#' @param sum99 data table with stan output 
#' @param gene character vector with gene ID
#' @param pos column name to use for position for x axis, defaults to tag
#' @param geneStEnd data table  with gene_id, gene start and end columns, defaults to NULL
#' @param gwas whether to highlight the SNP in sum with the highest r2 with the gwas hit, when there is one
#' @param title if all it includes gene name plus source of gene (GWAS or fold change) plus median info plus min.p value for fisher test, otherwise just gene name
#' @param yaxis character with y-axis label
#' @param var categorical variable to compare output by, defaults to skin
#' @param ci numeric vector with extremes for ci, defaults to 95%
#' @param null name for null column, defaults to log2_aFC_null
#' @param yvar name of variable to plot in y axis, defaults to log2_aFC_mean
#' @param colvar name of variable to use for color, defaults to Signif
#' @param colors character vector with colors to use, names factor in colvar, values colors, defaults to Signif
#' @param sizevar name of variable to use for size, defaults to NULL
#' @param size character vector with size to use, names factor in sizevar, values size, defaults to NULL
#' @param size.leg whether to show size legend, defaults to TRUE
#' @param shapevar name of variable to use for shape, defaults to NULL
#' @param shape character vector with shape to use, names factor in shapevar, values shape, defaults to NULL
#' @param shape.leg whether to show shape legend, defaults to TRUE
#' @param geneCoord to use with d for setting x lim, same data table as geneStEnd
#' @param d distance from gene to plot (gene start-d, gene end+d), defaults to NULL
#' @param rsid file with variant information formatted as ensembl ftp.ensembl.org/pub/grch37/current/variation/gvf/homo_sapiens/homo_sapiens-chr22.gvf.gz for appropiate built and chromosome, defaults to none
#' @keywords stan plot gene-snps skin
#' @export
#' @return ggplot object
#' gene_plot_skin()
gene_plot_skin <- function(sum99, gene, pos="tag", geneStEnd=NULL, gwas=NULL, title="all", yaxis="eQTL effect" ,var="skin",  ci=c(2.5,97.5), null=NULL, yvar="log2_aFC_mean" , colvar="Signif", colors= setNames(c("steelblue4","gold3"), c("Yes","No")) , sizevar=NULL, size=NULL, size.leg="legend", shapevar=NULL, shape=NULL, shape.leg="legend", d=NULL,geneCoord=NULL){

    dt <- sum99[Gene_id==gene,]
   
    ## create col for position
    dt[, Position:=as.numeric(gsub(":.*","", get(pos)))/10^6]
           
    ## create col Signif instead of log2_aFC_null
    if(colvar=="Signif"){
        if(! "Signif" %in% names(dt)){
            if('log2_aFC_null' %in% names(dt)) {
                dt[, log2_aFC_null:=factor(log2_aFC_null, levels=c("yes", "no"))]
            } else {
                dt[, log2_aFC_null:=factor(get(null), levels=c("yes", "no"))]
            }
            dt[, Signif:=ifelse(log2_aFC_null=="no","Yes","No")]
        } else {
            
            ## change Signif to Yes, No if ti was already in dt
            dt[, Signif:=ifelse(Signif=="no","No","Yes")]

        }
    }

    if(!is.null(ci)){
        ## get cred int 
        ci.low <-paste0("log2_aFC_",ci[1], "%")
        ci.high <- paste0("log2_aFC_",ci[2], "%")
        
        
        dt[,ci.low:= get(ci.low)]
        dt[,ci.high:=get(ci.high)]


    }
    ## set colors for consistency across plots
    ##colors <- c("steelblue4","gold3") #c("#D55E00","#009E73")
    ##names(colors) <- c("Yes","No")
    man.col <- colors[names(colors) %in% dt[[colvar]]]

    ## get gene name for ggtitle
    geneName <- unique(dt[!is.na(Gene_name), Gene_name])

    ## start plot

    if(is.null(shapevar)) {
        if(is.null(sizevar)) {
            p <- ggplot(dt, aes(x=Position, y=get(yvar))) +
                guides(colour=guide_legend(override.aes= list(shape=1)))
        } else {
            man.size <- size[names(size) %in% dt[[sizevar]]]
            p <- ggplot(dt, aes(x=Position, y=get(yvar), size=get(sizevar))) +         
                scale_size_manual(values=man.size) +
                guides(colour=guide_legend(override.aes= list(shape=1)),
                       size=size.leg) +
                geom_text_repel(aes(label=get(sizevar)))
        } 
            
    } else {
        man.shape <- shape[names(shape) %in% dt[[shapevar]]]
        if(is.null(sizevar)) {
            p <- ggplot(dt, aes(x=Position, y=get(yvar), shape=get(shapevar))) +
                scale_shape_manual(values=man.shape)  +
                guides(colour=guide_legend(override.aes= list(shape=1)),
                       shape=shape.leg)
            
        } else {
            man.size <- size[names(size) %in% dt[[sizevar]]]
            p <- ggplot(dt, aes(x=Position, y=get(yvar), size=get(sizevar),shape=get(shapevar))) +         
                scale_size_manual(values=man.size) +
                scale_shape_manual(values=man.shape)+
                guides(size=size.leg,
                       colour=guide_legend(override.aes= list(shape=1)),
                       shape=shape.leg) +
                geom_text_repel(aes(label=get(sizevar)))
        } 
        
    }
    
    if(!is.null(ci)) {
        p <- p + geom_errorbar(data=dt, aes(ymin=ci.low, ymax=ci.high),
                               color="grey83",linetype="dashed",width=0, size=.2)
    } else {
        p <- p + ylim(NA, max(dt[[yvar]]) + 0.2*max(dt[[yvar]]) )
    }
    
    
        
    p <-  p + scale_colour_manual(values=man.col) +
        geom_hline(yintercept=0, linetype="dashed", color="gray") +
        geom_point(aes(colour=get(colvar))) +
        labs(color = colvar,
             shape=shapevar) +
        ylab(label=yaxis) +
        xlab(label="Position (MB)") +
        
        theme_bw()

    if(!is.null(gwas)){
        ## gwas hit
        if(any(!is.na(dt$r2TagGwas))){
            dt[, r2:=""][r2TagGwas==max(r2TagGwas,na.rm=T), r2:= paste0("GWAS r2= ",as.character(round(r2TagGwas,2)))]
            ## only keep one max per var, in case of repeated values
            row.max <- lapply(levels(as.factor(dt[[var]])), function(i) dt[r2!="" & get(var)==i, which=T])
            if (any(lapply(row.max, length) >1)){
                rem <- unlist(lapply(row.max, function(i) {
                    if( length(i) > 1){
                        i[2:length(i)]
                    }
                }))
                dt <- dt[!rem,]
                }                
                
            p <- p + geom_text_repel(aes(label=r2), size=3 ,nudge_y=min(dt$log2_aFC_mean),
                                     direction="x",
                                     segment.color="red",
                                     data=dt) 
        }
    }
    if(!is.null(geneStEnd)) {
        stend <- unlist(geneStEnd[gene_id==gene, .(start,end)])/10^6
        p <- p + geom_vline(xintercept=stend, linetype="dotted", color="#009E73")
    }

    if(!(is.null(geneCoord) & is.null(d))){
        xl <-  (unlist(geneCoord[gene_id==gene, .(start,end)]) + c(-d,d)) /10^6
        p <- p + xlim(xl)
    }
    

    if(title=="all"){
        ## get source of gene
        sourceGene <- ifelse(all(is.na(unique(dt$gene_dist))), paste("DRG/FC",round(unique(dt[!is.na(FC),FC]),2)), paste("GWAS",unique(dt$rs_id) ))
        p <- p + ggtitle(paste(unique(dt[,get(var)]), geneName,"from", sourceGene, "\nmedian info=", round(median(dt$info),2),", min.p=",formatC(min(dt$min.p.fsnp),0)))
    } else {
        p <- p + ggtitle(geneName)
    }

   
    
    if(length(unique(dt[,get(var)])) > 1){
            
         p <- p + facet_grid(as.formula(paste(var, '~.'))) +
            theme(plot.title = element_text(size = 12),axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5), strip.background=element_rect(fill="white"))
    } else {
        p <- p + theme(plot.title = element_text(size = 12),axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5), strip.background=element_rect(fill="white"))
    }
    
  
    return(p)
    
}


#' helper function to add rsid for top snp to use shape and size options in gene_plot_skin
#'
#' @param x data table with the data to produce p
#' @param gene character vector with gene ID
#' @param var categorical variable to compare output by, defaults to skin
#' @param null name for null column, defaults to log2_aFC_null
#' @param rsid file with variant information formatted as ensembl ftp.ensembl.org/pub/grch37/current/variation/gvf/homo_sapiens/homo_sapiens-chr22.gvf.gz for appropiate built and chromosome, defaults to none
#' @export
#' @return data table with rsid for top SNP by skin as extra column
#' rsid2plot()

rsid2plot <- function(x, gene, var="skin", null="log2_aFC_null" , rsid){
    x <- x[Gene_id==gene,]
    ## get top tag by var
    top.tag <- lapply(levels(as.factor(x[[var]])), function(i) x[get(null)=="no" & get(var) ==i, .(tag, log2_aFC_d)][log2_aFC_d == max(log2_aFC_d), tag])
    ## if any (significant), get rsid
    if(length(unlist(top.tag))){
            rs <- get.rsid(rsid,snps= unique(unlist(top.tag)))
    }

    tmp <- mapply(function(i,j) {
            x[get(var) ==j, rsid:=""]
            if(length(i)){
                x[tag==i & get(var)==j, rsid:=rs[SNP==i,rsid]]
            }},          
            i=top.tag,
            j=levels(as.factor(x[[var]])), SIMPLIFY=F)
    
    ## lev.snp <- lapply(skin, function(i) levels(as.factor(x[skin==i, rsid])))
    ## shape <- lapply(lev.snp, function(i) {
    ##     if(length(i) ==2){
    ##         shape <- rep(1,2)
    ##         shape[which(i != "")] <- 23
    ##     } else {
    ##         shape <- 1
    ##     }
    ##     names(shape) <- i
    ##     return(shape)
    ## })

    ## size.p <- lapply(lev.snp, function(i) {
    ##     if(length(i) ==2){
    ##         shape <- rep(1.5,2)
    ##         shape[which(i != "")] <- 3
    ##     } else {
    ##         shape <- 1.5
    ##     }
    ##     names(shape) <- i
    ##     return(shape)
    ## })

    return(x)
   

}

#' Plot to compare gene-snps associations by a categorical variable run in 1 MODEL
#'
#' @param sum99 data table with stan output from 2 tissues model
#' @param gene character vector with gene ID
#' @param geneStEnd data table  with gene_id, gene start and end columns
#' @param gwas whether to highlight the SNP in sum with the highest r2 with the gwas hit, when there is one
#' @param coef character vector with coefficients to plot 
#' @param ci numeric vector with extremes for ci, defaults to 95%
#' @param sig variable to color plot by
#' @keywords stan plot gene-snps skin
#' @export
#' @return ggplot object
#' gene_plot_2T()
gene_plot_2T <- function(sum99, gene, geneStEnd, gwas=NULL, coef, ci=c(2.5,97.5), sig){
    dt <- sum99[Gene_id==gene,]
    stend <- unlist(geneStEnd[gene_id==gene, .(start,end)])/10^6
    ## create col for position, use tag.ngt because it has more entries
    dt[, Position:=as.numeric(gsub(":.*","", tag))/10^6]

    ## get cred int
    cred <- as.character(outer(ci, coef, function(x,y) paste0("log2_aFC_",x, "%.",y)))

    ## get y axis

    y <- paste0("log2_aFC_mean.",coef)
    names(y) <- coef

    ## set colors for consistency across plots
    colors <- c("#999999", "#F0E442", "#0072B2", "#D55E00")
    names(colors) <- levels(dt[[sig]])
    man.col <- colors[names(colors) %in% dt[[sig]]]

    ## get gene name for ggtitle
    geneName <- unique(dt[!is.na(Gene_name), Gene_name])

    ## make same y axis
    y.max=max(unlist(lapply(cred, function(i) dt[[i]])))
    y.min=min(unlist(lapply(cred, function(i) dt[[i]])))

    p <- lapply(seq_along(coef), function(i) ggplot(dt,  aes(x=Position, y=get(y[i]))) +
                                             scale_colour_manual(values=man.col,name=sig) +
                                             geom_hline(yintercept=0) +
                                             geom_errorbar(data=dt, aes(ymin=get(grep(coef[i], cred, value=T)[1]), ymax=get(grep(coef[i], cred, value=T)[2])), color="grey83",linetype="dashed",width=0.0) +
                                             geom_point(aes(colour=get(sig))) +
                                             ylab(label="log2(aFC)") +
                                             xlab(label="Position (MB)") +
                                             geom_vline(xintercept=stend, linetype="dotted", color="#009E73")+
                                             theme_bw()  +
                                             scale_y_continuous(limits = c(y.min,y.max)) +
                                             
                                             ggtitle(paste0(coef[i],"; ", geneName)) 
                                            
                )
    return(plot_grid(plotlist=p))
    
}

#' Facet plot to compare gene-snps associations by a categorical variable run in 1 MODEL
#'
#' @param sum data table with stan output from 2 tissues model
#' @param gene character vector with gene ID
#' @param geneStEnd data table  with gene_id, gene start and end columns
#' @param d distance from gene to plot (gene start-d, gene end+d), defaults to NULL
#' @param gwas whether to highlight the SNP in sum with the highest r2 with the gwas hit, when there is one
#' @param yaxis character with y-axis label
#' @param yvar name of variable to plot in y axis, defaults to log2_aFC_mean
#' @param colvar name of variable to use for color, defaults to Signif
#' @param colors character vector with colors to use, names factor in colvar, values colors, defaults to Signif
#' @param sizevar name of variable to use for size, defaults to NULL
#' @param size character vector with size to use, names factor in sizevar, values size, defaults to NULL
#' @param size.leg whether to show size legend, defaults to TRUE
#' @param shapevar name of variable to use for shape, defaults to NULL
#' @param shape character vector with size to use, names factor in shapevar, values shape, defaults to NULL
#' @param shape.leg whether to show shape legend, defaults to TRUE
#' @param coef character vector with coefficients to facet plot 
#' @param ci numeric vector with extremes for ci, defaults to 95%
#' @param null name for null column, defaults to log2_aFC_null
#' @param sig variable to color plot by
#' @keywords stan plot gene-snps skin
#' @export
#' @return ggplot object
#' facet_plot_2T()

facet_plot_2T <- function(sum, gene, geneStEnd=NULL, d=NULL, gwas=NULL, yaxis="eQTL effect" ,yvar="log2_aFC_mean",  coef, ci=c(2.5,97.5), sig="Signif", title=NULL, null=NULL,  colvar="Signif", colors= setNames(c("steelblue4","gold3"), c("Yes","No")) , sizevar=NULL, size=NULL, size.leg="legend", shapevar=NULL, shape=NULL, shape.leg="legend", geneCoord=NULL){
    
    dt <- sum[Gene_id==gene,]

    dt.long <- rbindlist(lapply(coef, function(i) {
        names.coef <- grep(i, names(dt), value=T)
        tmp <- dt[ , c("Gene_id","Gene_name",  "tag", names.coef), with=F]
        setnames(tmp, names.coef, gsub(paste0(".", i), "", names.coef))
        tmp[, coef:=i]
        return(tmp)
        }))
    ## in dt.long Signif.bp corresponds to Signif.bp.bn 
    
    p <- gene_plot_skin(dt.long, gene, pos="tag", geneStEnd=geneStEnd, gwas=gwas, title="no", yaxis=yaxis, yvar=yvar, var="coef", ci=ci, null=null, d=d, colvar=colvar, colors=colors, sizevar=sizevar,size=size,size.leg=size.leg, shapevar=shapevar, shape=shape, shape.leg=shape.leg, geneCoord=geneCoord)
    
    return(p+ggtitle(title))



}
    ## covert dt to long
  


#' Add Signif column based on two columns
#'
#' @param dt data table with individual significance or null columns
#' @param type whether the columns to use are coded as "Significant" or doesnt contain the "null"
#' @param cols name of columns to based signif col
#' @param lab labels to use in new col
#' @param newCol name for new column, defaults to Signif
#' @keywords stan plot gene-snps skin
#' @export
#' @return data table with new colum newCol coded as factor with order "None, lab, "Both"
#' add.sig()

add.sig <- function(dt, type=c("Signif", "null"), cols, lab, newCol="Signif"){

    
    
    if(type=="Signif"){
        code <- c(sig = "yes", nosig = "no")
    } else {
        code <- c(sig = "no", nosig = "yes") ## null
    }
        
    dt[,eval(newCol):="None"]
    dt[get(cols[1]) == code['nosig'] & get(cols[2]) == code['sig'], eval(newCol):=lab[2] ]
    dt[get(cols[1]) == code['sig'] & get(cols[2]) == code['nosig'], eval(newCol):=lab[1] ]
    dt[get(cols[1]) == code['sig'] & get(cols[2]) == code['sig'], eval(newCol):="Both" ]
    dt[ , eval(newCol):=factor(get(newCol), levels=c("None", lab, "Both"))]

    return(dt)
}


#' Make genomicRanges object and data table with annotations to use ggbio to plot genes within a region
#'
#' @param gene.id gene id to center cis-window
#' @param cis_w numeric with size of cis window for each side og the gene in base pairs
#' @param gene.info data table with Gene_id start, end, chrom columns with these names
#' @param color character vector with colors for gene labels for all genes expect the reference one and the reference one
#' @keywords annotations gene plot
#' @export
#' @return list, first element genomicRanges object with gene info to plot and second; data table with annotations to add to plot
#' data.ggbio()

data.ggbio <- function(gene.id, cis_w, gene.info, color=NULL){

    gene <- add.name(data.table(Gene_id=gene.id))
    
    cis.w <- unlist(gene.info[gene_id ==  gene.id,  .(start, end)]) + c(-cis_w -10, cis_w+10) 

    chr <- gene.info[gene_id == gene.id, chrom]

    data(genesymbol, package = "biovizBase")
    
    wh <- range(genesymbol[seqnames(genesymbol) == paste0("chr", chr) & start(genesymbol)  > cis.w[1] & end(genesymbol) < cis.w[2]],
                ignore.strand = T)
    
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    gene_symbol<- data.table(AnnotationDbi::select(org.Hs.eg.db, keys=genes(txdb)$gene_id, 
                                                   
                                    columns="SYMBOL", keytype="ENTREZID"))

    temp <- crunch(txdb, which = wh)
    temp <- temp[temp$gene_id != ""]

    entrez <- data.table(ENTREZID=temp$gene_id)
    entrez <- merge(entrez, gene_symbol, by = "ENTREZID", all.x = T, sort=F)
              
    temp$symbol = entrez$SYMBOL

    tmp <- split(temp, temp$symbol)

    ## get start for first exon to place legend
    all.st <- as.data.table(temp)
    all.st <- all.st[type=="exon",][, .SD[1], symbol]
    if(!is.null(color)){
        
        all.st[, col:=color[1]][symbol==gene$Gene_name, col:=color[2]]
    }
    
    setkey(all.st, start)
    return(list(granges=tmp, annot=all.st))
}

#' Combine facet plot with ggbio plot
#'
#' @param p1 facet plot or list of facet plots
#' @param p2 ggbio plot
#' @param xlab labe for x axis
#' @param ylab label for y axis
#' @param title title for plot
#' @param rel.h vector with relative heights for facet and ggbio plots, defaults to 1,0.5
#' @keywords combine plots
#' @export
#' @return combined ggplot
#' plot.ggbio()

plot.ggbio <- function(p1, p2, xlab=NULL, ylab=NULL, title=NULL, rel.h=c(1,0.5)){
    ## remove titles from p1
    p <- lapply(p1, function(i) i + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                                          axis.text.x=element_blank(),
                                          plot.title = element_blank(),
                                          axis.ticks.x = element_blank()))
    
    p1 <- lapply(p, function(i) ggplot_gtable(ggplot_build(i)))
    ##p1 <- ggplot_gtable(ggplot_build(p))
    p2 <- ggplot_gtable(ggplot_build(p2))

   
    p2$widths <- p1[[1]]$widths

    if(!is.null(title)){
        title.grob <- textGrob(title, gp=gpar(fontsize=10))
    } else{
        title.grob <- NULL
    }
    

    if(!is.null(ylab)){

        y.grob <- textGrob(ylab, 
                           gp=gpar( fontsize=11), y = 0.65, rot=90)
    } else{
        y.grob <- NULL
    }
    

    if(!is.null(xlab)){

        x.grob <- textGrob(xlab,
                           gp=gpar( fontsize=11), y=1)
    } else {
        x.grob <- NULL
    }
    
    
    if(class(p1)=="list"){
        n <- length(p1)
        p1[[n+1]] <- p2
        plots <- plot_grid(plotlist=p1, ncol=1,
                           rel_heights=rel.h)
    } else {
        plots <- plot_grid(p1,p2, ncol=1,
                           rel_heights=rel.h)
    }
    
    

    a <- plot_grid(arrangeGrob(plots,
                               left = y.grob,
                               bottom = x.grob,
                               top=title.grob))

    return(a)
}

#' combine facet plot  with ggbio plot, wrapper for plot.ggbio and autoplot
#'
#' @param plot facet plot for 1 gene or list of plots all output from gene_plot_skin 
#' @param gene character with the gene name of plot
#' @param x character with gene id for gene
#' @param yaxis character with y-axis label, defaults to eQTL effect
#' @param y.text y coord to add gene names to autoplot
#' @param y.rect y coord for rectangule
#' @param just  hjustification for gene names in autoplot
#' @param geneStEnd, data table with gene_id, start, end, chrom 
#' @param d distance to plot defaults to 1e5
#' @param w width to add to rectangule to separate from gene, defaults to 2500 bp
#' @param relh vector with relative heights for facet and ggbio plots, defaults to 1,0.6
#' @keywords combine plots
#' @export
#' @return plot of combined ggplots
#' mix.ggbio()

mix.ggbio <- function(plot, gene, x, yaxis="eQTL effect", y.text, y.rect, just, geneStEnd, d=10^5, w=2500, relh=c(1,0.6)){
    data(genesymbol, package = "biovizBase")
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    print(gene)
    g.id <-x
    cis.w <- unlist(geneStEnd[gene_id == g.id, .(start, end)]) + c(-d -10, d+10)
    chrom <- geneStEnd[gene_id == g.id, chrom]
    wh <- range(genesymbol[seqnames(genesymbol) == paste0("chr", chrom)  & start(genesymbol)  > cis.w[1] & end(genesymbol) < cis.w[2]], ignore.strand = T)
    
    gene_symbol<- data.table(AnnotationDbi::select(org.Hs.eg.db, keys=genes(txdb)$gene_id,
                                                   columns="SYMBOL", keytype="ENTREZID"))

    temp <- crunch(txdb, which = wh)

    ## remove missing gene_id

    temp <- temp[temp$gene_id != ""]

    entrez <- data.table(ENTREZID=temp$gene_id)
    entrez <- merge(entrez, gene_symbol, by = "ENTREZID", all.x = T, sort=F)

    ## add gene symbol to temp, but some regions are not in gene_symbol (chr6), if it is the case try alternative way (I kept original way for backwards compatibility)
    if(any(is.na(entrez$SYMBOL))) {
        temp$symbol <-  unlist(lapply(as.character(temp$gene_id), function(i) lookUp(i, 'org.Hs.eg', 'SYMBOL')))
    } else {
        temp$symbol = entrez$SYMBOL
    }
    
    tmp <- split(temp, temp$symbol)

    ## get start for first exon to place legend
    all.st <- as.data.table(temp)
    all.st <- all.st[type=="exon",][, .SD[1], symbol]
    all.st[, col:="springgreen4"][symbol==gene, col:="black"]
    setkey(all.st, start)
    all.st[, y:= y.text]
    all.st[, hjust:= just]
    ##all.st[strand=="+", y:= c(0.6,1.5)][strand=="-", y:=c(2.5,1.8)]

    ## window for rectangle
    st.end <- unlist(geneStEnd[gene_id == g.id, .(start, end)])

    

    ## need to add gene names manually as otherwise they overlap or arent the right size

    dt <- data.table(x1=st.end[1]-w, x2=st.end[2] + w, y1=y.rect[1], y2=y.rect[2])
    
    if(nrow(all.st[is.na(symbol),])){
        return(plot)
    } else {
            
    p.txdb <- autoplot(tmp, aes(type=type), label=F, label.color = "black", color = "springgreen4",fill = "springgreen4")+ xlim(cis.w) + scale_x_sequnit("Mb") +
        ggplot2::geom_rect(data=dt, mapping=aes(xmin=x1 , xmax=x2, ymin=y1, ymax=y2), fill=NA, colour="gray") +
        ## geom_rect(mapping=aes(xmin= st.end[1]-w, xmax=st.end[2] + w, ymin=1.7, ymax=2.3), fill=NA, colour="gray") +
        theme(##legend.position = "none",
            ##panel.grid = element_blank(),
            ##axis.title = element_blank(),
            axis.text = element_text(size=10) ##,
            ##axis.ticks.x = element_blank())
        ) +
        
        annotate("text", label = all.st[["symbol"]], x=all.st[["start"]], y=all.st[["y"]] , colour=all.st[["col"]], hjust = all.st[['hjust']], size=2)

    p.gene <- plot.ggbio(p1=plot, p2=p.txdb@ggplot, ylab=yaxis, title=gene, rel.h=relh)

        return(p.gene)
    }
    

}
    
    
#' Plot to compare gene-snps associations by PIP and faceting by a categorical variable
#'
#' @param sum99  data table with stan output in long format
#' @param gene character vector with gene ID
#' @param pos column name to use for position for x axis, defaults to tag
#' @param geneStEnd data table  with gene_id, gene start and end columns, defaults to NULL
#' @param gwas whether to highlight the SNP in sum with the highest r2 with the gwas hit, when there is one
#' @param title if all it includes gene name plus source of gene (GWAS or fold change) plus median info plus min.p value for fisher test, otherwise just gene name
#' @param var categorical variable for facetting, defaults to Param
#' @param null name for null column, defaults to log2_aFC_null
#' @param colvar name of variable to use for color, defaults to Signif
#' @param colors character vector with colors to use, names factor in colvar, values colors, defaults to Signif
#' @param geneCoord to use with d for setting x lim, same data table as geneStEnd
#' @param d distance from gene to plot (gene start-d, gene end+d), defaults to NULL
#' @keywords stan plot gene-snps skin
#' @export
#' @return ggplot object
#' gene_plot_PIP()
gene_plot_PIP <- function(sum99, gene, pos="tag", geneStEnd=NULL, gwas=NULL, title="all", var="Param",  null=NULL, colvar="Signif", colors= setNames(c("steelblue4","gold3"), c("Yes","No")) , d=NULL,geneCoord=NULL){
    dt <- sum99[Gene_id==gene,]
   
    ## create col for position
    dt[, Position:=as.numeric(gsub(":.*","", get(pos)))/10^6]
           
    ## create col Signif instead of log2_aFC_null
    if(! "Signif" %in% names(dt)){
        if('log2_aFC_null' %in% names(dt)) {
            dt[, log2_aFC_null:=factor(log2_aFC_null, levels=c("yes", "no"))]
        } else {
            dt[, log2_aFC_null:=factor(get(null), levels=c("yes", "no"))]
        }
        dt[, Signif:=ifelse(log2_aFC_null=="no","Yes","No")]
    } else {
    
        ## change Signif to Yes, No if ti was already in dt
        dt[, Signif:=ifelse(Signif=="no","No","Yes")]

    }

    ## get cred int 
    ## ci.low <-paste0("log2_aFC_",ci[1], "%")
    ## ci.high <- paste0("log2_aFC_",ci[2], "%")
 
        
    ## dt[,ci.low:= get(ci.low)]
    ## dt[,ci.high:=get(ci.high)]
    ## set colors for consistency across plots
    ##colors <- c("steelblue4","gold3") #c("#D55E00","#009E73")
    ##names(colors) <- c("Yes","No")
    man.col <- colors[names(colors) %in% dt[[colvar]]]

    ## get gene name for ggtitle
    geneName <- unique(dt[!is.na(Gene_name), Gene_name])
  
    p <- ggplot(dt, aes(x=Position, y=post.out, shape=Top)) +
        
        scale_colour_manual(values=man.col) +
        geom_hline(yintercept=0) +
        #geom_errorbar(data=dt, aes(ymin=ci.low, ymax=ci.high), color="grey83",linetype="dashed",width=0.0) +
        geom_point(aes(colour=get(colvar))) +
        labs(color = colvar) +
        ylab(label="PIP") +
        xlab(label="Position (MB)") +
        
        theme_bw()

    if(!is.null(gwas)){
        ## gwas hit
        if(any(!is.na(dt$r2TagGwas))){
            dt[, r2:=""][r2TagGwas==max(r2TagGwas,na.rm=T), r2:= paste0("GWAS r2= ",as.character(round(r2TagGwas,2)))]
            ## only keep one max per var, in case of repeated values
            row.max <- lapply(levels(as.factor(dt[[var]])), function(i) dt[r2!="" & get(var)==i, which=T])
            if (any(lapply(row.max, length) >1)){
                rem <- unlist(lapply(row.max, function(i) {
                    if( length(i) > 1){
                        i[2:length(i)]
                    }
                }))
                dt <- dt[!rem,]
                }                
                
            p <- p + geom_text_repel(aes(label=r2), size=3 ,nudge_y=min(dt$log2_aFC_mean),
                                     direction="x",
                                     segment.color="red",
                                     data=dt) 
        }
    }
    if(!is.null(geneStEnd)) {
        stend <- unlist(geneStEnd[gene_id==gene, .(start,end)])/10^6
        p <- p + geom_vline(xintercept=stend, linetype="dotted", color="#009E73")
    }

    if(!(is.null(geneCoord) & is.null(d))){
        xl <-  (unlist(geneCoord[gene_id==gene, .(start,end)]) + c(-d,d)) /10^6
        p <- p + xlim(xl)
    }
    

    if(title=="all"){
        ## get source of gene
        sourceGene <- ifelse(all(is.na(unique(dt$gene_dist))), paste("DRG/FC",round(unique(dt[!is.na(FC),FC]),2)), paste("GWAS",unique(dt$rs_id) ))
        p <- p + ggtitle(paste(unique(dt[,get(var)]), geneName,"from", sourceGene, "\nmedian info=", round(median(dt$info),2),", min.p=",formatC(min(dt$min.p.fsnp),0)))
    } else {
        p <- p + ggtitle(geneName)
    }

   
    
    if(length(unique(dt[,get(var)])) > 1){
            
         p <- p + facet_grid(as.formula(paste(var, '~.'))) +
            theme(plot.title = element_text(size = 12),axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5), strip.background=element_rect(fill="white"))
    } else {
        p <- p + theme(plot.title = element_text(size = 12),axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5), strip.background=element_rect(fill="white"))
    }
    
  
    return(p)
    
}

#' get unique levels from 2 columns in data table to ease plotting
#'
#' @param cols vector with columns names
#' @param dt name of datatable
#' @export
#' @return vector with unique levels in any column
#' lev.cols()
#'
lev.cols <- function(cols, dt){
    return(unique(unlist(lapply(cols, function(i) levels(as.factor(dt[[i]]))))))
}


#' format dta table to ease plot for joint or independent models
#'
#' @param sum data table with stan sum
#' @param sig.col whether significant column is null.99 or Signif
#' @param col name of column to select top SNP, defaults to log2_aFC_d
#' @param model character with type of model run, either indpendent or joint
#' @param var character with name of column to make long table
#' @param rs if NULL there is a new col in sum named rsid with the top tag, otherwise a vector with names of files with rsids per chorm is provided (all or only relevant chrom, example:"/mrc-bsu/scratch/ev250/reference_genome/built37/variation_homo_sapiens-chr5.gvf.gz"
#' @param r list of matrices with correlation (or r2) between fSNPs, each element is for one gene
#' @param r.name name for r column, defaults to r
#' @export
#' @return data table to input into gene_plot_skin, new columns Top, rsid, r.col , r.col.cat and info.cat
#' for.plot()
#'
for.plot <- function(sum, sig.col=c("null.99", "Signif"), col="log2_aFC_d", var, rs=NULL, r , r.col="r"){

    ## if max(get(col)) is not unique for a gene/skin I need to select only one entry, I use gene.dist to select only one
    
    if (sig.col=="null.99"){
        tmp <- sum[null.99=="no",]
        
    } else {
        tmp <- sum[Signif=="yes",  ]
    }
    ## select tag with strongest effect size of more than one with max value -log10PEP
    tmp[, abs.mean:=abs(log2_aFC_mean)]
    setorderv(tmp, c(var, "Gene_id", "minuslog10PEP", "abs.mean"), order=c(1,1,-1,-1))
    
    top2add <- tmp[, .SD[1], by=c(var,"Gene_id")][, c(var, "Gene_id", "tag"), with=F]
    top2add[, Top:="x"]

    
    sum <- merge(sum, top2add, by=c(var, "Gene_id", "tag"), all.x=T)

    var.val <- unique(sum[[var]])
    ## add rsid/ tag to top
    u <- unique(sum[Top=="x", c("tag","chrom","Gene_id", var), with=F])
    if(!is.null(rs)){
        u.tags <- unique(u[, .(tag,chrom)])
        rs <- rbindlist(lapply(unique(u.tags$chrom), function(j)  get.rsid(grep(paste0("chr",j,".gvf.gz$"), rs, value=T),
                                                                  u.tags[chrom==j,tag])))
        u <- merge(u, rs, by.x="tag", by.y="SNP")
    } else {
        u[ ,rsid:=tag]
    }
    sum <- merge(sum, u, by=c("Gene_id", "tag", "chrom", var), all.x=T)
    ## when a combination of skin and gene doesnt have a sig snp, use the top tag for the other skin to get r2 with the snps tested
    u <- merge(u, expand.grid(v=var.val, Gene_id=unique(sum$Gene_id)),
               by.x=c("Gene_id", var),
               by.y=c("Gene_id","v"),
               all=T)

    
    ## convert Top to tag, 0 otherwise
    sum[Top=="x", Top:=tag][is.na(Top), Top:="0"]

    ## get gene_id for missing tags and then look for the tag for the same gene but other skin
    missing <- u[is.na(tag),]
    if(nrow(missing)){
        ## can be that both skins are not significant (using ba for top snp)
        m2.genes <- missing[,.N, Gene_id][N==2, Gene_id]
        if(length(m2.genes)){ ## ignore
            missing <- missing[!Gene_id %in% m2.genes,]
        }
        if(nrow(missing)){
            for(m in 1:nrow(missing)){
                top.other <- sum[Gene_id == missing[m,Gene_id] & get(var) != missing[m,get(var)] & Top !="0", Top]
                ##check if top.other was tested in non-sig skin:
                top.tested <- sum[Gene_id == missing[m,Gene_id] & get(var) == missing[m,get(var)], tag] == top.other
                if(any(top.tested)){
                    sum[Gene_id == missing[m,Gene_id] & get(var) == missing[m,get(var)] & tag==top.other, Top:=top.other]
                } else { ## if not tested select first tag as top
                    f.tag <- sum[Gene_id == missing[m,Gene_id] & get(var) == missing[m,get(var)], tag][1]
                    sum[Gene_id == missing[m,Gene_id] & get(var) == missing[m, get(var)] & tag==f.tag, Top:=top.other]
                }
                
            }
        }
        
    }
    
    ##add r/r2 based on top snp for genes that have at least one significant association in a skin type
    genes.ok <- as.character(unique(u[!is.na(tag),Gene_id]))
    sumr <- rbindlist(lapply(var.val,
                            function(i) rbindlist(lapply(genes.ok,
                                                         function(j) r2.stan(sum=sum[get(var) == i,], gene=j, r2=r[[j]], col=r.col)))))

    ## for failed genes rbind to sumr using r with 0
    if(exists("m2.genes")){
        sumr <- rbind(sumr, sum[Gene_id %in% m2.genes,][ , r:=1])
    }
    ## make r and info categorical by ranges

    sumr[ ,paste0(r.col, ".cat"):=cut(x=abs(get(r.col)), breaks=c(seq(0, 0.8, 0.2), 1.01)), Gene_id]
    
    sumr[, info.cat:=cut(info, breaks=c(0.3, .5, .8, 1.1, max(sumr$info))) ]


    ## format rsid

    sumr[is.na(rsid), rsid:=""]

    ## sort by Gene_id
    setkeyv(sumr, c("Gene_id", var))

    return(sumr)


}
