## Load R packages and functions:
library(data.table)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)

source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R", verbose=FALSE)
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R', verbose=FALSE)
source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")


#' Figures for eQTL no genotypes paper: psoriasis vs normal skin
#'
#' @param sum full name to data table with stan output 
#' @param sum2T full name to data table with stan output for 2 tissues in one model
#' @param gene_coord full name to data table  with gene_id, gene start and end columns
#' @param gene_names character vector with the names of genes to focus
#' @param y list with y coordinates for gene labels, each element for each gene
#' @param w numeric to adjust the size of the window for making rectangle for gene, defaults to 2000
#' @param tagsdir path to dir with tags for simple model
#' @param out dir to save files
#' @keywords figures eQTL no-genotypes
#' @export
#' @return save ggplots to files
#' 
#' figures4paper()

figures4paper <- function(sum, sum2T, gene_coord, gene_names, y, w=2000, tagsdir, out){
    
    sum99 <- fread(sum)

    ## add colum for fSNP being tagged
    sum99[, fSNP:="yes"][is.na(fSNP.id), fSNP:="no"]

    ## sum99 long to wide
    skin=unique(sum99$skin)

    

    sum99wide <- Reduce(function(x,y) merge(x,y, by=c("Gene_id","tag", "tag.EAF", "Gene_name","CHROM","POS","Alleles","rs_id", "gene_dist","table","RankP","CaseMedian", "FC","tag.pos","r2TagGwas","chrom", "fSNP.id", "fSNP.EAF", "tag.fsnp.op"), suffixes=paste0(".",skin), all=F), lapply(skin, function(i) sum99[skin==i,]))

    sum99wide <- add.sig(sum99wide, type="null", cols=c("log2_aFC_null.Psoriasis_skin", "log2_aFC_null.normal_skin"), lab=c("Pso.", "Norm."), newCol="Signif")

    table <- sum99wide[, .N,Signif][order(N, decreasing =T),]
    names(table) <- c("Signif", "SNPs")

  

######################################################################################################
############### Figure normal vs psoriasis skin all associations ####################################

    norPso <- btrecase.plot(sum99wide[Signif!="None",],
                            x1=paste(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'log2_aFC_null'), skin[1], sep="."),
                            x2=paste(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'log2_aFC_null'), skin[2], sep="."),
                            xl="eQTL-effect psoriasis skin",
                            yl="eQTL-effect normal skin",
                            col=c("Psoriasis", "Normal"),
                            axis.title=10,
                            axis.text=8,
                            legend.title=10,
                            legend.text=8,
                            legend.symbol=2,
                            point.size=1,
                            title="Normal vs psoriasis skin",
                            title.size=11) + 
        annotation_custom(tableGrob(table, rows=NULL, theme=ttheme_minimal(base_size = 8,padding = unit(c(2, 1.5), "mm"))), xmin=-1.3, xmax=-.5, ymin=0.5, ymax=1) 

    ggsave(file=paste0(out,'/allNormalPso.png'), norPso, width=4, height=4, dpi=100)


#########################################################################################################
######################## individual genes #############################################################
           
    geneStEnd <- fread(gene_coord)
    
    ## get tags for  gene_names, tags are the same for normal or psoriasis skin
    names(gene_names) <- sapply(gene_names, function(i) unique(sum99wide[Gene_name == i, Gene_id]))

    ## tags.all <- comb.files(paste0(tagsdir,"/", names(gene_names)), pattern="normal_skin.noGT.tags.lookup.txt")
    tags.all <- rbindlist(lapply(names(gene_names), function(i) comb.files(tagsdir, pattern=paste0("^refbias.", i, ".normal_skin.noGT.tags.lookup.txt"))))

    ##tags.all <- comb.files(paste0('/mrc-bsu/scratch/ev250/psoriasis/Btrecase/output99','/', names(plots.sub)), pattern="normal_skin.noGT.tags.lookup.txt")

    sum99.all <-  merge(tags.all, sum99, by= c("Gene_id","tag"))

    ## print(names(sum99.all))
    
    plots.all <- lapply(names(gene_names), function(i) gene_plot_skin(sum99.all,gene=i, gwas="yes", pos="SNP", title="no", var="skin",ci=c(0.5,99.5), colvar="Signif"))

    names(plots.all) <- names(gene_names)


    ## get genes within cis-window of a specif gene 
    cis_w=10^5

    ## combine plot with ggbio autoplot for each gene
    comb <- list()
    ## annot <- list()
    for(i in seq_along(gene_names)){
        
               
        gene.bio <- data.ggbio(gene.id=names(gene_names)[i], cis_w, gene.info=geneStEnd, color=c("springgreen4", "black"))
               
        all.st <- gene.bio$annot
               
        cis.w <- unlist(geneStEnd[gene_id == unique(sum99wide[Gene_name == gene_names[i], Gene_id]), .(start, end)]) + c(-cis_w -10, cis_w+10)

               ## make basic autoplot
        p.txdb <- autoplot(gene.bio$granges, aes(type=type), label=F, color = "springgreen4",fill = "springgreen4")+ xlim(cis.w) + scale_x_sequnit("Mb")

        ## modify all.st to format autoplot
        all.st[, y:= y[[i]] ]
        all.st[, hjust:= 0]

        ## window for rectangle to enclose SERPINB3
        
        st.end <-cis.w - + c(-10^5 -10, 10^5+10)
        d=data.frame(x1=st.end[1]-w,
                     x2=st.end[2]+w,
                     y1=1.7,
                     y2=2.3)
      

        p.txdb <- p.txdb + geom_rect(data=d, mapping=aes( xmin= x1,
                                                 xmax=x2,
                                                 ymin=y1,
                                                 ymax=y2
                                                 ), fill=NA, colour="gray") +
            theme(axis.text.x=element_text(size=7)) +
            
            annotate("text", label = all.st[["symbol"]], x=all.st[["start"]], y=all.st[["y"]], colour=all.st[["col"]], hjust = all.st[['hjust']], size=2)

       
        

        ## Combine eQTL plot with gene, first remove titles from plots.sub to use

        p.gene <- plot.ggbio(p1=plots.all[[ names(gene_names)[i] ]], p2=p.txdb@ggplot, ylab="eQTL-effect", title=gene_names[i], rel.h=c(1,0.6))

        comb[[gene_names[i] ]]  <- p.gene
    }
    
   
    p.corr <- plot_grid(norPso, NULL, ncol=1, rel_heights=c(1, 0.6))
    
    fig3 <- plot_grid(p.corr, plotlist=comb , nrow=1 , labels = 'auto')

    ggsave(file=paste0(out, '/Figure3.png'), fig3, width=12,
           height=4)
}


figures4paper(sum=snakemake@input[['stan_sum']],
              sum2T=snakemake@input[['stan_sum2T']],
              gene_coord=snakemake@input[['gene_coord']],
              gene_names=snakemake@params[['genes2follow']],
              y=list(snakemake@params[['y_ERAP1']], snakemake@params[['y_SERPINB3']]),
              w=as.numeric(snakemake@params[['w']]),
              tagsdir=snakemake@params[['tags_dir']],
              out=snakemake@params[['out_dir']]
              )
