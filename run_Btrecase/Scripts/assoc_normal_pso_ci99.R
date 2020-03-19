#' ---
#' title: Comparing eQTL in normal vs psoriatic skin
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule compare_btrecase_normal_pso from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Snakefile

## Load R packages and functions:
library(data.table)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)


source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R", verbose=FALSE)
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R', verbose=FALSE)
source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")

#' ## Summary of associations
#'
#' Get output for psoriatic and normal skin for 99% CI and min p-value for fisher test 10^-4
#'

sum99 <- fread(snakemake@input[['stan_sum']])

#' Number of genes and distribution of the number of tag snps per gene run by skin type
unique(sum99[,Gene_id, skin])[,.N,skin]
sum99[, .N, .(Gene_id, skin)][,.(Distribution = summary(N)), skin]

#' Number of genes with non-null associations and distribution of the number of tag snps with non-null associations

unique(sum99[log2_aFC_null== "no",Gene_id,skin])[,.N,skin]
sum99[log2_aFC_null== "no", .N, .(Gene_id,skin)][,.(Distribution= summary(N)), skin]

#' ## QC
#'
## add colum for fSNP being tagged
sum99[, fSNP:="yes"][is.na(fSNP.id), fSNP:="no"]

## add column to check if fSNP and tag in same or opposite direction to correct effect size
sum99[!is.na(fSNP.id) , tag.fsnp.op:="no"][(fSNP.EAF < 0.5 & tag.EAF > 0.5) | (fSNP.EAF > 0.5 & tag.EAF < 0.5), tag.fsnp.op:="yes"]

#' * Density of log2 fold change for all SNPs tested
ggplot(sum99, aes(x=log2_aFC_mean, color=skin)) +  geom_density() +
    scale_colour_manual(values=c("#999999", "#E69F00")) + theme_bw()

## same but only sig effects by fSNP being rSNP or not

## create log2_fSNP to express eQTL effect relative to alt allele of fSNP, if no fSNP, use tag

sum99[,log2_fSNP:=log2_aFC_mean][tag.fsnp.op=="yes", log2_fSNP:=-log2_aFC_mean]

ggplot(sum99[log2_aFC_null== "no",], aes(x=log2_fSNP, color=skin)) +  geom_density() +
    scale_colour_manual(values=c("#999999", "#E69F00")) + theme_bw() + geom_vline(xintercept=0, linetype="dashed", colour ="blue") +
    facet_wrap(~fSNP, labeller=labeller(fSNP =c(no="fSNP is not rSNP", yes="fSNP is rSNP" )))

ggsave(file="/mrc-bsu/scratch/ev250/psoriasis/Btrecase/results/tables_fig/psoSigAssocbyFsnp.png")


#' * Density of log2 fold change by gene

## For each gene get the mean of log2 allelic fold change
sum99.mean <- sum99[, .(Mean.log2FC=mean(log2_aFC_mean)), by=c("Gene_id","skin")]
ggplot(sum99.mean, aes(x=Mean.log2FC, color=skin)) +  geom_density() +
    scale_colour_manual(values=c("#999999", "#E69F00")) + theme_bw()

## same but only with sig effects

sum99.mean.sig <- sum99[log2_aFC_null== "no", .(Mean.log2FC=mean(log2_aFC_mean)), by=c("Gene_id","skin")]
ggplot(sum99.mean.sig, aes(x=Mean.log2FC, color=skin)) +  geom_density() +
    scale_colour_manual(values=c("#999999", "#E69F00")) + theme_bw() 


#' * Look at "flat lines"
#' * For each gene look at the variance on log2FC vs median(r2 betwen rsnps)

## Get median r2 by skin type for each gene and var log2FC
r2 <- readRDS(snakemake@input[['r2']])

med.r2 <- lapply(unique(sum99$skin), function(j) {
    dt <- sum99[skin==j,]
    tmp <- rbindlist(lapply(1:length(r2), function(i) {
        tags <- dt[Gene_id == names(r2)[i], tag]
        m=median(r2[[i]][upper.tri(r2[[i]][tags,tags])])
        var.logaFC=var(dt[Gene_id == names(r2)[i],log2_aFC_mean])
        dt <- data.table(Gene_id=names(r2)[i], Median.r2=m, var.log2aFC= var.logaFC)
        return(dt)
    }))
    ## remove NA, not all genes were run in both conditions
    tmp <- tmp[!is.na(Median.r2),]
    return(tmp)
    })
        
names(med.r2) <- unique(sum99$skin)              

med.r2 <- rbindlist(med.r2, idcol="skin")

## Plot

ggplot(med.r2, aes(log(var.log2aFC), Median.r2, color=skin)) + geom_point(shape=1) +
    geom_smooth(method='lm', se=FALSE) 


#' * Compare all associations Normal vs psoriasis skin

## sum99 long to wide
skin=unique(sum99$skin)

sum99wide <- Reduce(function(x,y) merge(x,y, by=c("Gene_id","tag", "tag.EAF", "Gene_name","CHROM","POS","Alleles","rs_id", "gene_dist","table","RankP","CaseMedian", "FC","tag.pos","r2TagGwas","chrom", "fSNP.id", "fSNP.EAF", "tag.fsnp.op"), suffixes=paste0(".",skin), all=F), lapply(skin, function(i) sum99[skin==i,]))

sum99wide <- add.sig(sum99wide, type="null", cols=c("log2_aFC_null.Psoriasis_skin", "log2_aFC_null.normal_skin"), lab=c("Pso.", "Norm."), newCol="Signif")

table <- sum99wide[, .N,Signif][order(N, decreasing =T),]
names(table) <- c("Signif", "SNPs")

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



#' Plot associations by gene
gene_coord <- snakemake@input[['gene_coord']]
#gene_coord <-"/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt"
geneStEnd <- fread(gene_coord)

plots <- lapply(unique(sum99$Gene_id), function(i) gene_plot_skin(sum99,gene=i,pos="tag", geneStEnd,
                                                                  gwas="yes", var="skin",ci=c(0.5,99.5)))

names(plots) <- unique(sum99$Gene_id)

## ggsave(file='Btrecase/results/tables_fig/ENSG00000170465PsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/SERPINB3PsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/GSTP1PsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/TIMELESSPsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/PI3PsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/IKBKEPsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/RASI1PPsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/FUT2PsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/RMI2PsoNor.png',width=4.26, height=3.69, units="in")
## ggsave(file='Btrecase/results/tables_fig/ERAP1PsoNor.png',width=4.26, height=3.69, units="in")

## poster
## ggsave(file='Btrecase/results/tables_fig/KRT6CPsoNor.png', width=4.26, height=2.65, units="in")



## pso99 <- sum99[skin=="Psoriasis_skin",]
## plots2 <- lapply(unique(pso99$Gene_id), function(i) gene_plot_skin(pso99,gene=i,geneStEnd,
##                                                                   gwas="yes", var="skin" ,ci=c(0.5,99.5)))

## names(plots2) <- unique(pso99$Gene_id)

##ggsave(file='Btrecase/results/tables_fig/ERAP1Pso.png')

## Example plot with Median.r2 == 0.65 in pso skin

plots[['ENSG00000074800']]

##Example plot with Median.r2 0.04 in normal skin
plots[['ENSG00000107201']]

## some interesting genes to follow up from GWAS side:
## 1- IKBKE, NFKB1Z and REL ( this one a bit unlucky with r2, sig) but the 3 are components of the NFkB pathway implicated in psoriasis.
## 2- ERAP1, process peptides for HLA classI pathway
## 3-CAST has been implicated in skin diseases
## 4-RASIP1 ras pathway always important for signalling, but not super convincing
## 5-FUT2 more convincing but unclear link to psoriasis


#' ## Plots per gene
#' * Blue dotted lines indicate gene boundaries
#' * min.p the minimun p value for the Fisher test of proportion of het fSNPs in sample and reference panel

## plots

## Follow up:
gene_names=c("FUT2", "ERAP1", "RASIP1", "IKBKE", "TIMELESS", "RMI2", "KRT6C", "SERPINB3", "GSTP1")

## for these genes look at the direction of effects for fsnps

sum99.sub <- sum99[Gene_name %in% gene_names & !is.na(fSNP.id),]

## Is the direction of effects between fSNP and tag the same? yes
sum99.sub[,.N,tag.fsnp.op]


plots.sub <-  lapply(unique(sum99[Gene_name %in% gene_names, Gene_id]), function(i) gene_plot_skin(sum99,gene=i, title="no", var="skin",ci=c(0.5,99.5), colvar="Signif"))

names(plots.sub) <- unique(sum99[Gene_name %in% gene_names, Gene_id])


#####################################################################################
## Make Figure 3 for paper, select ERAP1 and SERPINB3, expand tags and add gene plot.
#####################################################################################

## get tags for plots.sub genes, tags are the same for normal or psoriasis skin

tags.all <- comb.files(paste0(snakemake@params[['output_dir']],"/", names(plots.sub)), pattern="normal_skin.noGT.tags.lookup.txt")

##tags.all <- comb.files(paste0('/mrc-bsu/scratch/ev250/psoriasis/Btrecase/output99','/', names(plots.sub)), pattern="normal_skin.noGT.tags.lookup.txt")

sum99.all <-  merge(tags.all, sum99, by= c("Gene_id","tag"))


plots.all <- lapply(unique(sum99.all[Gene_name %in% gene_names, Gene_id]), function(i) gene_plot_skin(sum99.all,gene=i, gwas="yes", pos="SNP", title="no", var="skin",ci=c(0.5,99.5), colvar="Signif"))

names(plots.all) <- unique(sum99.all[Gene_name %in% gene_names, Gene_id])


## get genes within cis-window of a specif gene SERPINB3
cis_w=10^5
gene= "SERPINB3"

serp <- data.ggbio(gene.id=unique(sum99wide[Gene_name == gene, Gene_id]), cis_w, gene.info=geneStEnd, color=c("springgreen4", "black"))
all.st <- serp$annot
cis.w <- unlist(geneStEnd[gene_id == unique(sum99wide[Gene_name == gene, Gene_id]), .(start, end)]) + c(-cis_w -10, cis_w+10)

## make basic autoplot
p.txdb <- autoplot(serp$granges, aes(type=type), label=F, color = "springgreen4",fill = "springgreen4")+ xlim(cis.w) + scale_x_sequnit("Mb")

## modify all.st to format autoplot
all.st[, y:= c(0.6, 1.4, 0.6, 1.6, 0.6)]
all.st[, hjust:= c(0,0,0,0,0)]

## window for rectangle to enclose SERPINB3

st.end <-cis.w - + c(-10^5 -10, 10^5+10)

w <- 2000
#print(st.end)

## p.txdb <- p.txdb + geom_rect(mapping=aes(xmin= st.end[1]-w,
##                                                                                                                                                       xmax=st.end[2] + w,
##                                                                                                                                                       ymin=1.7,
##                                                                                                                                                       ymax=2.3
##                                                                                                                                                       ), fill=NA, colour="gray") +
##     theme(axis.text.x=element_text(size=10)) +
    
##     annotate("text", label = all.st[["symbol"]], x=all.st[["start"]], y=all.st[["y"]], colour=all.st[["col"]], hjust = all.st[['hjust']], size=2.5)

## ## Combine eQTL plot with gene, first remove titles from plots.sub to use

## p.serp <- plot.ggbio(p1=plots.all[[1]], p2=p.txdb@ggplot, ylab="eQTL-effect", title=gene, rel.h=c(1,0.6))

###################################################################################################

## gene="ERAP1"
## gene.id=unique(sum99wide[Gene_name == gene, Gene_id])

## erap <- data.ggbio(gene.id, cis_w, gene.info=geneStEnd, color=c("springgreen4", "black"))
## all.st <- erap$annot

## cis.w <- unlist(geneStEnd[gene_id == unique(sum99wide[Gene_name == gene, Gene_id]), .(start, end)]) + c(-cis_w -10, cis_w+10)

## ## make basic autoplot
## p.txdb <- autoplot(erap$granges, aes(type=type), label=F, color = "springgreen4",fill = "springgreen4")+ xlim(cis.w) + scale_x_sequnit("Mb")

## ## modify all.st to format autoplot
## all.st[, y:= c(0.6, 1.6)]
## all.st[, hjust:= c(0,0)]

## ## window for rectangle to enclose SERPINB3

## st.end <-cis.w - + c(-10^5 -10, 10^5+10)

## w <- 2000

## p.txdb <- p.txdb + geom_rect(mapping=aes(xmin= st.end[1]-w,
##                                                                                                                                                       xmax=st.end[2] + w,
##                                                                                                                                                       ymin=1.7,
##                                                                                                                                                       ymax=2.3
##                                                                                                                                                       ), fill=NA, colour="gray") +
##     theme(axis.text.x=element_text(size=10)) +
    
##     annotate("text", label = all.st[["symbol"]], x=all.st[["start"]], y=all.st[["y"]], colour=all.st[["col"]], hjust = all.st[['hjust']], size=2.5)

## ## Combine eQTL plot with gene, first remove titles from plots.sub to use

## p.erap <- plot.ggbio(p1=plots.all[[gene.id]], p2=p.txdb@ggplot, ylab="eQTL-effect", title=gene, rel.h=c(1,0.6))

## p.corr <- plot_grid(norPso, NULL, ncol=1, rel_heights=c(1, 0.6))
## fig3 <- plot_grid(p.corr, p.erap, p.serp, nrow=1 , labels = 'auto')

## ggsave(file='/mrc-bsu/scratch/ev250/psoriasis/Btrecase/results/tables_fig/Figure3.png', fig3, width=30,
##        height=10, dpi=100)
