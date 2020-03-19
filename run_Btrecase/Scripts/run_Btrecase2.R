## log <- file(snakemake@log[[1]], open="wt")
## sink(log)
## sink(log, type="message")


source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/Btrecase.R")


btrecase.nogt.rna(gene=snakemake@params[['gene']],
                  chr=as.numeric(snakemake@params[['chrom']]),
                  snps=as.numeric(snakemake@params[['snps']]),
                  counts.f=snakemake@input[['counts']],
                  covariates=snakemake@input[['libsize']],
                  e.snps=snakemake@input[['eSNPs']],
                  gene.coord=snakemake@input[['genecoord']],
                  vcf=snakemake@input[['vcf']],
                  le.file=snakemake@input[['leRef']],
                  h.file=snakemake@input[['hapRef']],
                  population=snakemake@params[['pop']],
                  maf=as.numeric(snakemake@params[['maf']]),
                  min.ase=as.numeric(snakemake@params[['minAse']]),
                  min.ase.snp=as.numeric(snakemake@params[['minAseSnp']]),
                  min.ase.n=as.numeric(snakemake@params[['minAseN']]),
                  tag.threshold=as.numeric(snakemake@params[['tag']]),
                  q.test="no",
                  info=as.numeric(snakemake@params[['info']]),
                  out=snakemake@params[['out']],
                  prefix=snakemake@params[['prefix']],
                  model="/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.eff2.stan",
                  ex.fsnp=as.numeric(snakemake@params[['pfsnp']])
                  )


## e.snps='/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/fSNP/chr1.fSNP.genes.txt'
## #gene="ENSG00000233645"
## gene="ENSG00000000457"
## counts.f <-  "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts/Psoriasis_skin.txt"
## covariates <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts/Psoriasis_skin_gc_lib_size.rds"
## vcf <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/GT/chr1.ASE.Psoriasis_skin.vcf.gz"
## chr <- 1
## gene.coord <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt"
## le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz'

## h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz'

## model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.eff2.stan'

## population="EUR"
## maf=0.05
## nhets=5
## min.ase=5
## min.ase.snp=5
## min.ase.n=5
## tag.threshold=.9
## q.test="no"
## info=0.3
## snps=5*10^5
## prob=NULL
## ex.fsnp=NULL
