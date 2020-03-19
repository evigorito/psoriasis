library(tictoc)



ufsnps=snakemake@input[['ueSNPS']]
if(!exists("ufsnps")){
    ufsnps <- NULL
}

AI_estimate=snakemake@input[['AI']]
if(!exists("AI_estimate")){
    AI_estimate <- NULL
}

prior=snakemake@params[['prior']]
if(!exists("prior")) {
    prior <- NULL
} else {
    k=length(prior)/3 ## number of gaussians
    s <- seq(1,length(prior),k)
    l <- lapply(1:3, function(i) as.numeric(prior[s[i]: (s[i]+k-1)]))
    names(l) <- c("mean", "sd", "mix")
    prior <- l
    
}


pretotalReads=snakemake@params[['pretotalReads']]
if(!exists("pretotalReads")){
    pretotalReads <- NULL
} else {
    pretotalReads <- as.numeric(pretotalReads)
}

mod=snakemake@input[['model']]



source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/Btrecase2T.R")

tic("Start running")

if(exists("mod")){
    btrecase2T.nogt.rna(gene=snakemake@params[['gene']],
                        chr=as.numeric(snakemake@params[['chrom']]),
                        snps=as.numeric(snakemake@params[['snps']]),
                        counts.f=snakemake@input[['counts']],
                        covariates=snakemake@input[['libsize']],
                        e.snps=snakemake@input[['eSNPs']],
                        u.esnps=ufsnps,
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
                        info=as.numeric(snakemake@params[['info']]),
                        out=snakemake@params[['out']],
                        model=mod,
                        prob=snakemake@params[['prob']],
                        prior=prior,
                        ex.fsnp=as.numeric(snakemake@params[['pfsnp']]),
                        AI_estimate=AI_estimate,
                        pretotalReads = pretotalReads,
                        skin=snakemake@params[['skin']]
                        )

} else {
    btrecase2T.nogt.rna(gene=snakemake@params[['gene']],
                        chr=as.numeric(snakemake@params[['chrom']]),
                        snps=as.numeric(snakemake@params[['snps']]),
                        counts.f=snakemake@input[['counts']],
                        covariates=snakemake@input[['libsize']],
                        e.snps=snakemake@input[['eSNPs']],
                        u.esnps=ufsnps,
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
                        info=as.numeric(snakemake@params[['info']]),
                        out=snakemake@params[['out']],
                        prob=snakemake@params[['prob']],
                        prior=prior,
                        ex.fsnp=as.numeric(snakemake@params[['pfsnp']]),
                        AI_estimate=AI_estimate,
                        pretotalReads = pretotalReads,
                        skin=snakemake@params[['skin']]
                        )
}

toc()
