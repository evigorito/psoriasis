##############################################################
## Apply bayesian trecase to psoriasis data, input preparation 
##############################################################

shell.prefix("source ~/.bashrc; ")

configfile: "config.yaml"

localrules: all

import os
import pandas as pd
from tempfile import TemporaryDirectory

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

def info_samples(s=config['sample_meta']):
    """
    Gets sample info including names and fastq paths from the meta data file in array express. 
    input: full file name
    output: data frame with sample meta data
    """
    data = pd.read_csv(s, sep="\t")
    return data

def dic_samples(var):
    """ Get dict with keys sample names from metadata and values any column from metadata.
    input: var is the variable name to use for dictionary values"""
    values = list(info_samples()[var])
    keys= list(info_samples()['Scan Name'])
    dic = dict(zip(keys,values))
    return dic
     




# 3 samples failed alignment ['SRR1146116','SRR1146118','SRR1146229'], I need to exclude them from rules to avoid recomputing star, get pandas data frame removing samples

pd_ex=info_samples()[~info_samples()['Comment [ENA_RUN]'].isin(['SRR1146116','SRR1146118','SRR1146229'])]
#pd_ex=info_samples()[info_samples()['Comment [ENA_RUN]'].isin(['SRR1146102'])]


rule all:
    input:
        #expand(config['output_dir'] + "/ASE/{sample}.chr{chrom}.ASE.vcf.gz", chrom=list(range(1,23)) , sample=pd_ex['Comment [ENA_RUN]'])
        #expand(config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fSNP.genes.txt", chrom=list(range(1,23)) )
        expand(config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fSNP.unique.genes.txt", chrom=list(range(1,23)) )


        
rule down:
    """ Downloads metadata for psoriasis RNAseq."""
    input:
        HTTP.remote("www.ebi.ac.uk/arrayexpress/files/E-GEOD-54456/E-GEOD-54456.sdrf.txt", keep_local=True)
    run:
        outputName =os.path.join(config['output_dir'] + "/sample_info", os.path.basename(input[0]))
        shell("mv {input} {outputName} ")


        
rule fastq_down:
    """ Downloads fastq files. Input is not directly used but is necessary to have the rule working
    """
    input:
        config['sample_meta']
    log:
        "logs/fastq_down/{sample}.log"
    params:
        dirout=config['output_dir'] + "/RNAseq/"       
    output:
        config['output_dir'] + "/RNAseq/{sample}.fastq.gz"     
    run:
        ftp= dic_samples('Comment[FASTQ_URI]')[wildcards.sample]
        shell("cd {params.dirout} ; "
              "wget {ftp} -A {wildcards.sample} -O {output}")
              

rule fastqc:
    """ Run fastqc. Since there can be race conditions if multiple jobs
    use the same fastqc dir, we create a temp dir."""
    input:
        config['output_dir'] + "/RNAseq/{sample}.fastq.gz"
    output:
        html=config['output_dir'] + "/FASTQC/{sample}.html",
        zipf=config['output_dir'] + "/FASTQC/{sample}.zip"
    params: ""
    run:       
        with TemporaryDirectory() as tempdir:
            html_path = os.path.join(tempdir, wildcards.sample + "_fastqc.html")
            zip_path = os.path.join(tempdir, wildcards.sample + "_fastqc.zip")
            shell(" fastqc {params} "
                  "--outdir {tempdir} {input} ")
            shell(" mv {html_path} {output.html} ;"
                  " mv {zip_path} {output.zipf} ")

rule filter_qual:
    """ Filter out reads with Ns and low quality (based on output from fastqc)"""
    input:
        config['output_dir'] + "/RNAseq/{sample}.fastq.gz"
    output:
         config['output_dir'] + "/QCfastq/{sample}.fastq.gz"
    params:
        MeanQual=25,
        Ns_per=1
    log:
        "logs/filter_qual/{sample}.log"
    shell:       
        "gzip -dc {input} | perl /mrc-bsu/scratch/ev250/bin/prinseq-lite.pl "
        " -fastq stdin -ns_max_p {params.Ns_per} "
        " -log {log} "
        " -out_format 3 -out_good stdout -out_bad null |  "
        " gzip > {output} "
        

rule star:
    input:
        config['output_dir'] + "/QCfastq/{sample}.fastq.gz"
    output:
        config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        index=config['indices'],
        read="zcat",
        out_dir=config['output_dir'] + "/STAR/{sample}/"
    threads: 8
    shell:
        "{config[STAR]} "
        " --runThreadN {threads} "
        " --genomeDir {params.index} "
        " --readFilesIn {input} "
        " --readFilesCommand {params.read} "
        " --outSAMtype BAM SortedByCoordinate "
        " --outFileNamePrefix {params.out_dir} "

       
rule total_gene_counts:
    """ Calculate total gene counts from RNA-seq BAM files, based on: counts for RNA-seq: http://www.bioconductor.org/help/workflows/rnaseqGene/#transcript-abundances-and-the-tximport-pipeline"""
    input:
        ebg=config['ebg'] ,
        bam=config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        mode="Union",
        ignore_strand="TRUE",
        singleEnd="TRUE"
    output:       
        config['output_dir'] + "/Btrecase/inputs/TotalGeneCounts/{skin}.{sample}.txt"
    script:
        "Scripts/total_gene_counts.R"

rule gene_info:
    """Get per gene info: start, end, chrom, longest transcript length and GC percentage """
    input:
        config['ref_gtf']
    output:
        config['output_dir'] + "/Btrecase/inputs/gene_inputs/gene_info.txt"
    script:
        "Scripts/gene_info.R"
        

rule group_gene_counts:
    """ Group total gene counts by skin type and make matrix with GC adjusted log(library size)"""
    input:
        counts=lambda wildcards: expand(config['output_dir'] + "/Btrecase/inputs/TotalGeneCounts/"+ wildcards.skin +".{sample}.txt",
                      sample=pd_ex[pd_ex['Comment [Sample_source_name]'].isin([wildcards.skin])]['Comment [ENA_RUN]']),
        geneInfo=config['output_dir'] + "/Btrecase/inputs/gene_inputs/gene_info.txt"
    output:   
        config['output_dir'] + "/Btrecase/inputs/Counts/{skin}.txt",
        config['output_dir'] + "/Btrecase/inputs/Counts/{skin}_gc_lib_size.rds"
    script:
         "Scripts/group_gene_counts.R"

rule list_bam:
    """ Prepare list of bam files to rename callVar vcf below. Format is full path to bam and sample name in each line space delimited. Also prepare bamPath, file with full name to bam files, input for calling variants"""
    input:
        bam=expand(config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam",
                   sample=pd_ex['Comment [ENA_RUN]'])
    output:
        bamSampNames=config['output_dir'] + "/STAR/bamSampleNames.txt",
        bamPath=config['output_dir'] + "/STAR/bamPath.txt"
    script:
        "Scripts/bamRename.py"

rule index:
    """ Index bam files """
    input:
        bam=config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index -b {input} {output}"
        
        
rule RNA_callVar:
    """Call variants from RNA, 4 shell commands to call variants, rename samples, add genomic annotations for QC, compress and index """
    input:
        bam=expand(config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam",
                   sample=pd_ex['Comment [ENA_RUN]']),
        bai=expand(config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai",
                   sample=pd_ex['Comment [ENA_RUN]']),
        fasta=config['fasta_ref'],
        bamList=config['output_dir'] + "/STAR/bamSampleNames.txt",
        bamPath=config['output_dir'] + "/STAR/bamPath.txt"
    output:
        samp=config['output_dir'] + "/call_vars/all/chr{chr}_Q20.vcf.gz",
        indx=config['output_dir'] + "/call_vars/all/chr{chr}_Q20.vcf.gz.tbi"
    params:
        q=255,
        QUAL=20
    threads: 16
    shell:
        "tmp=$(echo {config[output_dir]}/call_vars/all/chr{wildcards.chr}_Q20.vcf) ; "
        "tmp2=$(echo {config[output_dir]}/call_vars/all/chr{wildcards.chr}_Q20_v2.vcf) ; "
        "bcftools mpileup "
        "-r {wildcards.chr} -Ou "
        "-f {input.fasta} "
        "-b {input.bamPath} "
        "-q {params.q} "
        "-a DP -aDP -I "
        "--threads {threads} | "
        "bcftools call -mv -Ou --threads {threads} | "
        "bcftools view -i '%QUAL>={params.QUAL}' --threads {threads} -Ov -o $tmp ; "
        "bcftools reheader -s {input.bamList} $tmp > $tmp2; "
        "java -Xmx4g -jar {config[snpEFF]} "
        "-v GRCh37.75 -canon -ud 10000 -t -noStats "
        "$tmp2 >  $tmp ; "
        "bgzip -c $tmp > {output.samp} ; "
        "tabix -p vcf {output.samp} ; "
        

rule intersect_RP_psoriasis:
    """ Extracts and write records from input[0] shared by both input[0] and input[1] using exact allele match. In this case we extract from the reference panel the variants that are present in the psoriasis data"""
    input:
        RP=config['RP_alt'] + "/RP_chr{chrom}_alt_added.bcf" ,
        RPind=config['RP_alt'] + "/RP_chr{chrom}_alt_added.bcf.csi",
        samp=config['output_dir'] + "/call_vars/all/chr{chrom}_Q20.vcf.gz",     
        sampInd=config['output_dir'] + "/call_vars/all/chr{chrom}_Q20.vcf.gz.tbi"
    output:
        filt=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.vcf.gz",        
        filtInd=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.vcf.gz.tbi"
    shell:
        "bcftools isec -n=2 -w1 {input.samp} {input.RP} -Oz -o {output.filt} ; "         
        "tabix -p vcf {output.filt} "


rule select_variants_DP_1:
    """ Prepare inputs to subset variants called above user selected depth threshold"""
    input:
        samp=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.vcf.gz"    
    output:
        header=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.header.txt",
        body=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.txt"
    shell:
         "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%DP\\t%ANN[\\t%GT\\t%DP]\\n' "
         "{input.samp} > {output.body} ; "
         "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%DP\\t%ANN[\\t%GT\\t%DP]\\n' "
         "{input.samp} -H | "
         "head -1 > {output.header} "


rule select_variants_DP_2:
    """ Subset variants called above user selected depth threshold"""
    input:
        header=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.header.txt",
        body=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.txt"
    params:
        DP=10   
    output:
        body=config['output_dir'] + "/call_vars/DP10/chr{chrom}_Q20_filtRP_GTonly.txt"
    script:
        "Scripts/selVarsDP.R"

        
rule shapeit_in:
    """Prepare inputs for shapeit: 1) Need to add vcf body from rule select_variants_DP_2 into header, this file will be use for following rule  """
    input:
        vcf=config['output_dir'] + "/call_vars/RPvar/chr{chrom}_Q20_filtRP.vcf.gz",
        body=config['output_dir'] + "/call_vars/DP10/chr{chrom}_Q20_filtRP_GTonly.txt"
    output:
        vcf=config['output_dir'] + "/shapeit_in/chr{chrom}.allsamples.vcf"
    shell:
         "bcftools view -h {input.vcf} | "
         "bcftools annotate -x INFO,^FORMAT/GT  > {output.vcf} ; "
         "cat {input.body} >> {output.vcf} ; "
         
rule shapeit_out:
    """Phase GT using shapeit, I need to make a tmp vcf w/o missing values for each sample, convert output to vcf  """
    input:
        vcf=config['output_dir'] + "/shapeit_in/chr{chrom}.allsamples.vcf",
        refMap=config['RPpath'] + "/genetic_map_chr{chrom}_combined_b37.txt",
        refhap=config['RPpath'] + "/1000GP_Phase3_chr{chrom}.hap.gz",
        refleg=config['RPpath'] + "/1000GP_Phase3_chr{chrom}.legend.gz",
        refsamp=config['RPpath'] + "/1000GP_Phase3.sample",
        group=config['shapeit_group']
    output:
        vcf=config['output_dir'] + "/shapeit/chr{chrom}.{sample}.phased.vcf.gz",
        vcfInd=config['output_dir'] + "/shapeit/chr{chrom}.{sample}.phased.vcf.gz.tbi"
    shell:
        "tmp=$(echo {config[output_dir]}/shapeit_in/chr{wildcards.chrom}.{wildcards.sample}.vcf) ;"
        "tmp2=$(echo {config[output_dir]}/shapeit_in/chr{wildcards.chrom}.{wildcards.sample}.phased); "
        "vcf=$(echo {config[output_dir]}/shapeit_in/chr{wildcards.chrom}.{wildcards.sample}.phased.vcf); "
        "bcftools view -I -s {wildcards.sample} {input.vcf} | "
        "grep  -e '#' -e '[0-1]/[0-1]' > $tmp ; "
        "shapeit --input-vcf $tmp "
        "-M {input.refMap} "
        "--input-ref {input.refhap} {input.refleg} {input.refsamp} "
        "--include-grp {input.group} "
        "-O $tmp2 -T 1 ; "
        "shapeit -convert "
        "--input-hap $tmp2 "
        "--output-vcf $vcf ; "
        "bgzip -c $vcf > {output.vcf} ; "
        "tabix -p vcf {output.vcf} ; "
        #"rm $tmp ;"
        #"rm $tmp2 ;"        
    
rule phaser:
    """Apply phaser to get haplotypic counts output """
    input:
        vcf=config['output_dir'] + "/shapeit/chr{chrom}.{sample}.phased.vcf.gz",
        bam=config['output_dir'] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam",
        blacklistBed=config['phaser_blacklist'],
        blackHapCounts=config['phaser_blackHapCounts']
    params:
        config['output_dir'] + "/phaser/chr{chrom}.{sample}"
    output:
        config['output_dir'] + "/phaser/chr{chrom}.{sample}.allele_config.txt"
    threads: 12
    shell:
        "module load bedtools2-2.26.0-gcc-5.4.0-tqj36mo ; "
        "module load python-2.7.13-gcc-5.4.0-yubmrmn ; "
        "source /mrc-bsu/scratch/ev250/bin/PYTHON2/bin/activate ; "
        "python {config[phaser]} --vcf {input.vcf} "
        "--bam {input.bam} "
        "--paired_end 0 "
        "--mapq 255 "
        "--baseq 10 "
        "--sample {wildcards.sample} "
        "--blacklist {input.blacklistBed} "
        "--haplo_count_blacklist {input.blackHapCounts} "
        "--threads {threads} "
        "--gw_phase_vcf 1 --o {params} "
        
rule ASE:
    """Correct haplotypic counts phaser output from double counting reads with two het variants"""
    input:
        hapCounts=config['output_dir'] + "/phaser/chr{chrom}.{sample}.haplotypic_counts.txt"
    params:
        sep="_",
        min_cov=0
    output:
        config['output_dir'] + "/ASE/chr{chrom}.{sample}.corrected_haplotypic_counts.txt"
    script:
        "Scripts/variant_ase.py"  


rule format_GT:
    """Extract phased GT to make Btrecase input"""
    input:
        GT=config['output_dir'] + "/phaser/chr{chrom}.{sample}.vcf.gz"
    output:
        tab=config['output_dir'] + "/ASE/{sample}.chr{chrom}.less.tab",
        head=config['output_dir'] + "/ASE/{sample}.chr{chrom}.ASE.vcf"
    shell:
        "bcftools annotate -x ^FORMAT/GT {input} | "
        "bcftools query -f "
        "'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t GT[\\t%GT]\\n' "
        "-o {output.tab} ; "
        "bcftools annotate -x ^FORMAT/GT {input} | "
        "bcftools view -h | "
        "sed -e '7i\\\
##FORMAT=<ID=AS,Number=2,Type=Integer,Description=\"Allele-specific expression counts from RNA-seq\">' "
        " > {output.head} "


rule GT_ASE:
    """Merge phased GT with ASE counts"""
    input:
        GT=config['output_dir'] + "/ASE/{sample}.chr{chrom}.less.tab",
        ASE=config['output_dir'] + "/ASE/chr{chrom}.{sample}.corrected_haplotypic_counts.txt"
    output:
        tab=config['output_dir'] + "/ASE/{sample}.chr{chrom}.for.AS.tab"
    script:
        "Scripts/GT_ASE.R"

rule GT_ASE_vcf:
    """Make vcf with GT:ASE field"""
    input:
        head=config['output_dir'] +"/ASE/{sample}.chr{chrom}.ASE.vcf" ,
        body=config['output_dir'] + "/ASE/{sample}.chr{chrom}.for.AS.tab"
    output:
        vcf=config['output_dir'] + "/ASE/{sample}.chr{chrom}.ASE.vcf.gz",
        idx=config['output_dir'] + "/ASE/{sample}.chr{chrom}.ASE.vcf.gz.tbi"
    shell:
        "cat {input.body} >> {input.head} ; "
        "bgzip {input.head} ; "
        "bcftools index -t {output.vcf} "

        
rule merge_vcf:
    """Merge vcf with GT:ASE field by skin type"""
    input:
        lambda wildcards: expand(config['output_dir'] + "/ASE/{sample}.chr" + wildcards.chrom + ".ASE.vcf.gz", sample=pd_ex[pd_ex['Comment [Sample_source_name]'].isin([wildcards.skin])]['Comment [ENA_RUN]'])
    output:
        vcf=config['output_dir'] + "/Btrecase/inputs/GT/chr{chrom}.ASE.{skin}.vcf.gz",
        idx=config['output_dir'] + "/Btrecase/inputs/GT/chr{chrom}.ASE.{skin}.vcf.gz.tbi"
    shell:
        "bcftools merge -m none {input}  "
        "-Oz -o {output.vcf} ; "
        "bcftools index -t {output.vcf} "
        
rule fSNPid:
    """Get fSNP coordinates per chromosome. Merge vcfs by skin and then extract snp info"""
    input:
        snp=lambda wildcards: expand(config['output_dir'] + "/Btrecase/inputs/GT/chr" + wildcards.chrom + ".ASE.{skin}.vcf.gz", skin=set(info_samples()['Comment [Sample_source_name]']))
    output:
        fsnps=temp(config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fsnpID.txt"),
        head=temp(config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fsnpID.header.txt")
    shell:
        "bcftools merge -m none {input.snp} -Ou |  "
        "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input.snp}  "
        "> {output.fsnps} ; "
        "bcftools merge -m none {input.snp} -Ou |  "
        "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input.snp} -H | "
        "head -1 > {output.head} "

rule fSNP_coord:
    """Get fSNP coordinates per gene"""
    input:
        ebg=config['ebg'],
        fsnps=config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fsnpID.txt",
        head=config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fsnpID.header.txt"
    output:
        fsnps=config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fSNP.genes.txt"
    script:
        "Scripts/fSNP_coord.R"

rule fSNP_unique:
    """Get unique fSNPs per gene when no strand info is available"""
    input:
        fsnps=config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fSNP.genes.txt"
    output:
        ufsnps=config['output_dir'] + "/Btrecase/inputs/fSNP/chr{chrom}.fSNP.unique.genes.txt"
    script:
        "Scripts/fSNP_unique.R"
           
        
# snakemake -j 100 -k --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  -n {cluster.n}   -t {cluster.time} --output {cluster.error} -J {cluster.job} "


#snakemake -p -n     # testing (printing) commands, can add --quiet
          
