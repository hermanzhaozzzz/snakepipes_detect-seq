SAMPLES = [
    'VEGFA-All-PD',
    'VEGFA-Vector-PD',
    'VEGFA-sgRNADel-PD'
]
REP = [
    'rep1',
    'rep2'
]

PYTHON2 = "/home/zhaohuanan/miniconda3/envs/snakepipes_detect-seq_fastq_bam_plot/bin/python"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_detect-seq_fastq_bam_plot/bin/bedtools"
SAMTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_detect-seq_fastq_bam_plot/bin/samtools"
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa.fai"
DB_SNP_HSA = "/home/zhaohuanan/zhaohn_HD/2.database/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"



SELECT_BASES_FROM = [
    'A',
    'T',
    'C',
    'G'
]
SELECT_BASES_TO = [
    'A',
    'T',
    'C',
    'G'
]
rule all:
    input:
        expand("../bam/293T-bat_{sample}_{rep}_hg38.MAPQ20.bam",sample=SAMPLES,rep=REP),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup",sample=SAMPLES,rep=REP),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat",sample=SAMPLES,rep=REP),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.bmat",sample=SAMPLES,rep=REP,fbase=SELECT_BASES_FROM),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.pmat",sample=SAMPLES,rep=REP,fbase=SELECT_BASES_FROM),
        expand("../mpmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}2{tbase}.pmat",sample=SAMPLES,rep=REP,fbase=SELECT_BASES_FROM,tbase=SELECT_BASES_TO)
rule bam2mpileup:
    input:
        "../bam/293T-bat_{sample}_{rep}_hg38.MAPQ20.bam"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup"
    log:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup.log"
    shell:
        """
        {SAMTOOLS} mpileup -q 20 -Q 20 \
        --reference {GENOME} \
        -o {output} \
        {input} >{log} 2>&1
        """
rule mpileup2bmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat"
    log:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat.log"
    shell:
        """
        {PYTHON2} ./program/parse-mpileup-V04.py \
        -i {input} \
        -o {output} \
        -p 24 -n 0 >{log} 2>&1
        """
rule select_bmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.bmat"
    params:
        awk = """{print $0}""",
        fbase = "{fbase}"
    shell:
        """
        awk '$3 == "{params.fbase}" {params.awk}' {input} > {output} 
        """
rule select_bmat2pmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.bmat"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.pmat"
    log:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.pmat.log"
    shell:
        """
        {PYTHON2} ./program/bmat2pmat-V02.py \
        -i {input} \
        -o {output} \
        --InHeader False \
        --InLikeBED False \
        --OutHeader True >{log} 2>&1
        """
rule pmat_merge:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.pmat"
    output:
        "../mpmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}2{tbase}.pmat"
    log:
        "../mpmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}2{tbase}.pmat.log"
    params:
        FB = '{fbase}',
        TB = '{tbase}'
    shell:
        """
        FB={params.FB}
        TB={params.TB}
        if [ $FB = $TB ]; then
        echo "from base == to base, skip..." > {log}
        touch {output}
        else
        {PYTHON2} ./program/pmat-merge-V04.py \
        -f C -t T \
        -r hg38.fa \
        --OutHeader False \
        -i {input} \
        -o {output} \
        -d 50 -D 100 \
        --NoMutNumCutoff 2 \
        --OmitTandemNumCutoff 2 \
        --SNP {DB_SNP_HSA} >{log} 2>&1
        fi
        """

