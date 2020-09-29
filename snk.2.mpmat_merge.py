SAMPLES = [
    'VEGFA-All-PD',
    'VEGFA-Vector-PD',
    'VEGFA-sgRNADel-PD'
]
REP = [
    'rep1',
    'rep2'
]


MPMAT_MERGE = [['C'],['T'],['G'],['A']]  # mpmat_merge list means C2T.mpmat + G2A.mpmat == CT_merge_GA.mpmat
# # MPMAT_MERGE = [['A'],['G'],['T'],['C']]  # mpmat_merge list means A2G.mpmat + T2C.mpmat == AG_merge_TC.mpmat





PYTHON2 = "/home/zhaohuanan/miniconda3/envs/snakepipes_detect-seq_fastq_bam_plot/bin/python"
BEDTOOLS = "/home/zhaohuanan/miniconda3/envs/snakepipes_detect-seq_fastq_bam_plot/bin/bedtools"
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa.fai"


rule all:
    input:
        expand("../mpmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{a}2{b}.mpmat",sample=SAMPLES, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1]),
        expand("../mpmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{c}2{d}.mpmat",sample=SAMPLES, rep=REP,c=MPMAT_MERGE[2],d=MPMAT_MERGE[3]),
        expand("../mpmat_merge/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{a}2{b}_merge_{c}2{d}.mpmat",sample=SAMPLES, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3]),
        expand("../mpmat_merge/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{a}2{b}_merge_{c}2{d}.sort.mpmat",sample=SAMPLES, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3]),
        # todo
        expand("../bam/293T-bat_{sample}_{rep}_hg38.MAPQ20.bam", sample=SAMPLES, rep=REP),
        expand("../json/293T-bat_{sample}_{rep}_hg38.MAPQ20.background.json", sample=SAMPLES, rep=REP)
rule mpmat_merge:
    input:
        "../mpmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{a}2{b}.mpmat",
        "../mpmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{c}2{d}.mpmat"
    output:
        "../mpmat_merge/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{a}2{b}_merge_{c}2{d}.mpmat"
    shell:
        """
        cat {input[0]} {input[1]} > {output}
        """
# todo
rule sort_mpmat:
    input:
        "../mpmat_merge/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{a}2{b}_merge_{c}2{d}.mpmat"
    output:
        "../mpmat_merge/293T-bat_{sample}_{rep}_hg38.MAPQ20.select.{a}2{b}_merge_{c}2{d}.sort.mpmat"
    shell:
        """
        {BEDTOOLS} sort -i {input} -g {GENOME_FAI} > {output}
        """   
rule count_genome_background:
    input:
        "../bam/293T-bat_{sample}_{rep}_hg38.MAPQ20.bam"
    output:
        "../json/293T-bat_{sample}_{rep}_hg38.MAPQ20.background.json"
    log:
        "../json/293T-bat_{sample}_{rep}_hg38.MAPQ20.background.json.log"
    shell:
        """
        {PYTHON2} ./program/calculate-mut-stats-V02.py \
        -i {input} \
        -r {GENOME} \
        -p 24 \
        -o {output} 2>{log}
        """