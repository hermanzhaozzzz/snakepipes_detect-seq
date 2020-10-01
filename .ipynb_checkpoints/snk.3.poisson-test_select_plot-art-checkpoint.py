CTRL = [
    'VEGFA-Vector-PD'
]
TREAT = [
    'VEGFA-All-PD',
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
HG38_JSON = "./lib/hg38.json"


rule all:
    input:
        # merge_sort_mpmat (treat)
        expand("../mpmat_merge/293T-bat_{treat}_{rep}_hg38.filtered.{a}2{b}_merge_{c}2{d}.sort.mpmat",treat=TREAT, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3]),
        # ctrl_sort_bam
        expand("../bam/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.bam", ctrl=CTRL, rep=REP),
        # treat_sort_bam
        expand("../bam/293T-bat_{treat}_{rep}_hg38.MAPQ20.bam", treat=TREAT, rep=REP),
        # bg_ctrl_json
        expand("../json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.background.json", ctrl=CTRL, rep=REP),
        # bg_treat_json
        expand("../json/293T-bat_{treat}_{rep}_hg38.MAPQ20.background.json", treat=TREAT, rep=REP),
        # output table of poisson test
        expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered.{a}2{b}_merge_{c}2{d}.tsv",ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3])
rule poisson_test:
    input:
        merge_sort_mpmat = "../mpmat_merge/293T-bat_{treat}_{rep}_hg38.filtered.{a}2{b}_merge_{c}2{d}.sort.mpmat",
        ctrl_sort_bam = "../bam/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.bam",
        treat_sort_bam = "../bam/293T-bat_{treat}_{rep}_hg38.MAPQ20.bam",
        bg_ctrl_json = "../json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.background.json",
        bg_treat_json = "../json/293T-bat_{treat}_{rep}_hg38.MAPQ20.background.json"
    output:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered.{a}2{b}_merge_{c}2{d}.tsv"
    log:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered.{a}2{b}_merge_{c}2{d}.tsv.log"
    params:
        mut_type = "{a}{b},{c}{d}"
    shell:
        """
        {PYTHON2} ./program/find-significant-mpmat-V02.py \
        -i {input.merge_sort_mpmat} \
        -c {input.ctrl_sort_bam} \
        -t {input.treat_sort_bam} \
        -r {GENOME} \
        -m {input.bg_ctrl_json} \
        -n {input.bg_treat_json} \
        -g {HG38_JSON} \
        -o {output} \
        --region_mutation_min_cutoff 2 \
        --query_mutation_type {params.mut_type} \
        --query_mutation_min_cutoff 2 \
        --query_mutation_max_cutoff 18 \
        --other_mutation_max_cutoff 12 \
        --total_mutation_max_cutoff 26 >{log} 2>&1
        """