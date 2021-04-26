TREAT = [
    "YE1-EMX1-PD",
    "YE1-VEGFA-PD",
    "33A-EMX1-PD",
    "33A-VEGFA-PD",
    "M2-EMX1-PD",
    "M2-VEGFA-PD",
    "All-EMX1-PD",
    "All-VEGFA-PD",
]

CTRL = [
    "Vector-VEGFA-PD",
]

REP = ["rep1","rep2"]

# mpmat_merge list means C2T.mpmat + G2A.mpmat == CT_merge_GA.mpmat
MPMAT_MERGE = [['C'],['T'],['G'],['A']]

import os
PATH = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_detect-seq/bin/"
PYTHON2 = os.path.join(PATH, 'python')
BEDTOOLS = os.path.join(PATH, 'bedtools')
SAMTOOLS = os.path.join(PATH, 'samtools')
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.fai"


SiteMutNum = [3]                 # M:  [影响较大]   建议多跑几个range, 1~10
SiteCoverNum = [6]              # C:  [影响较大]   不要太小就行，一般设置10，要call点准确的话可以加高
RegionPassNum = [1]              # R:  [影响较大]   region的pass的CtoT的个数，越高越严格，但是不能很大，2-3个就可以了，不要超过10
SiteMutRatio = [0.1]             #  :  [影响较小]   影响不如count大，先考虑count再考虑ratio，coverage非常高的时候可以试试卡ratio，call不出点可以小点，比如0.02，正常CBE就0.1
RegionToleranceNum = [2]         #  :  [影响较大]   Detect region可以容忍的false的C的个数，越低越严格, 影响比较大，探索时，可以设为False，False时取值为10000, 建议设置成2 




rule all:
    input:
        # mpmat_filter_merge faked
#         expand("../mpmat_merge_after_mpmat_filtered_generated_by_faker/293T-{treat}_{rep}_hg38.mpmat", treat=TREAT, rep=REP),
        # sort_mpmat
        expand("../mpmat_merge_after_mpmat_filtered_generated_by_faker_sort/293T-{treat}_{rep}_hg38.sort.mpmat", treat=TREAT, rep=REP),
        # input bam
        expand("../bam/293T-bat_{ctrl}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam", ctrl=CTRL, rep=REP),
        expand("../bam/293T-bat_{treat}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam", treat=TREAT, rep=REP),
        expand("../table_faker/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}_C6M3R1_SMR0.1RT2.tsv", ctrl=CTRL,treat=TREAT,rep=REP),
# rule sort_mpmat:
#     input:
#         "../mpmat_merge_after_mpmat_filtered_generated_by_faker/293T-{sample}_{rep}_hg38.mpmat"
#     output:
#         "../mpmat_merge_after_mpmat_filtered_generated_by_faker_sort/293T-{sample}_{rep}_hg38.sort.mpmat"

#     shell:
#         """
#         {BEDTOOLS} sort -i {input} -g {GENOME_FAI} > {output}
#         """
rule poisson_test:
    input:
        blocked_mpmat ="../mpmat_merge_after_mpmat_filtered_generated_by_faker_sort/293T-{treat}_{rep}_hg38.sort.mpmat",
        ctrl_bam = "../bam/293T-bat_{ctrl}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam",
        treat_bam = "../bam/293T-bat_{treat}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam"
    output:
        "../table_faker/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}_C6M3R1_SMR0.1RT2.tsv"
    log:
        "../table_faker/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}_C6M3R1_SMR0.1RT2.tsv.log"
    shell:
        """
        {PYTHON2} ./program/find-significant-mpmat-V05.py \
        -i {input.blocked_mpmat} \
        -o {output} \
        -c {input.ctrl_bam} \
        -t {input.treat_bam} \
        -r {GENOME} \
        -p 24 \
        --query_mutation_type CT,GA \
        --mpmat_filter_info_col_index -1 \
        --mpmat_block_info_col_index -1 \
        --seq_reads_length 150 \
        --poisson_method mutation > {log} 2>&1
        """