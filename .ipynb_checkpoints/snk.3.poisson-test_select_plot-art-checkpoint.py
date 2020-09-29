REP = [
    'rep1',
#     'rep2'
]

CTRL = [
    'VEGFA-Vector-PD'
]
TREAT = [
#     'VEGFA-All-PD',
    'VEGFA-sgRNADel-PD'
]


MPMAT_MERGE = [['C'],['T'],['G'],['A']]  # mpmat_merge list means C2T.mpmat + G2A.mpmat == CT_merge_GA.mpmat
# MPMAT_MERGE = [['A'],['G'],['T'],['C']]  # mpmat_merge list means A2G.mpmat + T2C.mpmat == AG_merge_TC.mpmat




PYTHON2 = "/home/zhaohuanan/miniconda3/envs/snakepipes_detect-seq_fastq_bam_plot/bin/python"
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa.fai"
HG38_JSON = "./program/hg38.json"

rule all:
    input:
        # in
        expand("../bam/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.bam", ctrl=CTRL, rep=REP),
        expand("../bam/293T-bat_{treat}_{rep}_hg38.MAPQ20.bam", treat=TREAT, rep=REP),
        expand("../json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.background.json", ctrl=CTRL, rep=REP),
        expand("../json/293T-bat_{treat}_{rep}_hg38.MAPQ20.background.json", treat=TREAT, rep=REP),
        # out
        expand("../table/293T-bat_ctrl-{ctrl}_treat-{treat}_{rep}_hg38.MAPQ20.detect_seq.StatsTest.table", ctrl=CTRL, treat=TREAT, rep=REP)

rule poisson_test:
    input:
        merge_sort_mpmat = expand("../mpmat_merge/293T-bat_{treat}_{rep}_hg38.MAPQ20.select.{a}2{b}_merge_{c}2{d}.sort.mpmat",treat=TREAT, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3])[0],
        ctrl_sort_bam = "../bam/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.bam",
        treat_sort_bam = "../bam/293T-bat_{treat}_{rep}_hg38.MAPQ20.bam",
        ctrl_json = "../json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.background.json",
        treat_json = "../json/293T-bat_{treat}_{rep}_hg38.MAPQ20.background.json",

        
    output:
        "../table/293T-bat_ctrl-{ctrl}_treat-{treat}_{rep}_hg38.MAPQ20.detect_seq.StatsTest.table"
    params:
        from_base_to_base =  expand('{a}{b},{c}{d}',a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3])
    shell:
        """
        {PYTHON2} ./program/find-significant-mpmat-V02.py \
        -i {input.merge_sort_mpmat} \
        -c {input.ctrl_sort_bam} \
        -t {input.treat_sort_bam} \
        -r {GENOME} \
        -m {input.ctrl_json} \
        -n {input.treat_json} \
        -g {HG38_JSON} \
        -o {output} \
        --region_mutation_min_cutoff 2 \
        --query_mutation_type {params.from_base_to_base} \
        --query_mutation_min_cutoff 2 \
        --query_mutation_max_cutoff 18 \
        --other_mutation_max_cutoff 12 \
        --total_mutation_max_cutoff 26
        """   
# 5. Select signicicant regions and run sgRNA alignment
# 5.1 Requirement
# FILE: detect_seq.StatsTest.table
# CMD: mpmat-to-art-V04.py
# 5.2 Select significant region
# Related to our experience, the regions with the following criterion can be considered as enriched off-targets region:

# FDR < 0.05;
# log2_FC_mut > 1;
# ctrl_mut_count < 3 OR ctrl_mut_count.norm < 3;
# After this select step, we can obtain detect_seq.StatsTest.Sign.table.


# 5.2 Select .mpmat file
# bedtools intersect -a detect_seq.merge.sort.mpmat -b detect_seq.StatsTest.Sign.table -wa > detect_seq.merge.sort.Sign.mpmat
# 5.3 Run sgRNA alignment
# You can run sgRNA alignment as the following command. Only one thing you should take care of, the --sgRNA sequence has to include PAM info.

# mpmat-to-art-V04.py \
# -i detect_seq.merge.sort.Sign.mpmat -r hg38.fa --sgRNA GACCCCCTCCACCCCGCCTCCGG > detect_seq.merge.sort.Sign.art



# 6. Plot sgRNA alignment results
# 6.1 Requirement
# FILE: detect_seq.merge.sort.Sign.art
# CMD: plot-art-V01.py
# 6.2 Run code
# plot-art-V01.py -i detect_seq.merge.sort.Sign.art --sgRNA GACCCCCTCCACCCCGCCTCCGG --out_figure_format png -o out.art.png


