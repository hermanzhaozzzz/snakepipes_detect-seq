# CTRL = [
#     'VEGFA-Vector-PD'
# ]
# TREAT = [
#     'VEGFA-All-PD',
#     'VEGFA-sgRNADel-PD'
# ]

CTRL = [
    "VEGFA-Vector-PD" # 这个是CBE！
]
TREAT = [
    "GBEmini-AP-RNF2-All-PD",
    "GBEmini-dU-RNF2-All-PD",
]


REP = [
    'rep1'
]


MPMAT_MERGE = [['C'],['G'],['G'],['C']]  # mpmat_merge list means C2T.mpmat + G2A.mpmat == CT_merge_GA.mpmat
# # MPMAT_MERGE = [['A'],['G'],['T'],['C']]  # mpmat_merge list means A2G.mpmat + T2C.mpmat == AG_merge_TC.mpmat



SiteMutNum = [1,2,3,4,5,6,7,8] # 师兄的M4 建议多跑几个range
SiteCoverNum = [5,10] # 师兄的C10
SiteMutRatio = [0.01,0.1] # 影响不大
RegionPassNum = [2] # 师兄的R2
RegionToleranceNum = [2] # 影响不大，师兄建议设置成2






defult_sgRNA_dict_for_plot = {
    "ABE":"GAATACTAAGCATAGACTCC", # 这里指的ABE-Site-7-On-Target
    'EMX1': "GAGTCCGAGCAGAAGAAGAAGGG", 
    'HEK2':'GAACACAAAGCATAGACTGC',
    'HEK3': "GGCCCAGACTGAGCACGTGATGG", 
    "HEK4":"GGCACTGCGGCTGGAGGTGGGGG", 
    "HEK5":'CTGGCCTGGGTCAATCCTTG',
    "MTND4P12":"TGCTAGTAACCACATTCTCCTGATCAAATATCACTCTCCTACTTACAGGA",
    "MTND5P11":"TAGCATTGGCAGGAATACCCTTCCTCACAGGTTTCTACTCCAAAGA",
    "MTND6P4":"TGACCCCCATGCCTCAGGATACTCCTCAATAGCCACCG",
    "PP2":"GGCACTCGGGGGCGAGAGGA",
    "PP6":"GGGGCTCAACATCGGAAGAG",
    "RNF2":"GTCATCTTAGTCATTACCTGAGG",
    'VEGFA': "GACCCCCTCCACCCCGCCTCCGG", 
}



import os
PATH = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_detect-seq/bin/"
PYTHON2 = os.path.join(PATH, 'python')
BEDTOOLS = os.path.join(PATH, 'bedtools')
SAMTOOLS = os.path.join(PATH, 'samtools')
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa.fai"
HG38_JSON = "./lib/hg38.json"


rule all:
    input:
        expand("../mpmat_filter_merge/293T-bat_{treat}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
               treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        # 
        expand("../bam/293T-bat_{ctrl}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam", 
               ctrl=CTRL, rep=REP),
        expand("../bam/293T-bat_{treat}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam", 
               treat=TREAT, rep=REP),
        # 
        expand("../json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.background.json", 
               ctrl=CTRL, rep=REP),
        expand("../json/293T-bat_{treat}_{rep}_hg38.MAPQ20.background.json", 
               treat=TREAT, rep=REP),
        # 
        expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.tsv",
               ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.tsv",
               ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.Sign.tsv",
               ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        # 
        expand("../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.mpmat",
               treat=TREAT, ctrl=CTRL, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
#         #
#         expand("../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.art",
#                treat=TREAT, ctrl=CTRL, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
#                SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
#         # 
#         expand("../plot_art/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.art_plot.pdf",
#                treat=TREAT, ctrl=CTRL, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
#                SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum)
        
rule poisson_test:
    input:
        merge_sort_mpmat = "../mpmat_filter_merge/293T-bat_{treat}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
        ctrl_sort_bam = "../bam/293T-bat_{ctrl}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam",
        treat_sort_bam = "../bam/293T-bat_{treat}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam",
        bg_ctrl_json = "../json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.background.json",
        bg_treat_json = "../json/293T-bat_{treat}_{rep}_hg38.MAPQ20.background.json"
    output:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.tsv"
    log:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.tsv.log"
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
rule del_first_row_of_table:
    input:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.tsv"
    output:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.tsv"
    shell:
        """
        sed '1d' {input} > {output}
        """
rule select_significant_region:
    input:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.tsv"
    output:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.Sign.tsv"
    params:
        """'{if($NF<0.05 && $13>1 && $10<3) print $0}'"""
#         """'{if($NF<0.05 && $13>1 && $10<3) print $0}'"""
    shell:
        """
        cat {input} | awk -F '\t' {params} > {output}
        """        
        

rule intersect_select_mpmat:
    input:
        "../mpmat_filter_merge/293T-bat_{treat}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.Sign.tsv"
    output:
        "../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.mpmat"
    shell:
        """
        if [[ `cat {input[1]} |wc -l` -eq 0 ]]; then 
        echo "table is empty" 
        echo "will touch a empty file" 
        touch {output}
        else
        echo "table is ok!"
        echo "run bedtools intersect"
        {BEDTOOLS} intersect -a {input[0]} -b {input[1]} -wa > {output}
        fi
        """
# rule sgRNA_alignment:
#     input:
#         "../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.mpmat"
#     output:
#         "../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.art"
#     params:
#         sgRNA_seq = lambda wildcards, input: defult_sgRNA_dict_for_plot[input[0].split("TREAT-")[1].split("_CTRL-")[0].split('-')[0]]
#     shell:
#         """
#         if [[ `cat {input} |wc -l` -eq 0 ]]; then 
#         echo "Sign mpmat is empty" 
#         echo "will touch a empty file" 
#         touch {output}
#         else
#         echo "Sign mpmat is ok!"
#         echo "run sgRNA_alignment/[mpmat to art]"
#         {PYTHON2} ./program/mpmat-to-art-V03-zhn.py -i {input} -r {GENOME} --sgRNA {params.sgRNA_seq} > {output}
#         fi
#         """
# rule art_plot:
#     input:
#         "../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.art"
#     output:
#         "../plot_art/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.art_plot.pdf"
#     params:
#         sgRNA_seq = lambda wildcards, input: defult_sgRNA_dict_for_plot[input[0].split("TREAT-")[1].split("_CTRL-")[0].split('-')[0]]
#     shell:
#         """
#         if [[ `cat {input} |wc -l` -eq 0 ]]; then 
#         echo "Sign art is empty" 
#         echo "will touch a empty file" 
#         touch {output}
#         echo "Formal output\nNo significant results" > {output}.empty.log
#         else
#         echo "Sign art is ok!"
#         echo "run art plot"
#         {PYTHON2} ./program/plot-art-V01.py \
#         -i {input} \
#         --sgRNA {params.sgRNA_seq} \
#         -o {output}
#         fi
#         """
        
        
        
        
