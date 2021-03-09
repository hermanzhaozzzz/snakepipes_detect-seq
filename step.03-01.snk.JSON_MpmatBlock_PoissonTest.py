TREAT = [
    "YE1-EMX1-PD",
    "YE1-VEGFA-PD",
    "33A-EMX1-PD",
    "33A-VEGFA-PD",
    "M2-EMX1-PD",
    "M2-VEGFA-PD",
    "All-EMX1-PD",
    "All-VEGFA-PD",
#     "test",
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
BEDTOOLS = os.path.join(PATH, 'bedtools')  # must be version 2.27.0 !
SAMTOOLS = os.path.join(PATH, 'samtools')
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/fasta_hg38/hg38_only_chromosome.fa.fai"



# filter
    # SiteMutNum: Detect region的C2T count 
    # SiteMutRatio: Detect region的C2T ratio 
    # SiteCoverNum: Detect region的coverage reads count 
    # RegionPassNum: Detect 
    # RegionToleranceNum: 
    
    

# 这里填step.02的heatmap选好的cutoff来跑第三步！
# SiteMutNum = [4]                 # M:  [影响较大]   建议多跑几个range, 1~10
# SiteCoverNum = [10]              # C:  [影响较大]   不要太小就行，一般设置10，要call点准确的话可以加高
# RegionPassNum = [1]              # R:  [影响较大]   region的pass的CtoT的个数，越高越严格，但是不能很大，2-3个就可以了，不要超过10
# SiteMutRatio = [0.1]             #  :  [影响较小]   影响不如count大，先考虑count再考虑ratio，coverage非常高的时候可以试试卡ratio，call不出点可以小点，比如0.02，正常CBE就0.1
# RegionToleranceNum = [2]         #  :  [影响较大]   Detect region可以容忍的false的C的个数，越低越严格, 影响比较大，探索时，可以设为False，False时取值为10000, 建议设置成2 
# explore

# c4 m3 r1
SiteMutNum = [3]         
SiteCoverNum = [4]       
RegionPassNum = [1]      
SiteMutRatio = [0.1]     
RegionToleranceNum = [2] 


rule all:
    input:
        ############################
        # block them if use json
        expand("../mpileup_pmat_bmat/293T-bat_{treat}_{rep}_hg38.MAPQ20.bmat", treat=TREAT, rep=REP),
        expand("../mpileup_pmat_bmat/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.bmat", ctrl=CTRL, rep=REP),
        # block them if use bmat
#         expand("../Split_Bmat_Json/293T-bat_{treat}_{rep}_hg38.MAPQ20.json", treat=TREAT, rep=REP),
#         expand("../Split_Bmat_Json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.json", ctrl=TREAT, rep=REP),
        ############################
        # the aim mpmat regions!
        ############################
        # mpmat
        # input for mpmat_block
        expand("../mpmat_filter_merge/293T-bat_{treat}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
               treat=TREAT, rep=REP, a=MPMAT_MERGE[0], b=MPMAT_MERGE[1], c=MPMAT_MERGE[2], d=MPMAT_MERGE[3],
               SMN=SiteMutNum, SCN=SiteCoverNum, SMR=SiteMutRatio, RPN=RegionPassNum, RTN=RegionToleranceNum),
        ############################
        # output for mpmat_block
#         expand("../mpmat_with_block_info/293T-bat_TREAT-{treat}_{rep}_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info.mpmat",
#                ctrl=CTRL, treat=TREAT, rep=REP, a=MPMAT_MERGE[0], b=MPMAT_MERGE[1], c=MPMAT_MERGE[2], d=MPMAT_MERGE[3],
#                SMN=SiteMutNum, SCN=SiteCoverNum, SMR=SiteMutRatio, RPN=RegionPassNum, RTN=RegionToleranceNum),
#         expand("../mpmat_with_block_info/293T-bat_TREAT-{treat}_{rep}_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info_rmHeader.mpmat",
#                ctrl=CTRL, treat=TREAT, rep=REP, a=MPMAT_MERGE[0], b=MPMAT_MERGE[1], c=MPMAT_MERGE[2], d=MPMAT_MERGE[3],
#                SMN=SiteMutNum, SCN=SiteCoverNum, SMR=SiteMutRatio, RPN=RegionPassNum, RTN=RegionToleranceNum),
#         ############################
#         # output for mpmat_merge
#         expand("../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.mpmat",
#                ctrl=CTRL, treat=TREAT, rep=REP, a=MPMAT_MERGE[0], b=MPMAT_MERGE[1], c=MPMAT_MERGE[2], d=MPMAT_MERGE[3],
#                SMN=SiteMutNum, SCN=SiteCoverNum, SMR=SiteMutRatio, RPN=RegionPassNum, RTN=RegionToleranceNum),
#         # output for mpmat_sort
#         expand("../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
#                ctrl=CTRL, treat=TREAT, rep=REP, a=MPMAT_MERGE[0], b=MPMAT_MERGE[1], c=MPMAT_MERGE[2], d=MPMAT_MERGE[3],
#                SMN=SiteMutNum, SCN=SiteCoverNum, SMR=SiteMutRatio, RPN=RegionPassNum, RTN=RegionToleranceNum),
#         # output for mpmat_merge2
#         expand("../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.rm-overlap.mpmat",
#                ctrl=CTRL, treat=TREAT, rep=REP, a=MPMAT_MERGE[0], b=MPMAT_MERGE[1], c=MPMAT_MERGE[2], d=MPMAT_MERGE[3],
#                SMN=SiteMutNum, SCN=SiteCoverNum, SMR=SiteMutRatio, RPN=RegionPassNum, RTN=RegionToleranceNum),
        ############################
        # poisson_test
        # input bam
        expand("../bam/293T-bat_{ctrl}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam", 
               ctrl=CTRL, rep=REP),
        expand("../bam/293T-bat_{treat}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam", 
               treat=TREAT, rep=REP),
        # tsv
        expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.tsv",
               ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
#         expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_call_peaks.tsv",
#        ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
#        SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.tsv",
               ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        expand("../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_noHeader.Sign.tsv",
               ctrl=CTRL,treat=TREAT,rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
#         Sign mpmat
        expand("../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.mpmat",
               treat=TREAT, ctrl=CTRL, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        #
#         expand("../mpmat_merge_Sign/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.art",
#                treat=TREAT, ctrl=CTRL, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
#                SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
#         # 
#         expand("../plot_art/293T-bat_TREAT-{treat}_CTRL-{ctrl}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.Sign_poisson.art_plot.pdf",
#                treat=TREAT, ctrl=CTRL, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
#                SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum)




# block if use json
rule input_bmat____mpmat_block:
    input:
        mpmat = "../mpmat_filter_merge/293T-bat_{treat}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
        ctrl_bmat = "../mpileup_pmat_bmat/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.bmat",
        treat_bmat = "../mpileup_pmat_bmat/293T-bat_{treat}_{rep}_hg38.MAPQ20.bmat"
    output:
        "../mpmat_with_block_info/293T-bat_TREAT-{treat}_{rep}_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info.mpmat"
    shell:
        """
        {PYTHON2} ./program/mpmat-block-V02.1.py \
        -i {input.mpmat} \
        -o {output} \
        -C {input.ctrl_bmat} \
        -T {input.treat_bmat} \
        -r {GENOME} \
        -p 24 \
        --query_mutation_type CT,GA
        """
rule rm_block_mpmat_header:
    input:
        "../mpmat_with_block_info/293T-bat_TREAT-{treat}_{rep}_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info.mpmat"
    output:
        "../mpmat_with_block_info/293T-bat_TREAT-{treat}_{rep}_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info_rmHeader.mpmat"
    shell:
        """
        sed '1d' {input} > {output}
        """

# block if use bmat
# rule input_json____mpmat_block:
#     input:
#         mpmat = "../mpmat_filter_merge/293T-bat_{treat}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
#         ctrl_json = "../Split_Bmat_Json/293T-bat_{ctrl}_{rep}_hg38.MAPQ20.json",
#         treat_json = "../Split_Bmat_Json/293T-bat_{treat}_{rep}_hg38.MAPQ20.json"
#     output:
#         "../mpmat_with_block_info/293T-bat_TREAT-{treat}_{rep}_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info.mpmat"
#     shell:
#         """
#         {PYTHON2} ./program/mpmat-block-V02.1.py \
#         -i {input.mpmat} \
#         -o {output} \
#         --ctrl_bmat_split_json {input.ctrl_json} \
#         --treat_bmat_split_json {input.treat_json} \
#         -r {GENOME} \
#         -p 24 \
#         --query_mutation_type CT,GA
#         """


# mpmat merge
# rule mpmat_merge:
#     input:
#         "../mpmat_with_block_info/293T-bat_TREAT-{treat}_rep1_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info_rmHeader.mpmat",
#         "../mpmat_with_block_info/293T-bat_TREAT-{treat}_rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info_rmHeader.mpmat"
#     output:
#         "../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.mpmat"
#     shell:
#         """
#         cat {input[0]} {input[1]} > {output}
#         """
# rule mpmat_sort:
#     input:
#         "../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.mpmat"
#     output:
#         "../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat"
#     shell:
#         """
#         {BEDTOOLS} sort -i {input} > {output}
#         """
# rule mpmat_merge2:
#     input:
#         "../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat"
#     output:
#         "../mpmat_filter_merge_rm-overlap/293T-bat_TREAT-{treat}_rep1-and-rep2_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.rm-overlap.mpmat"    
#     shell:
#         """
#         {BEDTOOLS} merge -i {input} > {output}
#         """



# next step
rule poisson_test:
    input:
        blocked_mpmat = "../mpmat_filter_merge/293T-bat_{treat}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
        ctrl_bam = "../bam/293T-bat_{ctrl}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam",
        treat_bam = "../bam/293T-bat_{treat}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam"
    output:
        "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.tsv"
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
        --mpmat_filter_info_col_index 13 \
        --mpmat_block_info_col_index -1 \
        --seq_reads_length 150 \
        --poisson_method mutation
        """
        
# usually block the rule below unless you need it to call peaks
# rule call_peaks_macs2:
#     input:
#         blocked_mpmat = "../mpmat_with_block_info/293T-bat_TREAT-{treat}_{rep}_CTRL-{ctrl}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.block_info.mpmat",
#         ctrl_bam = "../bam/293T-bat_{ctrl}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam",
#         treat_bam = "../bam/293T-bat_{treat}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam"
#     output:
#         "../table/detect_seq.StatsTest.table_CTRL-{ctrl}_TREAT-{treat}_{rep}.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}_call_peaks.tsv"
#     shell:
#         """
#         {PYTHON2} ./program/find-significant-mpmat-V05.py \
#         -i {input.blocked_mpmat} \
#         -o {output} \
#         -c {input.ctrl_bam} \
#         -t {input.treat_bam} \
#         -r {GENOME} \
#         -p 24 \
#         --query_mutation_type CT,GA \
#         --mpmat_filter_info_col_index 13 \
#         --mpmat_block_info_col_index 14 \
#         --seq_reads_length 150 \
#         --poisson_method all
#         """
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
