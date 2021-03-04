SAMPLES = [
    "YE1-EMX1-PD",
    "YE1-VEGFA-PD",
    "33A-EMX1-PD",
    "33A-VEGFA-PD",
    "M2-EMX1-PD",
    "M2-VEGFA-PD",
    "All-EMX1-PD",
    "All-VEGFA-PD",
    "Vector-VEGFA-PD"
#     "test",
]

REP = ["rep1","rep2"]




MPMAT_MERGE = [['C'],['T'],['G'],['A']]  # mpmat_merge list means C2T.mpmat + G2A.mpmat == CT_merge_GA.mpmat
# # MPMAT_MERGE = [['A'],['G'],['T'],['C']]  # mpmat_merge list means A2G.mpmat + T2C.mpmat == AG_merge_TC.mpmat





import os
PATH = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_detect-seq/bin/"
PYTHON2 = os.path.join(PATH, 'python')
BEDTOOLS = os.path.join(PATH, 'bedtools')
SAMTOOLS = os.path.join(PATH, 'samtools')
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.fai"
DB_SNP_HSA = "/home/zhaohuanan/zhaohn_HD/2.database/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"



# 尝试好的cutoff来跑第三步！


# solid
# SiteMutNum = [5]                 # M:  [影响较大]   建议多跑几个range, 1~10
# SiteCoverNum = [10]              # C:  [影响较大]   不要太小就行，一般设置10，要call点准确的话可以加高
# RegionPassNum = [2]              # R:  [影响较大]   region的pass的CtoT的个数，越高越严格，但是不能很大，2-3个就可以了，不要超过10
# SiteMutRatio = [0.1]             #  :  [影响较小]   影响不如count大，先考虑count再考虑ratio，coverage非常高的时候可以试试卡ratio，call不出点可以小点，比如0.02，正常CBE就0.1
# RegionToleranceNum = [2]         #  :  [影响较大]   Detect region可以容忍的false的C的个数，越低越严格, 影响比较大，探索时，可以设为False，False时取值为10000, 建议设置成2 



# explore
SiteMutNum = [1,2,3,4,5,6,7,8]         
SiteCoverNum = [4,5,6,7,8,9]       
RegionPassNum = [0,1,2,3,4,5]      
SiteMutRatio = [0.1]     
RegionToleranceNum = [2] 

rule all:
    input:
        # input
        expand("../mpmat/293T-bat_{sample}_{rep}_hg38.select.{a}2{b}.mpmat",
               sample=SAMPLES, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1]),
        expand("../mpmat/293T-bat_{sample}_{rep}_hg38.select.{c}2{d}.mpmat",
               sample=SAMPLES, rep=REP,c=MPMAT_MERGE[2],d=MPMAT_MERGE[3]),
        # mpmat_filter
        expand("../mpmat_filter/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}.mpmat",
               sample=SAMPLES, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        expand("../mpmat_filter/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{c}2{d}.mpmat",
               sample=SAMPLES, rep=REP,c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        # mpmat_filter_merge
        expand("../mpmat_filter_merge/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.mpmat",
               sample=SAMPLES, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum),
        # sort_mpmat
        expand("../mpmat_filter_merge/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat",
               sample=SAMPLES, rep=REP,a=MPMAT_MERGE[0],b=MPMAT_MERGE[1],c=MPMAT_MERGE[2],d=MPMAT_MERGE[3],
               SMN=SiteMutNum,SCN=SiteCoverNum,SMR=SiteMutRatio,RPN=RegionPassNum,RTN=RegionToleranceNum)
rule mpmat_filter:
    input:
        "../mpmat/293T-bat_{sample}_{rep}_hg38.select.{ac}2{bd}.mpmat"    
    output:
        "../mpmat_filter/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{ac}2{bd}.mpmat"
    params:
        ac = '{ac}',
        bd = '{bd}',
        SiteMutNum = "{SMN}",
        SiteCoverNum = "{SCN}",
        SiteMutRatio = "{SMR}",
        RegionPassNum = "{RPN}",
        RegionToleranceNum = "{RTN}"
    shell:
        """
        {PYTHON2} ./program/mpmat-select-V02.py \
        -i {input} -o {output} \
        -f {params.ac} -t {params.bd} \
        --SiteMutNum {params.SiteMutNum} \
        --SiteCoverNum {params.SiteCoverNum} \
        --SiteMutRatio {params.SiteMutRatio} \
        --RegionPassNum {params.RegionPassNum} \
        --RegionToleranceNum {params.RegionToleranceNum} \
        --InHeader True \
        --OutHeader False
        """

rule mpmat_filter_merge:
    input:
        "../mpmat_filter/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}.mpmat",
        "../mpmat_filter/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{c}2{d}.mpmat"
    output:
        "../mpmat_filter_merge/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.mpmat"
    shell:
        """
        cat {input[0]} {input[1]} > {output}
        """
rule sort_mpmat:
    input:
        "../mpmat_filter_merge/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.mpmat"

    output:
        "../mpmat_filter_merge/293T-bat_{sample}_{rep}_hg38.filtered_SMN-{SMN}_SCN-{SCN}_SMR-{SMR}_RPN-{RPN}_RTN-{RTN}.{a}2{b}_merge_{c}2{d}.sort.mpmat"

    shell:
        """
        {BEDTOOLS} sort -i {input} -g {GENOME_FAI} > {output}
        """