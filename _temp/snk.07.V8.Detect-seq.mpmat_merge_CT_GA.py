# _*_ coding: UTF-8 _*_

########################################################################
# Hao Wu
# 2022-09-02
# Li coop A3A DETEC-Seq EMX1
######################################################################## 
# run on abyss
# /lustre2/chengqiyi_pkuhpc/wuhao/02.Projects/01.LDL_coop_A3A/02.20220831_A3A_Detect-seq
# change settings ==========================================>>>>>>>>>>>>>>>>>


### 1.Unchanged Part ----------------------------->>>>>>>>>>>
BWA_HG38_IDX = "/lustre2/chengqiyi_pkuhpc/wuhao/00.Index_and_DB/align_index/hg38/bwa/hg38_only_chromosome.fa"
BWA_HG38_FAI = "/lustre2/chengqiyi_pkuhpc/wuhao/00.Index_and_DB/align_index/hg38/bwa/hg38_only_chromosome.fa.fai"
HG38_JSON = "/lustre2/chengqiyi_pkuhpc/wuhao/00.Index_and_DB/align_index/hg38/hg38.json"
SNP_VCF = "/lustre2/chengqiyi_pkuhpc/wuhao/00.Index_and_DB/mutation_files/293T/293T-Mock-Input-covaris_bwa_hg38_sort_rmdup.recall.merge.Genotype.filter.rmdup_signal.vcf"

PICARD = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/pkgs/picard-2.23.9-0/share/picard-2.23.9-0/picard.jar"
BISMARK = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/bismark"
CUTADAPT = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/cutadapt"
BWA = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/bwa"
SAMTOOLS = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/envs/samtools/bin/samtools"
SAMCLIP = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/samclip"
PYTHON = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/envs/whpy27/bin/python"
PYTHON3 = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/python"
BEDTOOLS = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/bedtools"
JAVA = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/envs/Java8/bin/java"
BOWTIE2 = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/bowtie2"
STAR = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/STAR" # index were builded with v2.7.1a
FASTQC = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/fastqc"
MULTIQC = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/multiqc"
TRIM_GALORE = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/trim_galore"
RSCRIPT = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/envs/R4.1.3/bin/Rscript"
SAMBAMBA = "/lustre2/chengqiyi_pkuhpc/wuhao/01.Biosoft/miniconda3/bin/sambamba"


### 2.Set params like samples, used in pipeline, changed with projects ------------>>>>>>>>>>>

ABY = ["/lustre2/chengqiyi_pkuhpc/wuhao/02.Projects/01.LDL_coop_A3A/02.20220831_A3A_Detect-seq"]

SAMPLES = [
    "A3A-Input-rep1",
    "A3A-PD-rep1",
    "VA-rep1",
    "YE1-rep1",
    "VEGFA-Vector-rep1",
    "Vector-merge"
]

CTRL_SAMPLE = ["Vector-merge"]

MUT_COMBINE = ["CT,GA"]
TYPE_ONE = str(MUT_COMBINE[0]).split(",")[0]
TYPE_TWO = str(MUT_COMBINE[0]).split(",")[1]

# COVER_NUM =["4","8","10","14","18","22"]
# MUT_NUM = ["8","10","12","14","16","18","20","22","24","26","28","30"]
# COVER_NUM =["4","10","18","25"]
# MUT_NUM = ["8","16","24","30"]
# PASS_NUM = ["1"]
# MUT_RATIO = ["0.01"]
COVER_NUM =["15"]
MUT_NUM = ["8"]
PASS_NUM = ["1"]
MUT_RATIO = ["0.01"]


### !!!!!! SMART = Sort, MQPQ20, Aligned, Rmdup, Rm_FiltClip => Final BAM
### !!!!!! FORMS = Filter, Merge, Sort, RmOverlap, RmOneSite => One big mpmat
# 3.Make rule input and ouput files ------------------------>>>>>>>>>>>>
rule all:
    input:
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.mpmat",
               sample=SAMPLES,type_one=TYPE_ONE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.mpmat",
               sample=SAMPLES,type_two=TYPE_TWO,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        
        ### [wh10_sort_merged_mpmat]: sort above aim analyzing merged (CT and GA) mpmat file for each sample, paras => "sample","type_one","type_two"
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}_{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.Merge.mpmat",
               sample=SAMPLES,type_one=TYPE_ONE,type_two=TYPE_TWO,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}_{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.Merge.sort.mpmat",
               sample=SAMPLES,type_one=TYPE_ONE,type_two=TYPE_TWO,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}_{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.FORMS.mpmat",
               sample=SAMPLES,type_one=TYPE_ONE,type_two=TYPE_TWO,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM)

# 4.run pipeline rules -------------------------------------->>>>>>>>>>   
rule wh9_merge_TYPE_ONEandTYPE_TWO_mpmat:
    input:
        type_one = "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.mpmat",
        type_two = "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.mpmat"
    output:
        "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}_{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.Merge.mpmat"
    shell:
        """
        cat {input.type_one} {input.type_two} > {output}
        """


rule wh10_sort_merged_mpmat:
    input:
        "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}_{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.Merge.mpmat"
    output:
        sort_mpmat = "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}_{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.Merge.sort.mpmat",
        rm_chr_mpmat = "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{type_one}_{type_two}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.FORMS.mpmat"
    shell:
        """
        {BEDTOOLS} sort -i {input} -g {BWA_HG38_FAI} | uniq > {output.sort_mpmat} && \
        cat {output.sort_mpmat} | grep -v chrY | grep -v chrM > {output.rm_chr_mpmat}
        """
