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

MUT_COMBINE = ["CT,GA"]
ANALYSIS_TYPE = MUT_COMBINE[0].split(",")
analysis_type = str(MUT_COMBINE[0]).split(",")[0]
TYPE_TWO = str(MUT_COMBINE[0]).split(",")[1]

### 2.Set params like samples, used in pipeline, changed with projects ------------>>>>>>>>>>>

ABY = ["/lustre2/chengqiyi_pkuhpc/wuhao/02.Projects/01.LDL_coop_A3A/02.20220831_A3A_Detect-seq"]

### 样本命中不要出现下划线“_"
SAMPLES = [
    "A3A-Input-rep1",
    "A3A-PD-rep1",
    "VA-rep1",
    "YE1-rep1",
    "VEGFA-Vector-rep1",
    "Vector-merge"
]

CTRL_SAMPLE = ["Vector-merge"]

COVER_NUM =["15"]
MUT_NUM = ["8"]
PASS_NUM = ["1"]
MUT_RATIO = ["0.01"]

# 3.Make rule input and ouput files ------------------------>>>>>>>>>>>>
rule all:
    input:
        ### [wh7_mpmat_filt_CT]: filter 16 mutation types mpmat file based on mut num. mut ratio, cover num, paras => "sample","ref_base","mut_base"
        expand("../3.Process/6.mpmat/{sample}/{sample}_{analysis_type}_bwa_hg38_sort_rmdup_MAPQ20.mpmat",sample=SAMPLES,analysis_type=ANALYSIS_TYPE),
        expand("../3.Process/6.mpmat/{sample}/Filt_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/{sample}_{analysis_type}_bwa_hg38_sort_rmdup_MAPQ20.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.filtered.mpmat",
               sample=SAMPLES,analysis_type=ANALYSIS_TYPE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        expand("../3.Process/6.mpmat/{sample}/Filt_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/{sample}_{analysis_type}_bwa_hg38_sort_rmdup_MAPQ20.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.filtered.mpmat.log",
               sample=SAMPLES,analysis_type=ANALYSIS_TYPE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        
        ### cat, sort big one merged mpmat and rm overlap and onesite ###
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.mpmat",
               sample=SAMPLES,analysis_type=ANALYSIS_TYPE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.mpmat",
               sample=SAMPLES,analysis_type=ANALYSIS_TYPE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.mpmat",
               sample=SAMPLES,analysis_type=ANALYSIS_TYPE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM),
        expand("../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.mpmat",
               sample=SAMPLES,analysis_type=ANALYSIS_TYPE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM)
        
        

# 4.run pipeline rules -------------------------------------->>>>>>>>>>>
rule wh7_mpmat_filt_CT_GA:
    input:
        "../3.Process/6.mpmat/{sample}/{sample}_{analysis_type}_bwa_hg38_sort_rmdup_MAPQ20.mpmat"
    output:
        "../3.Process/6.mpmat/{sample}/Filt_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/{sample}_{analysis_type}_bwa_hg38_sort_rmdup_MAPQ20.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.filtered.mpmat"
    log:
        "../3.Process/6.mpmat/{sample}/Filt_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/{sample}_{analysis_type}_bwa_hg38_sort_rmdup_MAPQ20.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.filtered.mpmat.log"
    params:
        ref = lambda wildcards, input: input[0].split("/")[-1].split("_")[-6][0],
        mut = lambda wildcards, input: input[0].split("/")[-1].split("_")[-6][1],
        filt = "../2.Codes/mpmat-select.py",
        mutnum = "{mut_num}",
        covernum = "{cover_num}",
        mutratio = "{mut_ratio}",
        passnum = "{pass_num}"
    shell:
        """
        {PYTHON3} {params.filt} \
        -i {input} -o {output} -f {params.ref} -t {params.mut} \
        -m {params.mutnum} -c {params.covernum} -r {params.mutratio} \
        --RegionPassNum {params.passnum} --RegionToleranceNum 10 \
        --RegionMutNum 2 --InHeader True --OutHeader False > {log} 2>&1
        """

        
rule Big_Merge:
    input:
        expand("../3.Process/6.mpmat/{sample}/Filt_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/{sample}_{analysis_type}_bwa_hg38_sort_rmdup_MAPQ20.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.filtered.mpmat",
               sample=SAMPLES,analysis_type=ANALYSIS_TYPE,mut_num=MUT_NUM,cover_num=COVER_NUM,mut_ratio=MUT_RATIO,pass_num=PASS_NUM)
    output:
        "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.mpmat"
    params:
        files = lambda wildcards, output: " ".join("/".join(output[0].split("/")[0:2] +  ["6.mpmat"] + [i] + ["Filt_mpmat"] + output[0].split("/")[3:6] + [i+"_"+output[0].split("/")[-1].split(".")[2]+"_bwa_hg38_sort_rmdup_MAPQ20."+output[0].split("/")[-1].split(output[0].split("/")[-1].split(".")[2]+".")[1].split(".OneBig")[0]+".filtered.mpmat"]) for i in SAMPLES)
    shell:
        """
        cat {params.files} > {output}
        """
        
        

rule Sort_big_merge:
    input:
        "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.mpmat"
    output:
        "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.mpmat"
    shell:
        """
        {BEDTOOLS} sort -i {input} -g {BWA_HG38_FAI} | uniq > {output}
        """
        
        
rule rm_overlap:
    input:
         "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.mpmat"
    output:
         "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.mpmat"
    shell:
        """
        {PYTHON3} ../2.Codes/remove-overlap-mpmat.py -i {input} -o {output} 
        """
        
rule rm_one_site:
    input:
        "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.mpmat"
    output:
        "../3.Process/7.Merge_One_big_mpmat/mut_num_{mut_num}/cover_num_{cover_num}/mut_ratio_{mut_ratio}/All.SMART.{analysis_type}.M{mut_num}C{cover_num}R{mut_ratio}P{pass_num}.OneBig.sort.Rmoverlap.RmOneSite.mpmat"
    shell:
        """
        awk '$5 > 1 {{print $0}}' {input} > {output}
        """ 
