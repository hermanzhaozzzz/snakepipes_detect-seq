# ——————————————————>>>>>>>>>>
# Project Name: Detect-seq
# Author: Hua-nan ZHAO
# E-mail: hermanzhaozzzz@gmail.com
# Update log:
#     2022-10-05: start project

# test run:
# ```shell
# snakemake -pr -j 4 -s step.02.Snakefile.smk.py -n
# ```
# ——————————————————>>>>>>>>>>
import os
import json


# ------------------------------------------------------------------->>>>>>>>>>
# FUNCTIONS
# ------------------------------------------------------------------->>>>>>>>>>
def print_head(SAMPLES, MODE):
    print('----------\nSAMPLES:')
    [print('\t' + i) for i in SAMPLES]
    print('----------\nMODE:')
    print('\t' + MODE)
    print('----------\n\n')


def check_cmd(x):
    return any(
        os.access(os.path.join(path,x),os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )


def check_read(x):
    if x == "PE":
        read = ['R1', 'R2']
    elif x == "SE":
        read = ['SE']
    else:
        raise ValueError()
    return read


# ------------------------------------------------------------------->>>>>>>>>>
# SAMPLE INFO
# ------------------------------------------------------------------->>>>>>>>>>
with open('./samples.json') as f:
    dt = json.loads(f.read())

SAMPLES = dt['samples']
MODE = dt['seq_mode']
THREAD = dt['thread']
READ = check_read(MODE)

print_head(SAMPLES,MODE)
print(READ)
# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
PLATFORM = dt['platform']
# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
GENOME_HISAT3N_INDEX = dt['genome_hisat3n_index']
BASE_CHANGE = dt['base_change']
SNP_LIST = dt['snp_list']
CTRL = dt['control']
CTRL_NAME = CTRL["name"]
CTRL_BAM = CTRL["bam"]
QMT = dt['query_mutation_type']
# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# check if cmd exists
assert check_cmd("cutadapt")
assert check_cmd("hisat-3n")
assert check_cmd("samtools")
assert check_cmd("bwa")
assert check_cmd("samclip")
assert check_cmd("sambamba")
assert check_cmd("bedtools")
# manually set cmd path
CUTADAPT = "cutadapt"
HISAT3N = "hisat-3n"
SAMTOOLS = "samtools"
BWA = "bwa"
SAMCLIP = "samclip"
SAMBAMBA = "sambamba"  # must > 0.7.1 on OSX
# SAMBAMBA = '/usr/local/Caskroom/miniconda/base/envs/sambamba/bin/sambamba'  # must > 0.7.1 on OSX
BEDTOOLS = "bedtools"





# ------------------------------------------------------------------->>>>>>>>>>
# test filter
# ------------------------------------------------------------------->>>>>>>>>>
# C
C = [
    '6', 
    '8', 
    '10'
]  # SITECOVERNUM
#                         Site-cutoff, total cover reads number, default=5 
# M
M = [
    '4', 
    '5', 
    '6', 
    '7'
]  # SITEMUTNUM
#                         Site-cutoff, mutation reads number, default=1
# R
R = [
    '1', 
    '2', 
    '3'
]  # REGIONPASSNUM
#                         Region-cutoff region should contain no-less than Pass site number,
#                         default=3
FILTER = []
for c in C:
    for m in M:
        for r in R:
            FILTER.append(f'C{c}M{m}R{r}')
print(f'FILTER = \n{FILTER}')

dt_samples = {}
for flt in FILTER:
    dt_samples[flt] = SAMPLES

# exit()
# ------------------------------------------------------------------->>>>>>>>>>
# rule all
# ------------------------------------------------------------------->>>>>>>>>>
rule all:
    input:
        expand("../mpmat/{sample}_merge_C{c}M{m}R{r}.select.sort_rmchrYM.mpmat.gz", c=C, m=M, r=R, sample=SAMPLES),
        expand("../mpmat/mpmat_all_{flt}_before_rmdup.mpmat.gz", flt=FILTER),
        expand("../mpmat/mpmat_all_{flt}_after_rmdup.mpmat.gz", flt=FILTER)
    
rule mpmat_filter_1:
    input:
        "../pmat/{sample}_CT.mpmat.gz"
    output:
        "../mpmat/{sample}_CT_C{c}M{m}R{r}.mpmat.gz"
    params:
        C=lambda wildcards: wildcards['c'],
        M=lambda wildcards: wildcards['m'],
        R=lambda wildcards: wildcards['r']
    log:
        "../mpmat/{sample}_CT_C{c}M{m}R{r}.mpmat.log"
    shell:
        """
        python program/detect_seq/mpmat-select.py \
            -i {input} \
            -o {output} \
            -f C \
            -t T \
            -m {params.M} \
            -c {params.C} \
            -r 0.01 \
            --RegionPassNum {params.R} \
            --RegionToleranceNum 10 \
            --RegionMutNum 2 \
            --InHeader True \
            --OutHeader False
        """
rule mpmat_filter_2:
    input:
        "../pmat/{sample}_GA.mpmat.gz"
    output:
        "../mpmat/{sample}_GA_C{c}M{m}R{r}.mpmat.gz"
    params:
        C=lambda wildcards: wildcards['c'],
        M=lambda wildcards: wildcards['m'],
        R=lambda wildcards: wildcards['r']
    log:
        "../mpmat/{sample}_GA_C{c}M{m}R{r}.mpmat.log"
    shell:
        """
        python program/detect_seq/mpmat-select.py \
            -i {input} \
            -o {output} \
            -f G \
            -t A \
            -m {params.M} \
            -c {params.C} \
            -r 0.01 \
            --RegionPassNum {params.R} \
            --RegionToleranceNum 10 \
            --RegionMutNum 2 \
            --InHeader True \
            --OutHeader False
        """   
# python program/detect_seq/mpmat-select.py -h
rule merge_mpmat_per_sample:
    input:
        "../mpmat/{sample}_CT_C{c}M{m}R{r}.mpmat.gz",
        "../mpmat/{sample}_GA_C{c}M{m}R{r}.mpmat.gz"
    output:
        "../mpmat/{sample}_merge_C{c}M{m}R{r}.select.mpmat.gz",
        "../mpmat/{sample}_merge_C{c}M{m}R{r}.select.sort.mpmat.gz",
        "../mpmat/{sample}_merge_C{c}M{m}R{r}.select.sort_rmchrYM.mpmat.gz",
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa.fai'
    shell:
        """
        cat {input[0]} {input[1]}  > {output[0]}

        {BEDTOOLS} sort -i {output[0]} -g {params.ref} | uniq | gzip > {output[1]}

        gunzip -d -c {output[1]} | grep -v chrY | grep -v chrM | gzip > {output[2]}
        """
rule merge_mpmat_all_samples:
    input:
        expand("../mpmat/{sample}_merge_{flt}.select.sort_rmchrYM.mpmat.gz", flt=FILTER, sample=SAMPLES)
    output:
        "../mpmat/mpmat_all_{flt}_before_rmdup.mpmat.gz"
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa.fai',
        inputs=lambda wildcards, input: [i for i in input if wildcards['flt'] in i]
    log:
        "../mpmat/mpmat_all_{flt}_before_rmdup.mpmat.log"
    shell:
        """
        echo "These files were merged:\n" > {log}
        echo {params.inputs} > {log}
        cat {params.inputs} > {output}_temp.gz
        {BEDTOOLS} sort -i {output}_temp.gz -g {params.ref} | uniq | gzip > {output}
        rm {output}_temp.gz
        """
rule mpmat_all_rmdup:
    input:
        "../mpmat/mpmat_all_{flt}_before_rmdup.mpmat.gz"
    output:
        "../mpmat/mpmat_all_{flt}_after_rmdup.mpmat.gz"
    shell:
        "python program/detect_seq/remove-overlap-mpmat.py -i {input} -o {output}"
