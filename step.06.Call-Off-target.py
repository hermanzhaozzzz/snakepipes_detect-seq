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
from glob import glob


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
MPMAT_SPLIT_NUMBERS = list(map(str, range(50)))
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
SAMBAMBA = "sambamba"  # must >= 0.7.1 on OSX
# SAMBAMBA = '/usr/local/Caskroom/miniconda/base/envs/sambamba/bin/sambamba'  # must >= 0.7.1 on OSX
BEDTOOLS = "bedtools"
# bioat must >= 0.1.4.0
rule all:
    input:
        expand("../poisson_res/{sample}_vs_ctrl_%s.pvalue_table" % CTRL_NAME, sample=SAMPLES),
        "../poisson_res/poisson_res_all.tsv.gz"


# ------------------------------------------------------------------->>>>>>>>>>
# run the Poisson enrichment test
# ------------------------------------------------------------------->>>>>>>>>>

rule find_significant_mpmat:
    input:
        # mpmat="../mpmat/mpmat_all_rmdup.mpmat.gz",
        mpmat="../mpmat/mpmat_all_C6M6R2_after_rmdup_rm_oldoverlap_rmblacklist_add_NatIndex_sort2_rmdup2.mpmat",
        bam="../bam/{sample}_final_rmdup.bam"
    output:
        "../poisson_res/{sample}_vs_ctrl_%s.pvalue_table" % CTRL_NAME
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa'
    log:
        "../poisson_res/{sample}_vs_ctrl_%s.log" % CTRL_NAME
    shell:
        """
        python program/detect_seq/find-significant-mpmat.py \
            -p {THREAD} \
            -i {input.mpmat} \
            -o {output} \
            -c {CTRL_BAM} \
            -t {input.bam} \
            -r {params.ref} \
            --query_mutation_type {QMT}  \
            --mpmat_filter_info_col_index -1 \
            --mpmat_block_info_col_index -1  \
            --region_block_mut_num_cutoff 2  \
            --query_mut_min_cutoff 2  \
            --query_mut_max_cutoff 16  \
            --total_mut_max_cutoff 16  \
            --other_mut_max_cutoff 6   \
            --seq_reads_length 150  \
            --lambda_method ctrl_max \
            --poisson_method mutation > {log} 2>&1
        """
# ------------------------------------------------------------------->>>>>>>>>>
# poisson_method
#  Can be set as 'mutation' OR 'all', default=mutation. 'mutation'
#                         means only use mutation alignments to run Poisson test,'all' means
#                         use all alignments to run Poisson, which similar to MACS2
# ------------------------------------------------------------------->>>>>>>>>>
rule merge_poisson_res_table:
    input:
        expand("../poisson_res/{sample}_vs_ctrl_%s.pvalue_table" % CTRL_NAME, sample=SAMPLES)
    output:
        "../poisson_res/poisson_res_all.tsv.gz"
    params:
        inputs=lambda wildcards, input: ",".join(input),
        tags=lambda wildcards, input: ",".join([i.split('/')[-1].split('_vs_ctrl')[0] for i in input])
    shell:
        "bioat table merge {params.inputs} {params.tags} {output} --input_header True --output_header True"
