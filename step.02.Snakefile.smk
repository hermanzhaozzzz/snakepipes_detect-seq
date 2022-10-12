# ——————————————————>>>>>>>>>>
# Project Name: Detect-seq
# Author: Hua-nan ZHAO
# E-mail: hermanzhaozzzz@gmail.com
# Update log:
#     2022-10-05: start project
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
        os.access(os.path.join(path, x), os.X_OK) 
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

print_head(SAMPLES, MODE)
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
# manually set cmd path
CUTADAPT = "cutadapt"
HISAT3N = "hisat-3n"
SAMTOOLS = "samtools"
BWA = "bwa"
SAMCLIP = "samclip"
SAMBAMBA = "sambamba"
# ------------------------------------------------------------------->>>>>>>>>>
# rule all
# ------------------------------------------------------------------->>>>>>>>>>
rule all:
    input:
        # raw input
        expand("../fastq/{sample}_%s.fastq.gz" % READ[0], sample=SAMPLES),
        # cutadapt
        expand("../fastq_temp/{sample}_%s.fastq.gz" % READ[0], sample=SAMPLES),
        # mapping by hisat3n
        expand("../bam/{sample}_hisat3n.bam", sample=SAMPLES),
        # remapped by bwa mem
        expand("../bam/{sample}_remapped.bam", sample=SAMPLES),
        expand("../bam/{sample}_final_rmdup.bam", sample=SAMPLES),
        expand("../pmat/{sample}.pmat.gz", sample=SAMPLES),
        expand("../mpmat/{sample}_merge.select.sort.mpmat.gz", sample=SAMPLES),
        expand("../mpmat/{sample}_merge.select.sort_rmchrYM.mpmat.gz", sample=SAMPLES),
        expand("../poisson_res/{sample}_vs_ctrl_%s.select.pvalue_table" % CTRL_NAME, sample=SAMPLES),
        
# ------------------------------------------------------------------->>>>>>>>>>
# rule cutadapt
# ------------------------------------------------------------------->>>>>>>>>>
rule cutadapt: 
    input: 
        "../fastq/{sample}_%s.fastq.gz" % READ[0]
    output: 
        R1="../fastq_temp/{sample}_R1.fastq.gz",
        R2="../fastq_temp/{sample}_R2.fastq.gz"
    params:
        SE="../fastq/{sample}_SE.fastq.gz",
        R1="../fastq/{sample}_R1.fastq.gz",
        R2="../fastq/{sample}_R2.fastq.gz"
    log:
        "../fastq_temp/{sample}.log"
    shell:
        """
        mkdir -p ../fastq_temp
        INPUT={input}
        if [[ $INPUT =~ .*SE.fastq.gz$ ]]; then 
            echo "[FATAL] find SE reads, raw reads should be PE reads in Detect-seq protocol" > {log}
        elif [[ $INPUT =~ .*R1.fastq.gz$ ]]; then
            echo "[DEBUG] find PE reads, go on CUTADAPT" > {log}
            {CUTADAPT} -j {THREAD} --times 1  -e 0.1  -O 3  --quality-cutoff 25 -m 55 \
                -a AGATCGGAAGAGCACACGT  \
                -A  AGATCGGAAGAGCGTCGTG \
                -o {output.R1} \
                -p {output.R2} \
                {params.R1} \
                {params.R2} >> {log} 2>&1
        else
            echo "[FATAL] fastq should be *.SE.fastq.gz or *.R[1,2].fastq.gz" > {log}
        fi
        echo "[DEBUG] CUTADAPT done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# rule hisat3n
# ------------------------------------------------------------------->>>>>>>>>>
# Map the Detect-seq sequencing reads to the reference genome by HISAT-3N,
# which uses a base conversion alignment strategy.
# ------------------------------------------------------------------->>>>>>>>>>
rule hisat3n:
    input:
        "../fastq_temp/{sample}_%s.fastq.gz" % READ[0]
    output:
        sam=temp("../bam/{sample}_hisat3n.sam"),
        unmapped_fq1=temp("../fastq_temp/{sample}_unmapped.fastq.1.gz"),
        unmapped_fq2=temp("../fastq_temp/{sample}_unmapped.fastq.2.gz"),
    params:
        SE="../fastq_temp/{sample}_SE.fastq.gz",
        R1="../fastq_temp/{sample}_R1.fastq.gz",
        R2="../fastq_temp/{sample}_R2.fastq.gz",
        unmapped_fq="../fastq_temp/{sample}_unmapped.fastq.gz",
        rg_id='{sample}',
        sample='{sample}'
    log:
        hisat3n="../bam/{sample}_hisat3n.hisat3n.log",
        smk="../bam/{sample}.smk.log"
    shell:
        """
        INPUT={input}
        if [[ $INPUT =~ .*SE.fastq.gz$ ]]; then
            echo "[FATAL] find SE reads, raw reads should be PE reads in Detect-seq protocol" > {log.smk}
        elif [[ $INPUT =~ .*R1.fastq.gz$ ]]; then
            echo "[DEBUG] find PE reads, go on hisat-3n mapping" > {log.smk}
            {HISAT3N} \
                --index {GENOME_HISAT3N_INDEX} \
                -1 {params.R1} \
                -2 {params.R2} \
                -S {output.sam} \
                -p {THREAD} \
                --sensitive \
                --base-change {BASE_CHANGE} \
                --unique-only \
                --repeat \
                --no-repeat-index \
                --repeat-limit 1000 \
                --no-spliced-alignment \
                -X 700 \
                --un-conc-gz {params.unmapped_fq} \
                --summary-file {log.hisat3n} \
                --rg-id "{params.rg_id}" \
                --rg "PL:{PLATFORM}" \
                --rg "ID:{params.rg_id}" \
                --rg "SM:{params.rg_id}"
        else
            echo "[FATAL] fastq should be *.SE.fastq.gz or *.R[1,2].fastq.gz" > {log.smk}
        fi
        echo "[DEBUG] mapping done" >> {log.smk}
        """
rule sam2bam:
    input: "../bam/{sample}_hisat3n.sam"
    output: temp("../bam/{sample}_hisat3n.bam")
    shell: "{SAMTOOLS} view -@ {THREAD} -Shb {input} > {output}"
# ------------------------------------------------------------------->>>>>>>>>>
# rule filter_low_mapq_reads
# ------------------------------------------------------------------->>>>>>>>>>
# select low mapping quality BAM
# ------------------------------------------------------------------->>>>>>>>>>
rule filter_low_mapq_reads:
    input:
        "../bam/{sample}_hisat3n.bam"
    output:
        temp("../bam/{sample}_hisat3n.LowerMAPQ20.bam")
    params:
        awk="""'$1~"@" || $5 <= 20  {print $0}'"""
    shell:
        """
         {SAMTOOLS} view -h -@ {THREAD} {input} \
             | awk {params.awk} \
             |  {SAMTOOLS} view -@ {THREAD} -hb > {output}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# BAM sort by reads name
# ------------------------------------------------------------------->>>>>>>>>>
rule low_mapq_reads_sortn:
    input: "../bam/{sample}_hisat3n.LowerMAPQ20.bam"
    output: temp("../bam/{sample}_hisat3n.LowerMAPQ20.sortn.bam")
    shell:
        """
        mkdir -p ../temp_file
        {SAMTOOLS} sort -O BAM -o {output} -T ../temp_file -@ {THREAD} -n {input}
            # -m 2G \
        """
# ------------------------------------------------------------------->>>>>>>>>>
# fetch low mapping quality reads from BAM file
# ------------------------------------------------------------------->>>>>>>>>>
rule fetch_low_mapq_reads_from_bam:
    input: "../bam/{sample}_hisat3n.LowerMAPQ20.sortn.bam"
    output:
        R1=temp("../fastq_temp/{sample}_low_mapq_R1.fastq.gz"),
        R2=temp("../fastq_temp/{sample}_low_mapq_R2.fastq.gz")
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa'
    shell:
        """
        {SAMTOOLS} fastq -@ {THREAD} -0 /dev/null -s /dev/null -n -F 0x900 \
            -1 {output.R1} -2 {output.R2} --reference {params.ref} {input}
        """

# ------------------------------------------------------------------->>>>>>>>>>
# merge unmapped reads and low mapping quality reads
# ------------------------------------------------------------------->>>>>>>>>>
rule merge_unmapped_and_lowmapq_reads:
    input:
        unmapped_fq1="../fastq_temp/{sample}_unmapped.fastq.1.gz",
        unmapped_fq2="../fastq_temp/{sample}_unmapped.fastq.2.gz",
        lowmapq_fq1="../fastq_temp/{sample}_low_mapq_R1.fastq.gz",
        lowmapq_fq2="../fastq_temp/{sample}_low_mapq_R2.fastq.gz"
    output:
        R1=temp("../fastq_temp/{sample}_unmapped_plus_lowmapq_R1.fastq.gz"),
        R2=temp("../fastq_temp/{sample}_unmapped_plus_lowmapq_R2.fastq.gz")
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa'
    shell:
        """
        cat {input.unmapped_fq1} {input.lowmapq_fq1} > {output.R1}
        cat {input.unmapped_fq2} {input.lowmapq_fq2} > {output.R2}
        """

# ------------------------------------------------------------------->>>>>>>>>>
# re-alignment with BWA MEM
# ------------------------------------------------------------------->>>>>>>>>>
rule remap_unmapped_and_lowmapq_reads:
    input:
        R1="../fastq_temp/{sample}_unmapped_plus_lowmapq_R1.fastq.gz",
        R2="../fastq_temp/{sample}_unmapped_plus_lowmapq_R2.fastq.gz"
    output:
        temp("../bam/{sample}_remapped.sam")
    log:
        "../bam/{sample}_remapped.log"
    params:
        # genome_ucsc_hg38.fa.bwa_index
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa.bwa_index',
        awk=r"""'@RG\tID:'{sample}'\tPL:%s\tSM:{sample}'""" % PLATFORM
    shell:
        """
        {BWA} mem {params.ref} {input.R1} {input.R2} -t {THREAD} -M \
            -R {params.awk} > {output} 2>{log}
        """
rule sam2bam_2:
    input: "../bam/{sample}_remapped.sam"
    output: temp("../bam/{sample}_remapped.bam")
    shell: "{SAMTOOLS} view -@ {THREAD} -Shb {input} > {output}"

# ------------------------------------------------------------------->>>>>>>>>>
# merge HISAT-3n BAM and BWA MEM BAM
# ------------------------------------------------------------------->>>>>>>>>>
# Merge the HISAT-3N alignment BAM file and the BWA MEM alignment BAM file 
# together, and then remove the duplication reads, non-concordant reads, 
# clipped mapping reads and the secondary alignment reads. 
# Only reads with a mapping quality larger than 20 are retained for 
# subsequent analysis.
# ------------------------------------------------------------------->>>>>>>>>>
rule merge_bams:
    input: 
        raw_bam = "../bam/{sample}_hisat3n.bam",
        remapped_bam = "../bam/{sample}_remapped.bam"
    output:
        temp("../bam/{sample}_merged.bam")
    shell:
        "{SAMTOOLS} cat -o {output} {input.raw_bam} {input.remapped_bam}"

# ------------------------------------------------------------------->>>>>>>>>>
# sort BAM by genome coordinate
# ------------------------------------------------------------------->>>>>>>>>>
rule sort_bam:
    input: "../bam/{sample}_merged.bam"
    output: temp("../bam/{sample}_merged_sorted.bam")
    shell: 
        """
        mkdir -p ../temp_file
        {SAMTOOLS} sort -O BAM -o {output} -T ../temp_file -@ {THREAD} {input}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# filter clip, non-concordant reads, low MAPQ reads and secondary alignment
# ------------------------------------------------------------------->>>>>>>>>>
rule filter_bam:
    input: "../bam/{sample}_merged_sorted.bam"
    output: temp("../bam/{sample}_merged_sorted_filter.bam")
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa.fai',
        awk="""'function abs(v) {return v < 0 ? -v : v} $1~"@" || ($7 == "=" && abs($9) <= 2500 ) {print $0}'"""
    shell: 
        """
        {SAMTOOLS} view -@ {THREAD} -h {input} -q 20 -f 3 -F 256 \
            | {SAMCLIP} --ref {params.ref} --max 3 --progress 0 \
            | awk {params.awk} | {SAMTOOLS} view -@ {THREAD} -hb > {output}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# remove duplications and build BAM index
# ------------------------------------------------------------------->>>>>>>>>>
rule sambamba_rmdup_and_build_index:
    input: "../bam/{sample}_merged_sorted_filter.bam"
    output: 
        "../bam/{sample}_final_rmdup.bam",
        "../bam/{sample}_final_rmdup.bam.bai"
    params:
        "../bam/{sample}_final_rmdup.bam"
    shell:
        """
        {SAMBAMBA} markdup \
                --remove-duplicates \
                --nthreads={THREAD} \
                --show-progress \
                --sort-buffer-size 8192 \
                {input} {output[0]}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# Tracing tandem C-to-T signals and perform statistical test
# ------------------------------------------------------------------->>>>>>>>>>
# Catch contiguous tandem mutation information from BAM files (i.e., tandem 
# C-to-T mutation signals for the Watson strand and tandem G-to-A mutation 
# signals for the Crick strand). 
# Then perform a comparison for those tandem mutation signals between control 
# sample and Detect-seq treatment sample by the Poisson statistical test. 
# After this step, one can obtain tables containing information of tandem 
# mutation signals (a mpmat file shown in the code part) and matched Poisson 
# test results (a TSV table shown in the code part).
# ------------------------------------------------------------------->>>>>>>>>>
# convert BAM to pmat format
# ------------------------------------------------------------------->>>>>>>>>>
rule bam2pmat:
    input:
        "../bam/{sample}_final_rmdup.bam"
    output:
        "../pmat/{sample}.pmat.gz"
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa'
    shell:
        """
        python program/detect_seq/bam2pmat.py \
            -i {input} \
            -r {params.ref} \
            -o {output} \
            -p {THREAD} \
            --out_format pmat \
            --bed_like_format True \
            --mut_type ALL \
            --block_size 1000000 \
            --cover_num_cutoff 0 \
            --mut_num_cutoff 0 \
            --mut_ratio_cutoff 0 \
            --keep_temp_file False \
            --out_header False
        """
# ------------------------------------------------------------------->>>>>>>>>>
# tracing tandem mutation signals
# merge tandem mutation signals and sort
# ------------------------------------------------------------------->>>>>>>>>>
# ------------------------------------------------------------------->>>>>>>>>>
# do  
#     in_CT_mpmat=pmat_and_mpmat/${sample}_hg38.MAPQ20.merge_d50_D100.CT.mpmat
#     in_GA_mpmat=pmat_and_mpmat/${sample}_hg38.MAPQ20.merge_d50_D100.GA.mpmat
#     ref_genome_fa_index=reference/hisat3n_hg38_CT/hg38_only_chromosome.fa.fai

#     out_CT_mpmat=pmat_and_mpmat/${sample}_hg38.MAPQ20.merge_d50_D100.CT.select.mpmat
#     out_GA_mpmat=pmat_and_mpmat/${sample}_hg38.MAPQ20.merge_d50_D100.GA.select.mpmat    
#     out_merge_mpmat=pmat_and_mpmat/${sample}_hg38.MAPQ20.merge_d50_D100.merge.select.mpmat
# out_merge_sort_mpmat=pmat_and_mpmat/${sample}_hg38.MAPQ20.merge.select.sort.mpmat
# out_rm_chr_mpmat=pmat_and_mpmat/${sample}_hg38.MAPQ20.merge.select.sort.RmChrYChrM.mpmat
#     # select CT
#     python mpmat-select.py -i {in_CT_mpmat} -o ${out_CT_mpmat} -f C -t T -m 4 -c 6 -r 0.01 --RegionPassNum 1 --RegionToleranceNum 10 --RegionMutNum 2 --InHeader True --OutHeader False

#     # select GA
#     python mpmat-select.py -i {in_GA_mpmat} -o ${out_GA_mpmat} -f G -t A -m 4 -c 6 -r 0.01 --RegionPassNum 1 --RegionToleranceNum 10 --RegionMutNum 2 --InHeader True --OutHeader False
#  # merge CT singal on the Watson strand and the Crick strand
#     cat ${out_CT_mpmat} ${out_GA_mpmat}  > ${out_merge_mpmat}

#     # sort by the genome coordinate
#     bedtools sort -i ${out_merge_mpmat} -g ${ref_genome_fa_index} | uniq > ${out_merge_sort_mpmat}

#     # remove chrY and chrM
#     cat ${out_merge_sort_mpmat} | grep -v chrY | grep -v chrM > ${out_rm_chr_mpmat}
# done
rule pmat_merge_1:
    input:
        "../pmat/{sample}.pmat.gz"
    output:
        raw=temp("../pmat/{sample}_CT.mpmat.gz"),
        filtered=temp("../pmat/{sample}_CT_select.mpmat.gz")
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa',
        f='C',
        t='T'
    log:
        step1="../pmat/{sample}_CT.mpmat.gz.log",
        step2="../pmat/{sample}_CT_select.mpmat.gz.log"
    shell:
        """
        python program/detect_seq/pmat-merge.py \
            --Input {input} \
            --FromBase {params.f} \
            --ToBase {params.t} \
            --reference {params.ref} \
            --MaxSiteDistance 50 \
            --MaxRegionDistance 100 \
            --NoMutNumCutoff 2 \
            --OmitTandemNumCutoff 2 \
            --Output {output.raw} \
            --SNP {SNP_LIST} > {log.step1} 2>&1
        
        python program/detect_seq/mpmat-select.py \
            -i {output.raw} \
            -o {output.filtered} \
            -f {params.f} \
            -t {params.t} \
            -m 4 \
            -c 6 \
            -r 0.01 \
            --RegionPassNum 1 \
            --RegionToleranceNum 10 \
            --RegionMutNum 2 \
            --InHeader True \
            --OutHeader False > {log.step2} 2>&1
        """
rule pmat_merge_2:
    input:
        "../pmat/{sample}.pmat.gz"
    output:
        raw=temp("../pmat/{sample}_GA.mpmat.gz"),
        filtered=temp("../pmat/{sample}_GA_select.mpmat.gz")
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa',
        f='G',
        t='A'
    log:
        step1="../pmat/{sample}_GA.mpmat.gz.log",
        step2="../pmat/{sample}_GA_select.mpmat.gz.log"
    shell:
        """
        python program/detect_seq/pmat-merge.py \
            --Input {input} \
            --FromBase {params.f} \
            --ToBase {params.t} \
            --reference {params.ref} \
            --MaxSiteDistance 50 \
            --MaxRegionDistance 100 \
            --NoMutNumCutoff 2 \
            --OmitTandemNumCutoff 2 \
            --Output {output.raw} \
            --SNP {SNP_LIST} > {log.step1} 2>&1
        
        python program/detect_seq/mpmat-select.py \
            -i {output.raw} \
            -o {output.filtered} \
            -f {params.f} \
            -t {params.t} \
            -m 4 \
            -c 6 \
            -r 0.01 \
            --RegionPassNum 1 \
            --RegionToleranceNum 10 \
            --RegionMutNum 2 \
            --InHeader True \
            --OutHeader False > {log.step2} 2>&1
        """
rule merge_mpmat:
    input:
        "../pmat/{sample}_CT_select.mpmat.gz",
        "../pmat/{sample}_GA_select.mpmat.gz"
    output:
        temp("../mpmat/{sample}_merge.select.mpmat.gz"),
        "../mpmat/{sample}_merge.select.sort.mpmat.gz",
        "../mpmat/{sample}_merge.select.sort_rmchrYM.mpmat.gz",
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa.fai'
    shell:
        """
        cat {input[0]} {input[1]}  > {output[0]}
        
        bedtools sort -i {output[0]} -g {params.ref} | uniq | gzip > {output[1]}
        
        zcat {output[1]} | grep -v chrY | grep -v chrM | gzip > {output[2]}
        """

# ------------------------------------------------------------------->>>>>>>>>>
# run the Poisson enrichment test
# ------------------------------------------------------------------->>>>>>>>>>
# ref_genome_fa=reference/hisat3n_hg38_CT/hg38_only_chromosome.fa

# # for BE4max sample
# python find-significant-mpmat.py -p 25 \
# -i pmat_and_mpmat/293T-BE4max-VEGFA-All-PD_hg38.MAPQ20.merge.select.sort.RmChrYChrM.mpmat \
# -o poisson_res/293T-BE4max-VEGFA-All-PD_vs_ctrl_hg38.select.pvalue_table \
# -c bam.hisat3n/293T-BE4max-mCherry-PD_hg38_merge_sort_rmdup.MAPQ20.bam \
# -t bam.hisat3n/293T-BE4max-VEGFA-All-PD_hg38_merge_sort_rmdup.MAPQ20.bam \
# -r ${ref_genome_fa} \
# --query_mutation_type CT,GA  \
# --mpmat_filter_info_col_index -1 \
# --mpmat_block_info_col_index -1  \
# --region_block_mut_num_cutoff 2  \
# --query_mut_min_cutoff 2  \
# --query_mut_max_cutoff 16  \
# --total_mut_max_cutoff 16  \
# --other_mut_max_cutoff 6   \
# --seq_reads_length 150  \
# --lambda_method ctrl_max \
# --poisson_method mutation \
# 2> poisson_res/293T-BE4max-VEGFA-All-PD_vs_ctrl_hg38_possion_test.log &

# # for DdCBE sample
# python find-significant-mpmat.py -p 25 \
# -i pmat_and_mpmat/293T-DdCBE-ND6-All-PD_hg38.MAPQ20.merge.select.sort.RmChrYChrM.mpmat \
# -o poisson_res/293T-DdCBE-ND6-All-PD_vs_ctrl_hg38.select.pvalue_table \
# -c bam.hisat3n/293T-DdCBE-GFP-PD_hg38_merge_sort_rmdup.MAPQ20.bam \
# -t bam.hisat3n/293T-DdCBE-ND6-All-PD_hg38_merge_sort_rmdup.MAPQ20.bam \
# -r ${ref_genome_fa} \
# --query_mutation_type CT,GA  \
# --mpmat_filter_info_col_index -1 \
# --mpmat_block_info_col_index -1  \
# --region_block_mut_num_cutoff 2  \
# --query_mut_min_cutoff 2  \
# --query_mut_max_cutoff 16  \
# --total_mut_max_cutoff 16  \
# --other_mut_max_cutoff 6   \
# --seq_reads_length 150  \
# --lambda_method ctrl_max \
# --poisson_method mutation \
# 2> poisson_res/293T-DdCBE-ND6-All-PD_vs_ctrl_hg38_possion_test.log &

# ------------------------------------------------------------------->>>>>>>>>>

rule find_significant_mpmat:
    input:
        mpmat="../mpmat/{sample}_merge.select.sort_rmchrYM.mpmat.gz",
        bam="../bam/{sample}_final_rmdup.bam"
    output:
        "../poisson_res/{sample}_vs_ctrl_%s.select.pvalue_table" % CTRL_NAME
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa'
    log:
        "../poisson_res/{sample}_vs_ctrl_%s.select.log" % CTRL_NAME
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