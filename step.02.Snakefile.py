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
# rule all
# ------------------------------------------------------------------->>>>>>>>>>
rule all:
    input:
        expand("../fastq/{sample}_%s.fastq.gz" % READ[0],sample=SAMPLES),
        expand("../bam/{sample}_final_rmdup.bam",sample=SAMPLES),
        expand("../bam/{sample}_final_rmdup.bam.flagstats.tsv",sample=SAMPLES),
        expand("../mpileup/{sample}_final_rmdup.mpileup.gz",sample=SAMPLES),
        expand("../bmat/{sample}_final_rmdup.bmat.gz",sample=SAMPLES),
        expand("../pmat/{sample}_final_rmdup.pmat.gz", sample=SAMPLES),
        expand("../pmat/{sample}_CT.mpmat.gz", sample=SAMPLES),
        expand("../pmat/{sample}_GA.mpmat.gz", sample=SAMPLES)

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
        {BWA} mem -t {THREAD} -M -R {params.awk} \
            {params.ref} {input.R1} {input.R2}  > {output} 2>{log}
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
        raw_bam="../bam/{sample}_hisat3n.bam",
        remapped_bam="../bam/{sample}_remapped.bam"
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
# to see mapping info use flagstats
# ------------------------------------------------------------------->>>>>>>>>>
rule flagstats:
    input: "../bam/{sample}_final_rmdup.bam"
    output: "../bam/{sample}_final_rmdup.bam.flagstats.tsv"
    shell:
        """
        {SAMTOOLS} flagstats -@ {THREAD} -O tsv {input} > {output}
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
rule bam2mpileup:
    input:
        "../bam/{sample}_final_rmdup.bam"
    output:
        "../mpileup/{sample}_final_rmdup.mpileup.gz"
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa',
    shell:
        """
        {SAMTOOLS} mpileup -q 20 -Q 20 --reference {params.ref} {input} | gzip > {output}
        """

rule mpileup2bmat:
    input:
        "../mpileup/{sample}_final_rmdup.mpileup.gz"
    output:
        "../bmat/{sample}_final_rmdup.bmat.gz"
    log:
        "../bmat/{sample}_final_rmdup.bmat.log"
    shell:
        """
        python program/detect_seq/parse-mpileup.py -i {input} -o {output} -p {THREAD} -n 0 > {log} 2>&1
        """

rule select_bmat2pmat:
    input:
        "../bmat/{sample}_final_rmdup.bmat.gz"
    output:
        "../pmat/{sample}_final_rmdup.pmat.gz"
    shell:
        """
        python program/detect_seq/bmat2pmat.py \
        -i {input} \
        -o {output} \
        --InHeader False \
        --InLikeBED False \
        --OutHeader False
        """
# ------------------------------------------------------------------->>>>>>>>>>
# tracing tandem mutation signals
# merge tandem mutation signals and sort
# ------------------------------------------------------------------->>>>>>>>>>
rule pmat_merge_1:
    input:
        "../pmat/{sample}_final_rmdup.pmat.gz"
    output:
        "../pmat/{sample}_CT.mpmat.gz"
        
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa',
        f='C',
        t='T'
    log:
        "../pmat/{sample}_CT.mpmat.gz.log"

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
            --Output {output} \
            --SNP {SNP_LIST} > {log} 2>&1
        """
rule pmat_merge_2:
    input:
        "../pmat/{sample}_final_rmdup.pmat.gz"
    output:
        "../pmat/{sample}_GA.mpmat.gz"
        
    params:
        ref=GENOME_HISAT3N_INDEX.split('.fa.')[0] + '.fa',
        f='G',
        t='A'
    log:
        "../pmat/{sample}_GA.mpmat.gz.log"

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
            --Output {output} \
            --SNP {SNP_LIST} > {log} 2>&1
        """
