SAMPLES = ['DD-gfp-pd', 'Vec-ab-pd-rep1']


THREADS = '24'


# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
BWA_HG38_IDX = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.bwa_index"
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"





# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
# make sure the programs used below are in you PATH
# get the application path

PICARD = "/home/zhaohuanan/zhaohn_HD/1.apps/picard/picard.jar"

with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
    print('PATH cutadapt:', CUTADAPT)
with os.popen("which bismark") as path:
    BISMARK = path.read().strip()
    print('PATH bismark:', BISMARK)
with os.popen("which bwa") as path:
    BWA = path.read().strip()
    print('PATH bwa:', BWA)
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
    print('PATH samtools:', SAMTOOLS)
with os.popen("which java") as path:
    JAVA = path.read().strip()
    print('PATH java:', JAVA)
with os.popen("which sambamba") as path:
    SAMBAMBA = path.read().strip()







rule all:
    input:
        expand("../fix.fastq/293T-{sample}_R1_cutadapt.fq.gz", sample=SAMPLES),
        expand("../fix.fastq/293T-{sample}_R2_cutadapt.fq.gz", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38.sam", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38.bam", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38_sort.bam", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.bam", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.FilterClip.bam", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup_MAPQ20.bam", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38.mapping_stats", sample=SAMPLES),
        expand("../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup_MAPQ20.mapping_stats", sample=SAMPLES),
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# FASTQ quality control
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule cutadapt_TrueSeq:
    input:
        "../fastq/{sample}_combined_R1.fastq.gz",
        "../fastq/{sample}_combined_R2.fastq.gz"
    output:
        "../fix.fastq/293T-{sample}_R1_cutadapt.fq.gz",
        "../fix.fastq/293T-{sample}_R2_cutadapt.fq.gz"
    log:
        "../fix.fastq/293T-{sample}_cutadapt.log"
    shell:
        """
        {CUTADAPT} -j {THREADS} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 55 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
        """         

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# hg38 mapping 
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule Bwa_mapping:
    input:
        "../fix.fastq/293T-{sample}_R1_cutadapt.fq.gz",
        "../fix.fastq/293T-{sample}_R2_cutadapt.fq.gz"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38.sam"
    log:
        "../bam.bwa/293T-{sample}_bwa_hg38.log"
    params:
        "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
    shell:
        "{BWA} mem {BWA_HG38_IDX} {input[0]} {input[1]} -t {THREADS} -M -R {params} > {output} 2>{log}"
rule Sam_to_bam_file:
    input:
        "../bam.bwa/293T-{sample}_bwa_hg38.sam"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38.bam"
    shell:
        "{SAMBAMBA} view -t {THREADS} -h -f bam -S -o {output} {input}"
rule BAM_sort_by_position:
    input:
        "../bam.bwa/293T-{sample}_bwa_hg38.bam"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort.bam"
    shell:
        "{SAMBAMBA} sort -t 24 -m 50G -o {output} --tmpdir=../bam.bwa {input}"

rule BAM_rmdup:
    input:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort.bam"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.bam",
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.matrix"
    log:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.log"
    shell:
        "{JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads={THREADS} -jar {PICARD} MarkDuplicates I={input} O={output[0]} M={output[1]} ASO=coordinate REMOVE_DUPLICATES=true 2>{log}"
# 268 = read unmapped + mate unmapped +  not primary alignment
# 1 = read paired
rule BAM_filter_clip:
    input:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.bam"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.FilterClip.bam"
    shell:
        """{SAMTOOLS} view -h -f 1 -F 268 {input} | samclip --ref {HG38_FA} --max 3 --progress 0 | awk 'function abs(v) {{return v < 0 ? -v : v}} $1~"@" || ($7 == "=" && abs($9) <= 2500 ) {{print $0}}' | samtools view -hb > {output}"""

# rule BAM_samtools_index:
#     input:
#         "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.bam",
#     output:
#         "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.bam.bai",
#     shell:
#         "{SAMTOOLS} index -@ {THREADS} {input} {output}"



rule MAPQ20_select_bam:
    input:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup.FilterClip.bam"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup_MAPQ20.bam"
    shell:
        '{SAMBAMBA} view -t {THREADS} -h -f bam -F "mapping_quality >= 20" -o {output} {input}'
rule get_mapping_stats:
    input:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort.bam"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38.mapping_stats"
    shell:
        "{SAMTOOLS} stats -@ {THREADS} --remove-overlaps --reference {HG38_FA} {input} > {output}"

rule final_get_mapping_stats:
    input:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup_MAPQ20.bam"
    output:
        "../bam.bwa/293T-{sample}_bwa_hg38_sort_rmdup_MAPQ20.mapping_stats"
    shell:
        "{SAMTOOLS} stats -@ {THREADS} --remove-overlaps --reference {HG38_FA} {input} > {output}"
