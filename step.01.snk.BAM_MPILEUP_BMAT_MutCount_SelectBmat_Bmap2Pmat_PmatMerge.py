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




import os
PATH = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_detect-seq/bin/"
PYTHON2 = os.path.join(PATH, 'python')
BEDTOOLS = os.path.join(PATH, 'bedtools')
SAMTOOLS = os.path.join(PATH, 'samtools')
GENOME = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
GENOME_FAI = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.fai"
DB_SNP_HSA = "/home/zhaohuanan/zhaohn_HD/2.database/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"



SELECT_BASES_FROM = [
    'A',
    'T',
    'C',
    'G'
]
SELECT_BASES_TO = [
    'A',
    'T',
    'C',
    'G'
]
rule all:
    input:
        expand("../bam/293T-bat_{sample}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam",sample=SAMPLES,rep=REP),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup",sample=SAMPLES,rep=REP),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat",sample=SAMPLES,rep=REP),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.bmat",sample=SAMPLES,rep=REP,fbase=SELECT_BASES_FROM),
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.pmat",sample=SAMPLES,rep=REP,fbase=SELECT_BASES_FROM),
        expand("../Split_Bmat_Json/293T-bat_{sample}_{rep}_hg38.MAPQ20.json",sample=SAMPLES,rep=REP),
#         expand("../bmat_info_tsv/293T-bat_{sample}_{rep}_hg38.MAPQ20.tsv",sample=SAMPLES,rep=REP),
        expand("../mpmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}2{tbase}.mpmat",sample=SAMPLES,rep=REP,fbase=SELECT_BASES_FROM,tbase=SELECT_BASES_TO)
rule bam2mpileup:
    input:
        "../bam/293T-bat_{sample}_{rep}_bwa_hg38_sort_rmdup_MAPQ20.bam"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup"
    log:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup.log"
    shell:
#         """
#         touch {output}
#         """
        """
        {SAMTOOLS} mpileup -q 20 -Q 20 \
        --reference {GENOME} \
        -o {output} \
        {input} >{log} 2>&1
        """
rule mpileup2bmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.mpileup"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat"
    log:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat.log"
    shell:
#         """
#         touch {output}
#         """
        """
        {PYTHON2} ./program/parse-mpileup-V04.py \
        -i {input} \
        -o {output} \
        -p 24 -n 0 >{log} 2>&1
        """
# rule bmat_mutation_count:
#     input:
#         "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat"
#     output:
#         "../bmat_info_tsv/293T-bat_{sample}_{rep}_hg38.MAPQ20.tsv"
#     shell:
#         """
#         {PYTHON2} ./program/bmat_mutation_count.py -i {input} -o {output} --header False
#         """
rule select_bmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.bmat"
    params:
        awk = """{print $0}""",
        fbase = "{fbase}"
    shell:
#         """
#         touch {output}
#         """  
        """
        awk '$3 == "{params.fbase}" {params.awk}' {input} > {output} 
        """ 
rule select_bmat2pmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.bmat"
    output:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.pmat"
    shell:
#         """
#         touch {output}
#         """
        """
        {PYTHON2} ./program/bmat2pmat-V02.py \
        -i {input} \
        -o {output} \
        --InHeader False \
        --InLikeBED False \
        --OutHeader True
        """
rule pmat_merge:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}.pmat"
    output:
        "../mpmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}2{tbase}.mpmat"
    log:
        "../mpmat/293T-bat_{sample}_{rep}_hg38.select.{fbase}2{tbase}.mpmat.log"
    params:
        FB = '{fbase}',
        TB = '{tbase}'
    shell:
        """
        touch {output}
        """
        """
        FB={params.FB}
        TB={params.TB}
        if [ $FB = $TB ]; then
        echo "from base == to base, skip..." > {log}
        touch {output}
        else
        {PYTHON2} ./program/pmat-merge-V04.py \
        -f {params.FB} -t {params.TB} \
        -r {GENOME} \
        --OutHeader False \
        -i {input} \
        -o {output} \
        -d 50 -D 100 \
        --NoMutNumCutoff 2 \
        --OmitTandemNumCutoff 2 \
        --SNP {DB_SNP_HSA} >{log} 2>&1
        fi
        """
rule Split_Bmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat"
    output:
        "../Split_Bmat_Json/293T-bat_{sample}_{rep}_hg38.MAPQ20.json"
    params:
        output_dir = "../Split_Bmat_Json"
    shell:
#         """
#         touch {output}
#         """
        """
        {PYTHON2} ./program/split-bmat-V01.py \
        -i {input} \
        -o {output} \
        -j {output}
        """

