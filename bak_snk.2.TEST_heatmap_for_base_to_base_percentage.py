SAMPLES = [
    "GBEmini-AP-RNF2-All-PD",
    "GBEmini-dU-RNF2-All-PD"
]

REP = [
    'rep1'
]

CHR = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']

rule all:
    input:
        expand("../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat",sample=SAMPLES,rep=REP),
        expand("../mpileup_pmat_bmat/split_chr/293T-bat_{sample}_{rep}_hg38.MAPQ20_chr{chr_}.bmat",sample=SAMPLES,rep=REP,chr_=CHR)
rule split_bmat:
    input:
        "../mpileup_pmat_bmat/293T-bat_{sample}_{rep}_hg38.MAPQ20.bmat"
    output:
        "../mpileup_pmat_bmat/split_chr/293T-bat_{sample}_{rep}_hg38.MAPQ20_chr{chr_}.bmat"
    params:
        """'{if($1=="chr{chr_}")print $0}'"""
    shell:
        """sed '1d' {input} | awk -F '\t' {params} > {output}"""
        