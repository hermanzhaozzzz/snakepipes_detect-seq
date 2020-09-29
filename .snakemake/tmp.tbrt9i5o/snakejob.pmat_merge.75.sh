#!/bin/bash
# properties = {"type": "single", "rule": "pmat_merge", "local": false, "input": ["../mpileup_pmat_bmat/293T-bat_VEGFA-All-PD_rep1_hg38.MAPQ20.select.T.pmat"], "output": ["../mpmat/293T-bat_VEGFA-All-PD_rep1_hg38.MAPQ20.select.T2C.mpmat"], "wildcards": {"sample": "VEGFA-All-PD", "rep": "rep1", "from_base": "T", "to_base": "C"}, "params": {"from_base": "T", "to_base": "C"}, "log": ["../mpmat/293T-bat_VEGFA-All-PD_rep1_hg38.MAPQ20.select.T2C.mpmat.log"], "threads": 1, "resources": {}, "jobid": 75, "cluster": {}}
 cd /gpfs/user/zhaohuanan/3.project/10.2020-09_0923_Practice_for_detect-seq-data/snakepipes_detect-seq_fastq_bam_plot && \
/home/zhaohuanan/miniconda3/bin/python \
-m snakemake ../mpmat/293T-bat_VEGFA-All-PD_rep1_hg38.MAPQ20.select.T2C.mpmat --snakefile /gpfs/user/zhaohuanan/3.project/10.2020-09_0923_Practice_for_detect-seq-data/snakepipes_detect-seq_fastq_bam_plot/snk.1.bam2mpmat.py \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/user/zhaohuanan/3.project/10.2020-09_0923_Practice_for_detect-seq-data/snakepipes_detect-seq_fastq_bam_plot/.snakemake/tmp.tbrt9i5o ../mpileup_pmat_bmat/293T-bat_VEGFA-All-PD_rep1_hg38.MAPQ20.select.T.pmat --latency-wait 60 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules pmat_merge --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

