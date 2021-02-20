#!/bin/bash
# properties = {"type": "single", "rule": "mpmat_filter", "local": false, "input": ["../mpmat/293T-bat_All-EMX1-PD_rep1_hg38.select.G2A.mpmat"], "output": ["../mpmat_filter/293T-bat_All-EMX1-PD_rep1_hg38.filtered_SMN-3_SCN-6_SMR-0.1_RPN-1_RTN-2.G2A.mpmat"], "wildcards": {"sample": "All-EMX1-PD", "rep": "rep1", "SMN": "3", "SCN": "6", "SMR": "0.1", "RPN": "1", "RTN": "2", "ac": "G", "bd": "A"}, "params": {"ac": "G", "bd": "A", "SiteMutNum": "3", "SiteCoverNum": "6", "SiteMutRatio": "0.1", "RegionPassNum": "1", "RegionToleranceNum": "2"}, "log": [], "threads": 1, "resources": {}, "jobid": 27, "cluster": {}}
 cd /gpfs/user/zhaohuanan/3.project/2021_CBE_variant_and_GBE_topic/20210202_DetectSeq_M2-YE1-33A-call-points_sgRNA-EMX1-VEGFA_redo/snakepipes_Detect-seq && \
PATH='/home/zhaohuanan/zhaohn_HD/miniconda3/bin':$PATH /home/zhaohuanan/zhaohn_HD/miniconda3/bin/python3.8 \
-m snakemake ../mpmat_filter/293T-bat_All-EMX1-PD_rep1_hg38.filtered_SMN-3_SCN-6_SMR-0.1_RPN-1_RTN-2.G2A.mpmat --snakefile /gpfs/user/zhaohuanan/3.project/2021_CBE_variant_and_GBE_topic/20210202_DetectSeq_M2-YE1-33A-call-points_sgRNA-EMX1-VEGFA_redo/snakepipes_Detect-seq/step.02.snk.MpmatFilter.py \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/user/zhaohuanan/3.project/2021_CBE_variant_and_GBE_topic/20210202_DetectSeq_M2-YE1-33A-call-points_sgRNA-EMX1-VEGFA_redo/snakepipes_Detect-seq/.snakemake/tmp._95jc29l ../mpmat/293T-bat_All-EMX1-PD_rep1_hg38.select.G2A.mpmat --latency-wait 60 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules mpmat_filter --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

