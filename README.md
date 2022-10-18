# snakepipes_detect-seq

## envs
```shell
conda install python=3.10 pysam pysamstats biopython
pip install bioat --upgrade
```

# conda environment setting
```
mamba create -n snakepipes_detect-seq python=2.7 bedtools=2.27 samtools biopython=1.72 pysam statsmodels pandas
```
# snk.1.bam2mpmat.py
这一步可以尽可能多的并行任务，如24CPU cores并行24任务
```
snakemake --jobs 24 --snakefile snk.1.bam2mpmat.py -pr
```
# snk.2.mpmat_merge.py
这一步一个node尽量并行4个以下的任务，bedtools sort容易内存崩溃，跑任务的时候多注意
且注意bedtools version 2.27，不然可能出问题
```
snakemake --jobs 4 --snakefile snk.2.mpmat_merge.py -pr
```


```
   446070 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
   453410 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
   261000 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
   265020 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.1_RPN-2_RTN-False.G2A.mpmat
   473592 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
   480956 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
   288521 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
   292566 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.1_RPN-2_RTN-False.G2A.mpmat
     5754 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
     5777 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
     2093 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
     2118 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.1_RPN-2_RTN-False.G2A.mpmat
     5820 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
     5853 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
     2159 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
     2193 293T-bat_GBEmini-AP-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.1_RPN-2_RTN-False.G2A.mpmat
    54403 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
    54698 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
    46667 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
    47060 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-10_SMR-0.1_RPN-2_RTN-False.G2A.mpmat
   142924 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
   143161 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
   135186 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
   135524 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-2_SCN-5_SMR-0.1_RPN-2_RTN-False.G2A.mpmat
     2305 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
     2296 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
      662 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
      734 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-10_SMR-0.1_RPN-2_RTN-False.G2A.mpmat
     2398 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.01_RPN-2_RTN-False.C2T.mpmat
     2378 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.01_RPN-2_RTN-False.G2A.mpmat
      755 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.1_RPN-2_RTN-False.C2T.mpmat
      816 293T-bat_GBEmini-dU-RNF2-All-PD_rep1_hg38.filtered_SMN-4_SCN-5_SMR-0.1_RPN-2_RTN-False.G2A.mpmat


```




# snk.3.poisson-test_select_plot-art.py
注意biopython版本一定是1.72