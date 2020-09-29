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