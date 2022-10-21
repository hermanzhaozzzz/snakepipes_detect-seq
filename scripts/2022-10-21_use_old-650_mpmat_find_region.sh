SAMPLES="DetectSeq_ATP8-DddA6_REP-1 DetectSeq_ATP8-DddA11_REP-1 DetectSeq_ATP8-DddAwt_REP-1 DetectSeq_JAK2-DddA11_REP-1 DetectSeq_SIRT6-DddA11_REP-1 test"

mkdir -p ../poisson_res_use_old650

for sample in $SAMPLES
do
    echo $sample
    
    sleep 1
    
    pkurun-cns 1 20 python program/detect_seq/find-significant-mpmat.py \
        -p 20 \
        -i ../mpmat/2022-10-21_nature_4_5.1_6_share_650-off-targets.mpmat \
        -o ../poisson_res_use_old650/${sample}_vs_ctrl_Vector-merge_hg38.select.pvalue_table \
        -c /lustre1/chengqiyi_pkuhpc/zhaohn/3.project/2022_DdCBE-3D-Genome_topic/2022-09-30_Detect-seq_batch-1/bam/Vector-merge_hg38_merge_sort_rmdup.MAPQ20.bam \
        -t ../bam/${sample}_final_rmdup.bam \
        -r /lustre1/chengqiyi_pkuhpc/zhaohn/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa \
        --query_mutation_type CT,GA \
        --mpmat_filter_info_col_index -1 \
        --mpmat_block_info_col_index -1 \
        --region_block_mut_num_cutoff 2 \
        --query_mut_min_cutoff 2 \
        --query_mut_max_cutoff 16 \
        --total_mut_max_cutoff 16 \
        --other_mut_max_cutoff 6 \
        --seq_reads_length 150 \
        --lambda_method ctrl_max \
        --poisson_method mutation \
        > ../poisson_res_use_old650/${sample}_vs_ctrl_Vector-merge_hg38.select.log 2>&1
done

