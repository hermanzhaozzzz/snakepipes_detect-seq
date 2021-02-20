library(tidyverse)
library(valr)
rm(list = ls())
setwd("~/zhaohn_HD/3.projects/2021_CBE_variant_and_GBE_topic/20210202_DetectSeq_M2-YE1-33A-call-points_sgRNA-EMX1-VEGFA_redo/20210204_all_scored_mpmat/")
# the col name setting
mpmat.standard.colname <- c(
  'chrom',
  'start',
  'end',
  'region_site_num',
  'region_mut_site_num',
  'region_SNP_mut_num',
  'region_site_index',
  'mut_base_num',
  'cover_base_num',
  'mut_ratio',
  'SNP_ann',
  'tandem_info',
  'PASS_state'
)

# load bed files
sgRNA.EMX1.bed.df <- read_tsv(
  file = "bed_20210202_CBE_DetectSeq_de-novo/293T-EMX1-off-target_plus_on-target.sorted.bed",
  col_names = c("chrom","start","end","name","score","strand")
)
sgRNA.VEGFA.bed.df <- read_tsv(
  file = "bed_20210202_CBE_DetectSeq_de-novo/293T-VEGFA-off-target_plus_on-target.sorted.bed",
  col_names = c("chrom","start","end","name","score","strand")
)








# load all-mpmat file and deal with it and write out
file.list <- list.files("mpmat_filter_merge/")
for (mpmat_name in file.list) {
  MPMAT.FILE <- stringr::str_c("mpmat_filter_merge/",mpmat_name)
  if (stringr::str_detect(MPMAT.FILE,".sort.mpmat")) {
  
    mpmat.df <- read_tsv(
      file = MPMAT.FILE,
      col_names = mpmat.standard.colname
    )
    

    if (stringr::str_detect(MPMAT.FILE,"-EMX1-")) {
      intersect.df <- bed_intersect(x = mpmat.df, sgRNA.EMX1.bed.df)
    } else if (stringr::str_detect(MPMAT.FILE,"-VEGFA-")) {
      intersect.df <- bed_intersect(x = mpmat.df, sgRNA.VEGFA.bed.df)
    } else {
      intersect.df <- NULL
    }

    colnames(intersect.df) <- c(mpmat.standard.colname, c("bed.start", "bed.end", "bed.name", "bed.score", "bed.strand", "overlap.info"))

    write_tsv(intersect.df, file = stringr::str_replace(MPMAT.FILE, ".sort.mpmat",".sort.intersected.mpmat"), col_names = FALSE)
  }
}


