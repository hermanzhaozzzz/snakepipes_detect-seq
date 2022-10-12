##################################################################################################################################
# Data 2020-04-08
# Author Wu Hao
# Pmat files show wth different cutoffs
##################################################################################################################################
rm(list=ls())

#------------------------------------------------------------------------------------------>>>>>>>
# 1. load table and make long format table 
#------------------------------------------------------------------------------------------>>>>>>>
library(ggplot2)
library(dplyr)
library(limma)

#setwd("~/menghw_HD/BE_project/24.APOBEC_mut_RNASeq/")

pmat_cutoffs.df = read.csv("All_samples_filtered_mpmat_diff_cutoffs_count_stats_deal.csv",header = T)

#### need to plot with heatmap !!!
library(ggplot2)
library(RColorBrewer)
RColorBrewer::display.brewer.all()
color_list = rev(colorRampPalette(brewer.pal(9,"RdYlBu"))(30))


# CT plot
df_CT <- subset(pmat_cutoffs.df, pmat_cutoffs.df$Mut_type ==  "CT" & pmat_cutoffs.df$Mut_ratio == 0.1)
df_GA <- subset(pmat_cutoffs.df, pmat_cutoffs.df$Mut_type ==  "GA" & pmat_cutoffs.df$Mut_ratio == 0.1)

ggplot(data=df_CT,aes(x=as.factor(Cover_num),y=as.factor(Mut_num))) +
  geom_tile(aes(fill= log2(Count+1))) +
  geom_text(aes(label= Count),check_overlap = TRUE) +
  scale_fill_gradientn(colours=color_list) + 
  #facet_grid(cols = vars(Pass_num)+vars(Mut_type), rows = vars(Sample)) + 
  facet_grid(Sample ~ Pass_num) +
  labs(title = "All mutation count from mpmat",x="Coverage Number Cutoff",y="Mutation Reads Number Cutoff",colour = "log2(count)") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "./20201112_All_mut_points.heatmap_diff_cutoffs_CT_mutratio0.1.pdf",width =20,height = 10,limitsize = FALSE)

ggplot(data=df_GA,aes(x=as.factor(Cover_num),y=as.factor(Mut_num))) +
  geom_tile(aes(fill= log2(Count+1))) +
  geom_text(aes(label= Count),check_overlap = TRUE) +
  scale_fill_gradientn(colours=color_list) + 
  #facet_grid(cols = vars(Pass_num)+vars(Mut_type), rows = vars(Sample)) + 
  facet_grid(Sample ~ Pass_num) +
  labs(title = "All mutation count from mpmat",x="Coverage Number Cutoff",y="Mutation Reads Number Cutoff",colour = "log2(count)") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "./20201112_All_mut_points.heatmap_diff_cutoffs_GA_mutratio0.1.pdf",width =20,height = 10,limitsize = FALSE)

