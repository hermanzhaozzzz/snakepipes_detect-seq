library(tidyverse)
rm(list=ls())

df = read_csv("./20211014-Fig1f.source_data.full_info.csv")
colnames(df)
info.df.select_dist <- df %>% 
  select(sgRNA_name, distance_to_sapce, site_name, rep_name, genome_base, mut_percent, Detect.A, Detect.T, Detect.total)

info.df.select_dist <- info.df.select_dist %>% filter(
  genome_base=="C" & Detect.T >= 10 | genome_base=="G" & Detect.A >= 10
)


info.ratio.df <- info.df.select_dist %>% mutate(
  target_seq.ratio = mut_percent,
  detect_seq.score = ifelse(
    genome_base=="C", 
    (Detect.T/Detect.total)^2 * Detect.T, 
    (Detect.A/Detect.total)^2 * Detect.A
  )
  ) %>% group_by(sgRNA_name, distance_to_sapce, site_name) %>% summarise(
  target_seq.ratio.mean = mean(target_seq.ratio),
  detect_seq.score.mean = mean(detect_seq.score)
)

df <- info.ratio.df
df$target_seq.ratio.mean <- df$target_seq.ratio.mean / 100


pearson_r = cor(log10(df$target_seq.ratio.mean), log10(df$detect_seq.score.mean), method = "pearson")
# pvalue:
# cor.test(log10(df$target_seq.ratio.mean), log10(df$detect_seq.score.mean), method = "pearson")
info.ratio.df.select <- df
# Merge TargetSeq
ggplot(data = info.ratio.df.select) +
  geom_point(mapping = aes(x = target_seq.ratio.mean, y=detect_seq.score.mean, label=site_name, colour = sgRNA_name)) +
  geom_smooth(mapping = aes(x = target_seq.ratio.mean, y=detect_seq.score.mean), method = "lm") +
  scale_color_manual(values=c("#457b9d","#ff9f1c","#e36414")) +
  scale_x_log10(
    name="Mutation ratio by targeted deep sequencing",
    limits=c(0.0001,1),
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
    labels = scales::percent,
    ) +
  scale_y_log10(
    name="Detect-seq Score",
    limits=c(0.01, 1200),
    breaks = c(0.01, 0.1, 1, 10, 100, 1000),
    labels = c(0.01, 0.1, 1, 10, 100, 1000),
  ) +
  ggtitle("Targeted deep sequencing v.s. Detect-seq") +
  theme_classic() + 
  # geom_abline(intercept = 0, slope = 1,linetype="dashed",color="#333333") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text",x=0.0002,y=500,hjust = 0,parse=TRUE,colour="black", label=sprintf("r == %.3f",pearson_r)) +
  annotate("text",x=0.0002,y=300,hjust = 0,parse=TRUE,colour="black", label="p-value < 2.2e-16")

ggsave(filename = "2021-10-14-293T-Merge-Target_seq__vs__Detect_seq.ScatterPlot.pdf",width = 5.4,height = 4)
