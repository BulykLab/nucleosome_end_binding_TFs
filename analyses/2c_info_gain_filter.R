source("utils.R")
source("config.R")

info_gain = read.csv("../data/info_gain.txt",row.names = 1)
# Filter TFs by HT information gain
thres_tfs = info_gain %>% filter(Experiment=="HT-SELEX_cycle-4") %>% 
  filter(InfoGain>=0.15) %>% select(TF)

# Keep TFs with 500 reads for statistical analyses
slope_avg_plot_df = slope_avg_plot_df %>% filter(tf %in% thres_tfs$TF)
low_read_tfs = slope_avg_plot_df %>% filter(n<500)
slope_avg_plot_df = slope_avg_plot_df %>% filter(! tf %in% low_read_tfs$tf)

write.csv(unique(slope_avg_plot_df$tf), "../data/ht_thres_TF.csv", 
          row.names = F)