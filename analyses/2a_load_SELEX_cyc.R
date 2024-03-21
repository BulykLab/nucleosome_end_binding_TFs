library(rstatix)
library(matrixStats)

source("utils.R")
source("config.R")

tf_list = read.table(SELEX_tf_list_file)$V1
data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")
k=7


complete_df = data.frame(matrix(nrow = 0, ncol=6))
for (tf_name in tf_list){ 
  tf_dir = paste0(data_dir, tf_name,"/cyc/")
  
  if (!(dir.exists(paste0(tf_dir, "sample_left/")) && 
        dir.exists(paste0(tf_dir, "control_left/")) &&
        dir.exists(paste0(tf_dir, "ht_left/")))){
    next
  }
  
  tf_sample_left = load_cyc_from_dir(paste0(tf_dir, "sample_left/"))
  tf_sample_right_rc = load_cyc_from_dir(paste0(tf_dir, "sample_right_rc/"))
  tf_sample = cbind(tf_sample_left, tf_sample_right_rc)
  tf_sample_plot_df = build_plot_df(tf_sample, tf_name, "NCAP-SELEX")
  
  tf_control_left = load_cyc_from_dir(paste0(tf_dir, "control_left/"))
  tf_control_right_rc = load_cyc_from_dir(paste0(tf_dir, "control_right_rc/"))
  tf_control = cbind(tf_control_left, tf_control_right_rc)
  tf_control_plot_df = build_plot_df(tf_control, tf_name, "Nucleosome-SELEX")
  
  tf_ht_left = load_cyc_from_dir(paste0(tf_dir, "ht_left/"))
  tf_ht_right_rc = load_cyc_from_dir(paste0(tf_dir, "ht_right_rc/"))
  tf_ht = cbind(tf_ht_left, tf_ht_right_rc)
  tf_ht_plot_df = build_plot_df(tf_ht, tf_name,"HT-SELEX")
  
  plot_df = rbind(tf_sample_plot_df, tf_control_plot_df, tf_ht_plot_df)
  complete_df = rbind(complete_df, plot_df)
}

write.csv(complete_df, "../data/cyc_results.csv")