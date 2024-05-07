library(rstatix)
library(matrixStats)

source("utils.R")
source("config.R")

tf_list = read.table(SELEX_tf_list_file)$V1
data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")
k=7

# Store cyclizability predictions
complete_df = data.frame(matrix(nrow = 0, ncol=6))
internal_df = data.frame(matrix(nrow = 0, ncol=6))
for (tf_name in tf_list){ 
  tf_dir = paste0(data_dir, tf_name,"/cyc/")
  
  # Only consider TFs with enough reads for statistical analyses
  if (!(file.exists(paste0(tf_dir, "sample_left_primer.csv")) && 
        file.exists(paste0(tf_dir, "control_left_primer.csv")) &&
        file.exists(paste0(tf_dir, "ht_left_primer.csv")))){
    next
  }
  
  # Load cyc predictions and build a data frame for plots
  tf_sample_left = read.csv(paste0(tf_dir, "sample_left_primer.csv"))
  tf_sample_right_rc = read.csv(paste0(tf_dir, "sample_right_rc_primer.csv"))
  tf_sample = rbind(tf_sample_left, tf_sample_right_rc)
  tf_sample_plot_df = build_plot_df(t(tf_sample[,2:ncol(tf_sample)]), tf_name, "NCAP-SELEX")
  
  tf_control_left = read.csv(paste0(tf_dir, "control_left_primer.csv"))
  tf_control_right_rc = read.csv(paste0(tf_dir, "control_right_rc_primer.csv"))
  tf_control = rbind(tf_control_left, tf_control_right_rc)
  tf_control_plot_df = build_plot_df(t(tf_control[,2:ncol(tf_control)]), tf_name, "Nucleosome-SELEX")
  
  tf_ht_left = read.csv(paste0(tf_dir, "ht_left_primer.csv"))
  tf_ht_right_rc = read.csv(paste0(tf_dir, "ht_right_rc_primer.csv"))
  tf_ht = rbind(tf_ht_left, tf_ht_right_rc)
  tf_ht_plot_df = build_plot_df(tf_ht[, 2:ncol(tf_ht)], tf_name,"HT-SELEX")
  
  plot_df = rbind(tf_sample_plot_df, tf_control_plot_df, tf_ht_plot_df)
  complete_df = rbind(complete_df, plot_df)
  
  
  if (!file.exists(paste0(tf_cyc_dir, "sample_left_1_5_primer.csv"))){
    next
  }
  
  for (internal_pos in seq(1, 20, 5)){
    left_pos = (internal_pos):(internal_pos+4)
    right_pos = (pos_col-internal_pos-3):(pos_col-internal_pos+1)
    
    left_sample_cyc_csv = paste0(tf_cyc_dir, "sample_left_",min(left_pos),
                                 "_", max(left_pos) ,"_primer.csv")
    right_sample_rc_cyc_csv = paste0(tf_cyc_dir, "sample_right_rc_",min(right_pos),
                                     "_", max(right_pos) ,"_primer.csv")
    left_cyc = read.csv(left_sample_cyc_csv)
    right_cyc = read.csv(right_sample_rc_cyc_csv)
    sub_plot_df = build_plot_df(t(rbind(left_cyc, right_cyc)[,2:ncol(left_cyc)]), tf_name, paste0(internal_pos+24))
    internal_df = rbind(internal_df, 
                    sub_plot_df
    )
  }
  
}

write.csv(complete_df, SELEX_cyc_result_csv)
write.csv(internal_df, SELEX_internal_cyc_result_csv)