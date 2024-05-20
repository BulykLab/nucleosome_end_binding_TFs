source("utils.R")
source("config.R")

tf_list = read.table(SELEX_tf_list_file)$V1
data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")



left_res = data.frame(matrix(nrow=0, ncol=4))
right_res = data.frame(matrix(nrow=0, ncol=4))

for (tf_name in tf_list){
  tf_dir <- paste0(data_dir, tf_name, "/")
  tf_cyc_dir <- paste0(tf_dir,"cyc/")
  if (!dir.exists(tf_cyc_dir)) {
    next
  }
  if ((!file.exists(paste0(tf_cyc_dir, "sample_left_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "sample_right_rc_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "control_left_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "control_right_rc_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "ht_left_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "ht_right_rc_primer.fasta")))){
    next
  }
  
  sample_left = readFasta(paste0(tf_cyc_dir, "sample_left_primer.fasta"))
  sample_right_rc = readFasta(paste0(tf_cyc_dir, "sample_right_rc_primer.fasta"))
  sample = DNAStringSet(c(sread(sample_left),sread(sample_right_rc)))
  
  control_left = readFasta(paste0(tf_cyc_dir, "control_left_primer.fasta"))
  control_right_rc = readFasta(paste0(tf_cyc_dir, "control_right_rc_primer.fasta"))
  control = DNAStringSet(c(sread(control_left),sread(control_right_rc)))
  
  ht_left = readFasta(paste0(tf_cyc_dir, "ht_left_primer.fasta"))
  ht_right_rc = readFasta(paste0(tf_cyc_dir, "ht_right_rc_primer.fasta"))
  ht = DNAStringSet(c(sread(ht_left),sread(ht_right_rc)))
  
  tf_fft_dir = paste0(tf_dir,"/fft/")
  if (!dir.exists(tf_fft_dir)){
    dir.create(tf_fft_dir)
  }
  
  left_freq=get_nuc_spectrum_x(left_region)
  right_freq=get_nuc_spectrum_x(right_region)
  
  left_sample_density= get_nuc_spectrum_from_sample(nuc, sample, left_region)
  right_sample_density= get_nuc_spectrum_from_sample(nuc, sample, right_region)
  
  left_control_density= get_nuc_spectrum_from_sample(nuc, control, left_region)
  right_control_density= get_nuc_spectrum_from_sample(nuc, control, right_region)
  
  left_ht_density= get_nuc_spectrum_from_sample(nuc, ht, left_region)
  right_ht_density= get_nuc_spectrum_from_sample(nuc, ht, right_region)
  
  left_res = rbind(left_res, 
                   data.frame(tf=tf_name, 
                                  freq=left_freq,
                                  density=left_sample_density,
                                  group="NCAP-SELEX"
                                  ),
                   data.frame(tf=tf_name, 
                              freq=left_freq,
                              density=left_control_density,
                              group="Nucleosome-SELEX"
                   ),
                   data.frame(tf=tf_name, 
                              freq=left_freq,
                              density=left_ht_density,
                              group="HT-SELEX"
                   ))
  right_res = rbind(right_res, 
                   data.frame(tf=tf_name, 
                              freq=right_freq,
                              density=right_sample_density,
                              group="NCAP-SELEX"
                   ),
                   data.frame(tf=tf_name, 
                              freq=right_freq,
                              density=right_control_density,
                              group="Nucleosome-SELEX"
                   ),
                   data.frame(tf=tf_name, 
                              freq=right_freq,
                              density=right_ht_density,
                              group="HT-SELEX"
                   ))
}

write.csv(left_res, SELEX_fft_left)
write.csv(right_res, SELEX_fft_right)
