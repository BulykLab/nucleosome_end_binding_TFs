source("utils.R")
library(ShortRead)
library(glue)


tf_list = read.table("/home/xil493/capstone/data/all_lig147_tfs.txt")$V1
data_dir <- "/n/data2/bch/medicine/bulyk/Katrina/capstone/data/lig147/"
k=7

control_dir <- paste0(data_dir, "Nucleosome-SELEX__Bound/Nucleosome-SELEX_cycle-4/")
control_fastq = readFastq(get_clean_file(control_dir))
control_seqs <- sread(control_fastq)

cyc_script = c()


for (tf_name in tf_list){
  tf_dir <- paste0(data_dir, tf_name, "/")
  tf_sample_dir <- paste0(tf_dir,"NCAP-SELEX_cycle-4/")
  if (!dir.exists(tf_sample_dir)){
    next
  }
  
  if (!file.exists(paste0(tf_sample_dir,"final_enriched_kmers.txt"))){
    next
  }
  
  tf_cyc_dir <- paste0(tf_dir,"cyc/")
  if (!file.exists(tf_cyc_dir)) {
    dir.create(tf_cyc_dir)
  }
  
  sample_fastq = readFastq(get_clean_file(tf_sample_dir))
  sample_seqs <- sread(sample_fastq)
  
  
  
  motifs = read.table(paste0(tf_sample_dir,"final_enriched_kmers.txt"))$V1
  sample_positions = to_get_many_kmer_match_reads(motifs, sample_seqs)
  control_positions = to_get_many_kmer_match_reads(motifs, control_seqs)
  pos_col = ncol(sample_positions)
  
  
  left_sample = filter_seq_indices(sample_positions,1:4)
  right_sample = filter_seq_indices(sample_positions, (pos_col-3):pos_col)
  if (1000> length(left_sample)+length(right_sample) && length(left_sample)+length(right_sample) >= 500){
    left_sample_seqs = sample_fastq[left_sample,]
    right_sample_seqs = sample_fastq[right_sample, ]
    right_sample_rc_seqs = reverseComplement(right_sample_seqs)
    writeFasta(left_sample_seqs, paste0(tf_cyc_dir, "sample_left.fasta"))
    writeFasta(right_sample_rc_seqs, paste0(tf_cyc_dir, 
                                            "sample_right_rc.fasta"))
    cyc_script = append(cyc_script, 
                        glue('sbatch run_cycp_tmp.sh {paste0(tf_cyc_dir, "sample_left.fasta")} {paste0(tf_cyc_dir, "sample_left/")}' ))
    cyc_script = append(cyc_script, 
                        glue('sbatch run_cycp_tmp.sh {paste0(tf_cyc_dir, "sample_right_rc.fasta")} {paste0(tf_cyc_dir, "sample_right_rc/")}' ))
  }
  

  left_control = filter_seq_indices(control_positions,1:4)
  right_control = filter_seq_indices(control_positions, (pos_col-3):pos_col)
  if (length(left_control)+length(right_control) >= 500){
    left_control_seqs = control_fastq[left_control,]
    right_control_seqs = control_fastq[right_control, ]
    right_control_rc_seqs = reverseComplement(right_control_seqs)
    writeFasta(left_control_seqs, paste0(tf_cyc_dir, "control_left.fasta"))
    writeFasta(right_control_rc_seqs, paste0(tf_cyc_dir,
                                             "control_right_rc.fasta"))
    cyc_script = append(cyc_script,
                        glue('sbatch run_cycp_tmp.sh {paste0(tf_cyc_dir, "control_left.fasta")} {paste0(tf_cyc_dir, "control_left/")}' ))
    cyc_script = append(cyc_script,
                        glue('sbatch run_cycp_tmp.sh {paste0(tf_cyc_dir, "control_right_rc.fasta")} {paste0(tf_cyc_dir, "control_right_rc/")}' ))

  }

  tf_ht_dir <- paste0(tf_dir,"HT-SELEX_cycle-4/")
  if (!dir.exists(tf_ht_dir)){
    next
  }



  ht_fastq = readFastq(get_clean_file(tf_ht_dir))
  ht_seqs <- sread(ht_fastq)
  ht_positions = to_get_many_kmer_match_reads(motifs, ht_seqs)


  left_ht = filter_seq_indices(ht_positions,1:4)
  right_ht = filter_seq_indices(ht_positions, (pos_col-3):pos_col)

  if (length(left_ht)+length(right_ht) >= 500){
    left_ht_seqs = ht_fastq[left_ht,]
    right_ht_seqs = ht_fastq[right_ht, ]
    right_ht_rc_seqs = reverseComplement(right_ht_seqs)
    writeFasta(left_ht_seqs, paste0(tf_cyc_dir, "ht_left.fasta"))
    writeFasta(right_ht_rc_seqs, paste0(tf_cyc_dir, "ht_right_rc.fasta"))
    cyc_script = append(cyc_script,
                        glue('sbatch run_cycp_tmp.sh {paste0(tf_cyc_dir, "ht_left.fasta")} {paste0(tf_cyc_dir, "ht_left/")}' ))
    cyc_script = append(cyc_script,
                        glue('sbatch run_cycp_tmp.sh {paste0(tf_cyc_dir, "ht_right_rc.fasta")} {paste0(tf_cyc_dir, "ht_right_rc/")}' ))

  }
  
}
# 
write.table(cyc_script, "selex_cyc_script_sup.sh",
            quote=F,row.names = F, col.names = F)
