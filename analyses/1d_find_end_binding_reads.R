source("utils.R")
source("config.R")

tf_list = read.table(SELEX_tf_list_file)$V1
data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")
k=7

# Load Nucleosome SELEX cycle 4 as control library with only nucleosome binding
control_dir <- paste0(data_dir, 
                      "Nucleosome-SELEX__Bound/Nucleosome-SELEX_cycle-4/")
control_fastq = readFastq(get_clean_file(control_dir))
control_seqs <- sread(control_fastq)

cyc_script = c() # Store cyc commands for running DNA cycP on slurm

perc_ends = c() # Store percentage of enriched k-mers on random region ends
good_tfs = c() # Store TFs with at least 500 end binding reads 

for (tf_name in tf_list){
  tf_dir <- paste0(data_dir, tf_name, "/")
  
  # Ignore TF if not both NCAP-SELEX and HT-SELEX exist
  tf_sample_dir <- paste0(tf_dir,"NCAP-SELEX_cycle-4/")
  if (!dir.exists(tf_sample_dir)){
    next
  }
  
  tf_ht_dir <- paste0(tf_dir,"HT-SELEX_cycle-4/")
  if (!dir.exists(tf_ht_dir)){
    next
  }
  
  # Ignore TF if no enriched k-mers in NCAP-SELEX library
  if (!file.exists(paste0(tf_sample_dir,"final_enriched_kmers.txt"))){
    next
  }
  
  # Creating directory to save end binding reads and cyclizability predictions
  tf_cyc_dir <- paste0(tf_dir,"cyc/")
  if (!file.exists(tf_cyc_dir)) {
    dir.create(tf_cyc_dir)
  }
  
  
  sample_fastq = readFastq(get_clean_file(tf_sample_dir))
  sample_seqs <- sread(sample_fastq)
  
  ht_fastq = readFastq(get_clean_file(tf_ht_dir))
  ht_seqs <- sread(ht_fastq)
  
  motifs = read.table(paste0(tf_sample_dir,"final_enriched_kmers.txt"))$V1
  
  # Find the positions of the top enriched motifs in library sequences
  sample_positions = to_get_many_kmer_match_reads(motifs, sample_seqs)
  control_positions = to_get_many_kmer_match_reads(motifs, control_seqs)
  ht_positions = to_get_many_kmer_match_reads(motifs, ht_seqs)
  
  pos_col = ncol(sample_positions)
  
  # Filter for the sequences with motifs only at ends
  left_sample = filter_seq_indices(sample_positions,1:4)
  right_sample = filter_seq_indices(sample_positions, (pos_col-3):pos_col)
  left_control = filter_seq_indices(control_positions,1:4)
  right_control = filter_seq_indices(control_positions, (pos_col-3):pos_col)
  left_ht = filter_seq_indices(ht_positions,1:4)
  right_ht = filter_seq_indices(ht_positions, (pos_col-3):pos_col)
  
  # Only store sequences from libraries containing at least 500 end binding reads
  if (length(left_sample)+length(right_sample) < 500 || 
      length(left_control)+length(right_control) < 500 || 
      length(left_ht)+length(right_ht) < 500){
    next
  }
  
  left_sample_seqs = sample_fastq[left_sample,]
  right_sample_seqs = sample_fastq[right_sample, ]
  # Get reverse complement of reads with end binding on exits
  right_sample_rc_seqs = reverseComplement(right_sample_seqs)
  writeFasta(left_sample_seqs, paste0(tf_cyc_dir, "sample_left.fasta"))
  writeFasta(right_sample_rc_seqs, paste0(tf_cyc_dir, 
                                          "sample_right_rc.fasta"))
  # Add slurm commands for predicting cyclizability with DNAcycP
  cyc_script = append(cyc_script, 
                      glue('sbatch pred_fasta_cpu.sh {paste0(tf_cyc_dir, "sample_left.fasta")} {paste0(tf_cyc_dir, "sample_left.csv")}' ))
  cyc_script = append(cyc_script, 
                      glue('sbatch pred_fasta_cpu.sh {paste0(tf_cyc_dir, "sample_right_rc.fasta")} {paste0(tf_cyc_dir, "sample_right_rc.csv")}' ))
  
  # Calculate end binding percentages
  perc_ends = c(perc_ends, sum(sample_positions[,c(1:4,(pos_col-3):pos_col)])/sum(sample_positions))
  good_tfs = c(good_tfs, tf_name)
  
  # Do the same for nucleosome control
  left_control_seqs = control_fastq[left_control,]
  right_control_seqs = control_fastq[right_control, ]
  right_control_rc_seqs = reverseComplement(right_control_seqs)
  writeFasta(left_control_seqs, paste0(tf_cyc_dir, "control_left.fasta"))
  writeFasta(right_control_rc_seqs, paste0(tf_cyc_dir,
                                           "control_right_rc.fasta"))
  cyc_script = append(cyc_script,
                      glue('sbatch pred_fasta_cpu.sh {paste0(tf_cyc_dir, "control_left.fasta")} {paste0(tf_cyc_dir, "control_left.csv")}' ))
  cyc_script = append(cyc_script,
                      glue('sbatch pred_fasta_cpu.sh {paste0(tf_cyc_dir, "control_right_rc.fasta")} {paste0(tf_cyc_dir, "control_right_rc.csv")}' ))

  # Do the same for free DNA
  left_ht_seqs = ht_fastq[left_ht,]
  right_ht_seqs = ht_fastq[right_ht, ]
  right_ht_rc_seqs = reverseComplement(right_ht_seqs)
  writeFasta(left_ht_seqs, paste0(tf_cyc_dir, "ht_left.fasta"))
  writeFasta(right_ht_rc_seqs, paste0(tf_cyc_dir, "ht_right_rc.fasta"))
  cyc_script = append(cyc_script,
                      glue('sbatch pred_fasta_cpu.sh {paste0(tf_cyc_dir, "ht_left.fasta")} {paste0(tf_cyc_dir, "ht_left.csv")}' ))
  cyc_script = append(cyc_script,
                      glue('sbatch pred_fasta_cpu.sh {paste0(tf_cyc_dir, "ht_right_rc.fasta")} {paste0(tf_cyc_dir, "ht_right.csv/")}' ))
  
}

write.table(cyc_script, "../Cyclizability-Prediction-Website/selex_cyc_script.sh",
            quote=F,row.names = F, col.names = F)

perc_end_df = data.frame(end_bind_perc = perc_ends, tf=good_tfs)
write.csv(perc_end_df, "../data/end_binding_percent.csv")
