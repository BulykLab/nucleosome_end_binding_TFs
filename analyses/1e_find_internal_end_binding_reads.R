source("utils.R")
source("config.R")

tf_list = read.table(SELEX_tf_list_file)$V1
data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")
k=7

cyc_script = c() # Store cyc commands for running DNA cycP on slurm


for (tf_name in tf_list){
  tf_dir <- paste0(data_dir, tf_name, "/")
  
  # Ignore TF if not both NCAP-SELEX and HT-SELEX exist
  tf_sample_dir <- paste0(tf_dir,"NCAP-SELEX_cycle-4/")
  if (!dir.exists(tf_sample_dir)){
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
  
  motifs = read.table(paste0(tf_sample_dir,"final_enriched_kmers.txt"))$V1
  sample_positions = to_get_many_kmer_match_reads(motifs, sample_seqs)
  pos_col = ncol(sample_positions)
  
  for (internal_pos in seq(1, 20, 5)){
    left_pos = (internal_pos):(internal_pos+4)
    right_pos = (pos_col-internal_pos-3):(pos_col-internal_pos+1)
    left_sample = filter_seq_indices(sample_positions,left_pos)
    right_sample = filter_seq_indices(sample_positions, right_pos)
    if (100<length(left_sample)+length(right_sample)){
      left_sample_seqs = sample_fastq[left_sample,]
      left_sample_id = left_sample_seqs@id
      left_sample_seqs = DNAStringSet(paste0(FW_primer, 
                                             as.character(sread(left_sample_seqs)),
                                             RV_primer))
      right_sample_seqs = sample_fastq[right_sample, ]
      right_sample_id = right_sample_seqs@id
      right_sample_seqs = DNAStringSet(paste0(FW_primer, 
                                              as.character(sread(right_sample_seqs)), 
                                              RV_primer))
      right_sample_rc_seqs = reverseComplement(right_sample_seqs)
      
      left_sample_fa_file = paste0(tf_cyc_dir, "sample_left_",min(left_pos),
                                   "_", max(left_pos) ,"_primer.fasta")
      left_sample_cyc_file = paste0(tf_cyc_dir, "sample_left_",min(left_pos),
                                   "_", max(left_pos) ,"_primer.csv")
      writeFasta(ShortRead(left_sample_seqs, left_sample_id), left_sample_fa_file)
      
      right_sample_rc_fa_file = paste0(tf_cyc_dir, "sample_right_rc_",min(right_pos),
                                       "_", max(right_pos) ,"_primer.fasta")
      right_sample_rc_cyc_file = paste0(tf_cyc_dir, "sample_right_rc_",min(right_pos),
                                       "_", max(right_pos) ,"_primer.csv")
      writeFasta(ShortRead(right_sample_rc_seqs, right_sample_id), right_sample_rc_fa_file)
      
      cyc_script = append(cyc_script, 
                          glue('sbatch pred_fasta_cpu.sh {left_sample_fa_file} {left_sample_cyc_file}' ))
      cyc_script = append(cyc_script, 
                          glue('sbatch pred_fasta_cpu.sh {right_sample_rc_fa_file} {right_sample_rc_cyc_file}' ))
    }
  }
  
}

write.table(cyc_script, "../Cyclizability-Prediction-Website/selex_internal_cyc_script.sh",
            quote=F,row.names = F, col.names = F)