source("utils.R")
source("config.R")

tf_list = read.table(SELEX_tf_list_file)$V1
exprs = c("NCAP-SELEX_cycle-4", "HT-SELEX_cycle-4")
k=7

rev_comp = function(s){
  c = chartr('ATGC', 'TACG',s)
  return (intToUtf8(rev(utf8ToInt(c))))
}

for (tf_name in tf_list){
  tf_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/", tf_name, "/")
  for (expr in exprs){
    tf_sample_dir <- paste0(tf_dir, expr, "/")
    if (!dir.exists(tf_sample_dir)){
      next
    }
    
    aff = read.csv(paste0(tf_sample_dir, "aff_mat_potential_", k, ".txt"))
    if (length(aff$Kmer) < 5){
      next
    }
    motifs = aff$Kmer[1:min(nrow(aff), 10)]
    rc_motifs = sapply(motifs, rev_comp)
    final_motifs = unique(c(motifs, rc_motifs))
    write.table(final_motifs, paste0(tf_sample_dir,"final_enriched_kmers.txt"), 
                quote=F,row.names = F, col.names = F)
  }
  
}