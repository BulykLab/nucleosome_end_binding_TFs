library(ShortRead)

# Get the clean sequence file from a sample directory
get_clean_file <- function(ddir) {
  all_files <- list.files(ddir)
  clean_files <- all_files[grepl("clean.fastq", all_files)]
  return(paste0(ddir, clean_files[1]))
}


# Find k-mer match positions from sequences
to_get_many_kmer_match_reads <- function(kmers, many_sequences,seq_len=101) { # } ,sequence_start_vec ){
  merged_single_sequence_from_all <- paste0(many_sequences, collapse = " ")
  pos_mat <- matrix(0, length(many_sequences), (seq_len+1)-nchar(kmers[1]))
  for (single_kmer in kmers){
    matches_kmer <- matchPattern(single_kmer, 
                                 merged_single_sequence_from_all)
    start_sites <- start(matches_kmer) %% (seq_len+1)
    read_index = start(matches_kmer) %/% (seq_len+1)
    for (i in 1:length(read_index)){
      pos_mat[read_index[i]+1, start_sites[i]] = 1
    }
  }
  return(pos_mat)
}


filter_seq_indices = function(seq_pos, col_range){
  return (which((rowSums(seq_pos[,col_range])>0)&(rowSums(seq_pos[,-col_range])==0)))
}


load_norm_cyc_from_file = function(filename){
  return (read.csv(filename)$c_score_norm)
}

load_cyc_from_dir = function(ddir){
  cyc_files = list.files(ddir, full.names = T)
  cyc_dfs = lapply(cyc_files, load_norm_cyc_from_file)
  return (do.call(cbind, cyc_dfs))
}

build_plot_df = function(df, tf_name, label){
  avg = rowMeans(df)
  n=ncol(df)
  plot_df = data.frame( tf=tf_name, n=n,
                        dist=1:(nrow(df)), avg=avg,
                        low=avg-qt(p=0.025, df=n-1,lower.tail=F)*colSds(t(df))/sqrt(n),
                        high=avg+qt(p=0.025, df=n-1,lower.tail=F)*colSds(t(df))/sqrt(n),
                        group=paste0(label))
  #group=paste0(label, " (n=", n, ")"))
  return(plot_df)
}

get_slopes = function(cyc, inpath, outpath=NULL){
  storage <- c()
  cyc_df = as.data.frame(cyc)
  colnames(cyc_df) = paste0("r", 1:ncol(cyc))
  cyc_df$dist = 1:nrow(cyc)
  cyc_df = cyc_df %>% filter((dist >=30)& (dist <= 62))
  for(i in colnames(cyc_df)[-ncol(cyc_df)]){
    storage <- append(storage, lm(get(i) ~ dist, cyc_df)$coefficients[["dist"]])
  }
  slope_ord = order(storage, decreasing=T)
  slope = list.files(inpath)[slope_ord]
  if (!is.null(outpath)){
    write.csv(cbind(slope_ord,slope), outpath)
  }
  return (data.frame(slope_ord=slope_ord,slope=slope, slope_val=storage[slope_ord]))
  
}
