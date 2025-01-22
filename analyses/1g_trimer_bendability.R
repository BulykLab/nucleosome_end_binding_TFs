source("utils.R")
source("config.R")

library(pheatmap)

trimer_plot_dir = paste0(data_dir, "tri_mer_plots/")
if (!dir.exists(trimer_plot_dir)){
  dir.create(trimer_plot_dir)
}

# Code adopted from stackoverflow
save_pheatmap_pdf <- function(x, filename, width=18, height=10) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

tf_list = c("SOX11", "CEBPB", "FOXA2")#read.table(SELEX_tf_list_file)$V1
data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")

tri_bend = read.table(tri_mer_bend_file)
tri_bend$RV = as.character(reverseComplement(DNAStringSet(tri_bend$Seq)))

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


  sample_dens = sapply(1:145, function(x) colMeans(trinucleotideFrequency(narrow(sample, x, x+2))))
  smaller_tri = tri_bend[which(tri_bend$Seq < tri_bend$RV), ]
  smaller_tri_sorted = smaller_tri[order(smaller_tri$BendDNasel), ]
  sample_combined_dens = sample_dens[smaller_tri_sorted$Seq, ] +
    sample_dens[smaller_tri_sorted$RV, ]

  smaller_tri_sorted$sample_freq = rowMeans(sample_combined_dens[,35:121])

  freq = spectrum(sample_combined_dens[1,35:121])$freq
  sample_periodicity = sapply(1:nrow(sample_combined_dens), function(i) return (sum(spectrum(sample_combined_dens[i,35:121])$spec[which((freq <= 0.11) | (freq >= 0.091 ))])))
  smaller_tri_sorted$FreqPeriodicityIntensity = sample_periodicity

  colnames(sample_combined_dens)= 1:(ncol(sample_combined_dens))
  sample_hm=pheatmap(sample_combined_dens, 
  cluster_cols=F, cluster_rows = F,#scale="column", 
  labels_row = paste0(smaller_tri_sorted$Seq, "/", smaller_tri_sorted$RV),
  annotation_row = smaller_tri_sorted[,c("BendDNasel", "BendConsensus", "FreqPeriodicityIntensity")],
  color=colorRampPalette(c("darkblue", "white","red"))(100),
  breaks = seq(0,0.06, length.out=101),
  show_colnames = T)


  control_dens = sapply(1:145, function (x) colMeans(trinucleotideFrequency(narrow(control, x, x+2 ))))
  control_combined_dens = control_dens[smaller_tri_sorted$Seq, ] + control_dens[smaller_tri_sorted$RV, ]
  smaller_tri_sorted$control_freq = rowMeans(control_combined_dens[,35:121])


  control_periodicity = sapply(1:nrow(control_combined_dens), function(i) return (sum(spectrum(control_combined_dens[i,35:121])$spec[which((freq <= 0.11) | (freq >= 0.091 ))])))
  smaller_tri_sorted$FreqPeriodicityIntensity = control_periodicity

  colnames(control_combined_dens)= 1:(ncol(control_combined_dens))
  control_hm=pheatmap(control_combined_dens, 
  cluster_cols=F, cluster_rows = F,#scale="column", 
  labels_row = paste0(smaller_tri_sorted$Seq, "/", smaller_tri_sorted$RV),
  annotation_row = smaller_tri_sorted[,c("BendDNasel", "BendConsensus", "FreqPeriodicityIntensity")],
  color=colorRampPalette(c("darkblue", "white","red"))(100),
  breaks = seq(0,0.06, length.out=101),
  show_colnames = T)


  ht_dens = sapply(1:145, function (x) colMeans(trinucleotideFrequency(narrow(ht, x, x+2 ))))
  ht_combined_dens = ht_dens[smaller_tri_sorted$Seq, ] + ht_dens[smaller_tri_sorted$RV, ]
  smaller_tri_sorted$ht_freq = rowMeans(ht_combined_dens[,35:121])


  ht_periodicity = sapply(1:nrow(ht_combined_dens), function(i) return (sum(spectrum(ht_combined_dens[i,35:121])$spec[which((freq <= 0.11) | (freq >= 0.091 ))])))
  smaller_tri_sorted$FreqPeriodicityIntensity = ht_periodicity


  colnames(ht_combined_dens)= 1:(ncol(ht_combined_dens))
  ht_hm=pheatmap(ht_combined_dens, 
  cluster_cols=F, cluster_rows = F,#scale="column", 
  labels_row = paste0(smaller_tri_sorted$Seq, "/", smaller_tri_sorted$RV),
  annotation_row = smaller_tri_sorted[,c("BendDNasel", "BendConsensus", "FreqPeriodicityIntensity")],
  color=colorRampPalette(c("darkblue", "white","red"))(100),
  breaks = seq(0,0.06, length.out=101),
  show_colnames = T)


  save_pheatmap_pdf(sample_hm, paste0(trimer_plot_dir, tf_name, "_sample_combined.pdf"))
  save_pheatmap_pdf(control_hm, paste0(trimer_plot_dir, tf_name, "_control_combined.pdf"))
  save_pheatmap_pdf(ht_hm,paste0(trimer_plot_dir, tf_name,  "_ht_combined.pdf"))


  all_tri_freq = smaller_tri_sorted[, c("sample_freq", "ht_freq", "control_freq")]

  get_num_AT = function(seq) sum(charToRaw(seq) == charToRaw("A"))+sum(charToRaw(seq) == charToRaw("T"))
  smaller_tri_sorted$ATnum = sapply(smaller_tri_sorted$Seq, get_num_AT)
  sorted_all_tri_freq = all_tri_freq[order(smaller_tri_sorted$ATnum, decreasing = T),]

  all_tri = pheatmap(all_tri_freq, 
  cluster_cols=F, cluster_rows = F,#scale="column", 
  labels_row = paste0(smaller_tri_sorted$Seq, "/", smaller_tri_sorted$RV),
  annotation_row = smaller_tri_sorted[,c("BendDNasel", "BendConsensus")],
  color=colorRampPalette(c("darkblue", "white","red"))(100),
  breaks = seq(0,0.06, length.out=101),
  show_colnames = T)


  save_pheatmap_pdf(all_tri, paste0(trimer_plot_dir, tf_name, "_combined.pdf"), 6, 10)

  sorted_all_tri = pheatmap(sorted_all_tri_freq, 
  cluster_cols=F, cluster_rows = F,#scale="column", 
  labels_row = paste0(smaller_tri_sorted[order(smaller_tri_sorted$ATnum, decreasing = T),]$Seq, "/", smaller_tri_sorted[order(smaller_tri_sorted$ATnum, decreasing = T),]$RV),
  color=colorRampPalette(c("darkblue", "white","red"))(100),
  breaks = seq(0,0.06, length.out=101),
  show_colnames = T)

  save_pheatmap_pdf(sorted_all_tri, paste0(trimer_plot_dir, tf_name, "_combined_AT.pdf"), 6,10)

}


