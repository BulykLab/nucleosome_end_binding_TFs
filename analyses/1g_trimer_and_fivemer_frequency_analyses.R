#################################################################
####################   Code Initialization  ####################
############################################################

if(TRUE)
{  
  rm(list=ls())
  
  # - R Libraries upload
  if(1)
  {
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggpubr)
  library(viridis)
  }
  
  # Directories
  if(1)
  {
    Parent_dir <-  $YOUR_PARENTAL_DIR  
    #for example 
    Parent_dir <- "/Volumes/My_Book_3/PIONEAR_seq_Project_Back_Up_2022/Manuscript/preparation/Github/DNA_Flexibility_Analysis_reproduction/NCAP_SELEX_analysis/"
    setwd(Parent_dir)
    system( "mkdir MS_files/") 
    
    source("analyses/utils.R")
    source("analyses/config.R")
    ligand_data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")
    
  }
  
  # - Variables and functions
  if(1)
    {
    save_pheatmap_pdf <- function(x, filename, width=18, height=10) 
      {
      stopifnot(!missing(x))
      stopifnot(!missing(filename))
      pdf(filename, width=width, height=height)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
      }

    end_binder = read.csv(paste0(   data_dir, "end_binding.csv"))$tf
    #Exclude RFX4 because of incorrect motif enrichment
    end_binder <- setdiff(end_binder, c("RFX4"))
    non_end_binder = read.csv("data/non_end_binding.csv")$tf
    ht_thres_tf = read.csv("data/ht_thres_TF.csv")$x
    }
  }
 

#################################################################################
################    Trimer frequency analysis in NCAP-SELEX data    ##############
###################################################################################

if(1)
  {
  tf_name_vec =  read.table(SELEX_tf_list_file)$V1
  
  trimer_FW_vec <- paste0( rep(rep( DNA_BASES, each=16), time=1),
                       rep(rep( DNA_BASES, each=4), time=4), 
                       rep(rep( DNA_BASES, each=1), time=16))
  trimer_RV_vec= as.character(reverseComplement(DNAStringSet(trimer_FW_vec )))
  trimer_FW_RV_lab = which(trimer_FW_vec < trimer_RV_vec) 
  smaller_trimer_FW_vec = trimer_FW_vec [trimer_FW_RV_lab]
  smaller_trimer_RV_vec  = trimer_RV_vec [trimer_FW_RV_lab]
  trimer_FW_RV_vec = paste0(smaller_trimer_FW_vec, "/", smaller_trimer_RV_vec )

  all_tfs_trimer_freq_sample_array = array(NA, dim=c(length(trimer_FW_RV_vec),145,length(tf_name_vec )),dimnames=list(trimer_FW_RV_vec,1:145,tf_name_vec ))
  all_tfs_trimer_freq_ht_array = all_tfs_trimer_freq_sample_array 
  all_tfs_trimer_freq_control_array = all_tfs_trimer_freq_sample_array 

  all_tfs_trimer_cumulative_freq_array = array(NA, dim=c(length(trimer_FW_RV_vec),3,length(tf_name_vec )),dimnames=list(trimer_FW_RV_vec,c("NCAP","HT", "Nucleosome"),tf_name_vec ))


  for (tf_name in tf_name_vec )
    {
    print(tf_name)
    tf_dir <- paste0(ligand_data_dir, tf_name, "/")
    tf_cyc_dir <- paste0(tf_dir,"cyc/")
    if (!dir.exists(tf_cyc_dir)) {next}
    if ((!file.exists(paste0(tf_cyc_dir, "sample_left_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "sample_right_rc_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "control_left_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "control_right_rc_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "ht_left_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "ht_right_rc_primer.fasta"))))
      {next}
    sample_left = readFasta(paste0(tf_cyc_dir, "sample_left_primer.fasta"))
    sample_right_rc = readFasta(paste0(tf_cyc_dir, "sample_right_rc_primer.fasta"))
    sample = DNAStringSet(c(sread(sample_left),sread(sample_right_rc)))
    sample_dens = sapply(1:145, function(x) colMeans(trinucleotideFrequency(narrow(sample, x, x+2))))
    sample_combined_dens = sample_dens[smaller_tri_bend_FW, ] + sample_dens[smaller_tri_bend_RV, ]
    colnames(sample_combined_dens)= 1:(ncol(sample_combined_dens))
    all_tfs_trimer_freq_sample_array[,,tf_name] = sample_combined_dens
    all_tfs_trimer_cumulative_freq_array[,"NCAP" ,tf_name] = rowMeans(sample_combined_dens[,35:121])
    
    control_left = readFasta(paste0(tf_cyc_dir, "control_left_primer.fasta"))
    control_right_rc = readFasta(paste0(tf_cyc_dir, "control_right_rc_primer.fasta"))
    control = DNAStringSet(c(sread(control_left),sread(control_right_rc)))
    control_dens = sapply(1:145, function (x) colMeans(trinucleotideFrequency(narrow(control, x, x+2 ))))
    control_combined_dens = control_dens[smaller_tri_bend_FW, ] + control_dens[smaller_tri_bend_RV, ]
    colnames(control_combined_dens)= 1:(ncol(control_combined_dens))
    all_tfs_trimer_freq_control_array[,,tf_name] = control_combined_dens
    all_tfs_trimer_cumulative_freq_array[,"Nucleosome" ,tf_name] = rowMeans(control_combined_dens[,35:121])
    
    ht_left = readFasta(paste0(tf_cyc_dir, "ht_left_primer.fasta"))
    ht_right_rc = readFasta(paste0(tf_cyc_dir, "ht_right_rc_primer.fasta"))
    ht = DNAStringSet(c(sread(ht_left),sread(ht_right_rc)))
    ht_dens = sapply(1:145, function (x) colMeans(trinucleotideFrequency(narrow(ht, x, x+2 )))) 
    ht_combined_dens = ht_dens[smaller_tri_bend_FW, ] + ht_dens[smaller_tri_bend_RV, ]
    colnames(ht_combined_dens)= 1:(ncol(ht_combined_dens))
    all_tfs_trimer_freq_ht_array[,,tf_name] = ht_combined_dens
    all_tfs_trimer_cumulative_freq_array[,"HT" ,tf_name] = rowMeans(ht_combined_dens[,35:121])
  
}

save(all_tfs_trimer_cumulative_freq_array, file=paste0(data_dir, "all_tfs_trimer_cumulative_freq_array.rdata"))
save(all_tfs_trimer_freq_sample_array, file=paste0(data_dir,"all_tfs_trimer_freq_sample_array.rdata") )
save(all_tfs_trimer_freq_control_array, file=paste0(data_dir,"all_tfs_trimer_freq_control_array.rdata") )
save(all_tfs_trimer_freq_ht_array, file=paste0(data_dir,"all_tfs_trimer_freq_ht_array.rdata") )

}
############################################################
################       Plots for Fig S10     ################
################################################################

if(1)
  {
  load( file=paste0(data_dir, "all_tfs_trimer_cumulative_freq_array.rdata"))
 
  tf_name_vec  <- dimnames(all_tfs_trimer_cumulative_freq_array)[[3]]
  tf_name_vec  <- intersect(tf_name_vec , ht_thres_tf)
  tf_end_binders <- intersect(  tf_name_vec, end_binder)
  tf_non_end_binders <- intersect(  tf_name_vec, non_end_binder)

  all_tfs_trimer_cumulative_freq_end_binders_array <- all_tfs_trimer_cumulative_freq_array[,, tf_end_binders ]
  all_tfs_trimer_cumulative_freq_non_end_binders_array <- all_tfs_trimer_cumulative_freq_array[,, tf_non_end_binders ]

  all_TFs_trimer_NCAP_end_binders <- t(all_tfs_trimer_cumulative_freq_end_binders_array[,"NCAP",])
  all_TFs_trimer_NCAP_non_end_binders <- t(all_tfs_trimer_cumulative_freq_non_end_binders_array[,"NCAP",])
  all_TFs_trimer_NCAP_both_binders <- rbind(all_TFs_trimer_NCAP_end_binders,all_TFs_trimer_NCAP_non_end_binders)

  NCAP_both_binders_cumulative_freq=pheatmap(all_TFs_trimer_NCAP_both_binders  ,  cluster_cols=F, cluster_rows = F,#scale="column", 
                                              color=colorRampPalette(c("darkblue", "white","red"))(100),
                                              breaks = seq(0,0.06, length.out=101),gaps_row = 69,
                                              show_colnames = T)
  figure_address <- paste0( "figures/NCAP_both_end_binders_cumulative_freq.pdf")
  save_pheatmap_pdf(NCAP_both_binders_cumulative_freq, figure_address)
  dev.off()
  system( paste( "cp", figure_address, "MS_files/Fig_S10a_middle.pdf"))
  write.table( cbind(TF=rownames(all_TFs_trimer_NCAP_both_binders), all_TFs_trimer_NCAP_both_binders), file="MS_files/Fig_S10a_middle_datasets.txt" ,  row.names=F, col.names = T, quote=F,       sep = "\t")
  
  all_TFs_trimer_Nucleosome_end_binders <- t(all_tfs_trimer_cumulative_freq_end_binders_array[,"Nucleosome",])
  all_TFs_trimer_Nucleosome_non_end_binders <- t( all_tfs_trimer_cumulative_freq_non_end_binders_array[,"Nucleosome",] ) 
  all_TFs_trimer_Nucleosome_both_binders <- rbind(all_TFs_trimer_Nucleosome_end_binders, all_TFs_trimer_Nucleosome_non_end_binders)

  Nucleosome_both_binders_cumulative_freq=pheatmap(all_TFs_trimer_Nucleosome_both_binders ,
                                                 cluster_cols=F, cluster_rows = F,#scale="column", 
                                                    color=colorRampPalette(c("darkblue", "white","red"))(100),
                                                    breaks = seq(0,0.06, length.out=101),gaps_row = 69,
                                                    show_colnames = T)
  figure_address <- "figures/Nucleosome_both_binders_cumulative_freq.pdf"
  save_pheatmap_pdf(Nucleosome_both_binders_cumulative_freq, figure_address)
  dev.off()

  system( paste( "cp", figure_address, "MS_files/Fig_S10a_top.pdf"))
  write.table( cbind(TF=rownames(all_TFs_trimer_Nucleosome_both_binders), all_TFs_trimer_Nucleosome_both_binders), file="MS_files/Fig_S10a_top_datasets.txt" ,  row.names=F, col.names = T, quote=F,       sep = "\t")
  
  all_TFs_trimer_HT_end_binders <- t(all_tfs_trimer_cumulative_freq_end_binders_array[,"HT",])
  all_TFs_trimer_HT_non_end_binders <- t(all_tfs_trimer_cumulative_freq_non_end_binders_array[,"HT",])
  all_TFs_trimer_HT_both_binders <- rbind(all_TFs_trimer_HT_end_binders, all_TFs_trimer_HT_non_end_binders) 

  HT_both_binders_cumulative_freq=pheatmap(all_TFs_trimer_HT_both_binders ,           cluster_cols=F, cluster_rows = F,#scale="column", 
                                                     color=colorRampPalette(c("darkblue", "white","red"))(100),
                                                    breaks = seq(0,0.06, length.out=101),gaps_row = 69,
                                                    show_colnames = T)
  figure_address <- "figures/HT_both_binders_cumulative_freq.pdf"
  save_pheatmap_pdf(HT_both_binders_cumulative_freq, figure_address)
  dev.off()
  
  system( paste( "cp", figure_address, "MS_files/Fig_S10a_bottom.pdf"))
  write.table( cbind(TF=rownames(all_TFs_trimer_HT_both_binders), 
                     all_TFs_trimer_HT_both_binders), file="MS_files/Fig_S10a_bottom_datasets.txt" ,  row.names=F, col.names = T, quote=F,       sep = "\t")
  


  AAA_both_binders <- cbind( Nucleosome=all_TFs_trimer_Nucleosome_both_binders[,"AAA/TTT"],NCAP=all_TFs_trimer_NCAP_both_binders[,"AAA/TTT"],HT=all_TFs_trimer_HT_both_binders[,"AAA/TTT"])
  AAA_both_binders_freq =pheatmap(AAA_both_binders  ,           cluster_cols=F, cluster_rows = F,#scale="column", 
                                         color=colorRampPalette(c("darkblue", "white","red"))(100),
                                         breaks = seq(0,0.06, length.out=101),gaps_row = 69,
                                         show_colnames = T)
 
  figure_address <- "figures/AAA_both_binders_cumulative_freq.pdf"
  save_pheatmap_pdf(AAA_both_binders_freq, figure_address)
  dev.off()
  
  system( paste( "cp", figure_address, "MS_files/Fig_S10b.pdf"))
  write.table( cbind(TF=rownames(AAA_both_binders), 
                     AAA_both_binders), file="MS_files/Fig_S10b_datasets.txt" ,  row.names=F, col.names = T, quote=F,       sep = "\t")
  
}

########################################################
##################### 5-mer frequency analysis  ###########
#########################################################

if(1)
  {

  tf_name_vec =  read.table(SELEX_tf_list_file)$V1
  all_tfs_fivemer_freq_sample_matrix = array(0, dim=c( 143,length(  tf_name_vec)),dimnames=list( 1:143,  tf_name_vec))
  print(tf_name_vec)
  for (tf_name in tf_name_vec)
    {
    print(tf_name)
    tf_dir <- paste0(ligand_data_dir, tf_name, "/")
    tf_cyc_dir <- paste0(tf_dir,"cyc/")
    if (!dir.exists(tf_cyc_dir)) 
      {    next}
    if ((!file.exists(paste0(tf_cyc_dir, "sample_left_primer.fasta"))) ||
      (!file.exists(paste0(tf_cyc_dir, "sample_right_rc_primer.fasta"))) )
      {next  }
  
    sample_left = readFasta(paste0(tf_cyc_dir, "sample_left_primer.fasta"))
    sample_right_rc = readFasta(paste0(tf_cyc_dir, "sample_right_rc_primer.fasta"))
    sample = DNAStringSet( c(sread(sample_left),sread(sample_right_rc)))
 
    fivemers_to_scan <- c("TAAAA", "TTAAA","TTTAA", "TTTTA") 
    fivemer_start_vec <- c()
    for(fivemer in fivemers_to_scan)
      {
      hits <- unlist(vmatchPattern(fivemer ,      sample , fixed = FALSE))@start
      fivemer_start_vec <- c(fivemer_start_vec ,hits)
      }
    sample_fivemer_dens_vec <- table( fivemer_start_vec)/length(sample)
    all_tfs_fivemer_freq_sample_matrix[names( sample_fivemer_dens_vec),tf_name] = sample_fivemer_dens_vec
    }
  save(  all_tfs_fivemer_freq_sample_matrix,file=paste0( data_dir,"all_tfs_fivemer_freq_sample_matrix.rdata") )
  }

############################################################
################       Plots for Fig 5c-e     ################
################################################################

if(1){
  
  #load all_tfs_fivemer_freq_sample_matrix
  load(  file=paste0( data_dir,"all_tfs_fivemer_freq_sample_matrix.rdata") )
  
  tf_name_vec <- dimnames( all_tfs_fivemer_freq_sample_matrix)[[2]]
  tf_name_vec <- intersect(tf_name_vec, ht_thres_tf)
  tf_end_binders <- intersect(tf_name_vec, end_binder)
  tf_non_end_binders <- intersect(tf_name_vec, non_end_binder)
  
  tf_end_binders_fivemer_freq_sample_matrix<- all_tfs_fivemer_freq_sample_matrix[,tf_end_binders ]
  tf_non_end_binders_fivemer_freq_sample_matrix <- all_tfs_fivemer_freq_sample_matrix[,tf_non_end_binders ]
  
  both_binder_freq_matrix <- t(cbind(  tf_end_binders_fivemer_freq_sample_matrix , 
                                        tf_non_end_binders_fivemer_freq_sample_matrix ))
    
  NCAP_both_binders_kmer=pheatmap(both_binder_freq_matrix,cluster_cols=F, cluster_rows = F, #scale="row", 
                                    color=colorRampPalette(c("darkblue", "white","red"))(100),
                                    breaks = seq(0,0.012, length.out=101), gaps_row=69,
                                    show_colnames = )
  figure_address <- paste0( "figures/NCAP_both_binders_", kmer ,"_full_reads.pdf")
  save_pheatmap_pdf(NCAP_both_binders_kmer,    figure_address,12,12)
  dev.off()
  system( paste( "cp", figure_address, "MS_files/Fig_5c.pdf"))
  write.table( cbind(TF=rownames(both_binder_freq_matrix), both_binder_freq_matrix), file="MS_files/Fig_5c_datasets.txt" ,  row.names=F, col.names = T, quote=F,       sep = "\t")
    
  end_binder_freq_submatrix <-  end_binder_freq_matrix[,36:119]
  non_end_binder_freq_submatrix <-  non_end_binder_freq_matrix[,36:119]
    
  right_ACF_end_binders <- rep(NA,length(TF_list_end_binders))
  names( right_ACF_end_binders ) <- TF_list_end_binders
  right_ACF_non_end_binders <-  rep(NA,length(TF_list_non_end_binders))
  names(  right_ACF_non_end_binders ) <- TF_list_non_end_binders
  start_bp <- 1
  mid_bp <- 40
  stop_bp <- 80
  for(tf_name in rownames(end_binder_freq_submatrix))
    {
    right_ACF_end_binders[tf_name]= acf(end_binder_freq_submatrix[tf_name,    mid_bp : stop_bp], pl=F)$acf[11]
    }
    
  for(tf_name in rownames(non_end_binder_freq_submatrix))
    {
     right_ACF_non_end_binders[tf_name]= acf(non_end_binder_freq_submatrix[tf_name,    mid_bp : stop_bp], pl=F)$acf[11]
     }
    
  left_freq_end_binders <-  rowMeans( end_binder_freq_submatrix[,start_bp:    mid_bp ]) 
  left_freq_non_end_binders <-  rowMeans( non_end_binder_freq_submatrix[,start_bp:mid_bp]) 
 
  x1 <- left_freq_end_binders
  x2 <- left_freq_non_end_binders
  x_def <- "left_freq"
  y1 <- right_ACF_end_binders
  y2 <- right_ACF_non_end_binders
  y_def <- "right_ACF"
  plot( NA,xlim=c(min( x1 ,x2), max( x1, x2)),
        ylim=c(min( y1 , y2), max( y1 , y2)),
        xlab=paste(kmer, x_def), ylab=paste(kmer,  y_def))
  points( x=x1, y= y1 ,pch=16, col="black")
  points( x=x2, y=y2 ,pch=16, col="red")
  figure_address <- paste0( "figures/plot_", kmer,'_', x_def , '_VS_' , y_def , '.png')
  dev.copy(png, figure_address,500, 500, res=144)
  dev.off()
  system( paste( "cp", figure_address, "MS_files/Fig_5d.pdf"))
  matrix_to_save <- cbind(TF=names(c(x1,x2)),Left_freq=c(x1,x2),Right_ACF=c(y1,y2), Type=c( rep("End_Binder", length(x1)), rep("Non_End_Binder", length(x2)) ))
  write.table(   matrix_to_save , file="MS_files/Fig_5d_datasets.txt" ,  row.names=F, col.names = T, quote=F,       sep = "\t")
    
  z1_df <- data.frame(names(x1),x1,y1,"End binders")
  z2_df <- data.frame(names(x2),x2,y2,"Non-end binders")
  names(z1_df) <-names(z2_df) <- c('tf','ACF','Freq','end_binding' ) 
  z_df <-rbind(z1_df,z2_df)
  z_df$rank <-  z_df$end_binding=="End binders"
  mylogit <- glm(  rank ~ ACF +  Freq , data = z_df, family = "binomial")
  end_pred <-  mylogit$fitted.values[TF_list_end_binders]
  non_end_pred <-  mylogit$fitted.values[TF_list_non_end_binders]
  z_df$logit <- c(end_pred, non_end_pred)
  wilcox.test( end_pred, non_end_pred )
  ggplot(z_df, aes(x=end_binding,  y=logit))+geom_boxplot(color="#7F3F98")+ geom_jitter(color="#7F3F98",size=0.5)+
        labs(  x="",   y="Odds" )+ylim(c(0,1))+theme_classic()+
        theme(title=element_text(size=7), text = element_text(size=7))
  figure_address <-paste0("figures/Boxplot_logit_",kmer,".png")
  ggsave(figure_address , width=50, height=50, units = "mm")
  dev.off()
  system( paste( "cp", figure_address, "MS_files/Fig_5d.pdf"))
  write.table( z_df, file="MS_files/Fig_5c_d_datasets.txt" ,  row.names=F, col.names = T, quote=F,       sep = "\t")
    
  }

#####################   End  ##################  
