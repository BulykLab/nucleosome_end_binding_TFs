####################################################################################
################################ To initialize the program #########################
####################################################################################
#
rm(list=ls())
#
if(1){
  
  
  # In this section, several features of the program are initialized, in a short a global manner, for example:
  
  # - Library Installation in Remote (optional)
  if(0){
    
    R_packages_dir <-   YOUR_R_PACKAGE_PATH
    
    install.packages("rJava", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages('seqinr', repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages('Cairo', repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages("reshape2", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages("ggplot2", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages("gplots", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages("readxl", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages( "abind", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages( "bedr", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages("stringr", repos="http://cran.r-project.org", lib=R_packages_dir)
    install.packages("qdap", repos="http://cran.r-project.org", lib=R_packages_dir)
    
    if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager", repos="http://cran.r-project.org", lib=R_packages_dir)
    BiocManager::install(version = "3.11")
    BiocManager::install("SELEX")
    BiocManager::install("seqLogo")
    BiocManager::install("ShortRead")
    
  }
  
  # - Libraries upload
  if(0){
    library(tidyr)
    library("rJava")
    options(java.parameters="-Xmx2048m")     # to try to avoid  "java.lang.OutOfMemoryError: Java heap space"
    
    # If not by default present, you may need 
    # 1 - to install the JDK, which you can download as .dmg file from Java SE Downloads page. 
    # https://www.oracle.com/java/technologies/javase-downloads.html
    # 2 -  and then to config java environment on your MAC/workstation by
    # opening the .bash_profile in your home directory
    # vi ~/.bash_profile
    # and add those 6 lines in it:
    #    JAVA_HOME=/Library/Java/JavaVirtualMachines/jdkXXXXXX.jdk/Contents/Home
    #   PATH=$JAVA_HOME/bin:$PATH:.
    # CLASSPATH=$JAVA_HOME/lib/tools.jar:$JAVA_HOME/lib/dt.jar:.
    # export JAVA_HOME
    # export PATH
    # export CLASSPATH
    #### Where in Line 1 you must write jdk version of your java
    # Then run
    #  source ~/.bash_profile
    # to use these config immediately then you can check by
    # java -version
    
    library(parallel)
    library(BiocGenerics)
    library(Biostrings)
    library(SELEX)
    library(grid)
    library(seqLogo)
    library(seqinr)

    suppressMessages(library(ShortRead))
 
    library(reshape2)
    library(ggplot2)
    library(gplots) 
    library(readxl)
    library(abind)
    library(bedr)
    library(stringr)
    library(plyr)
    library(qdap)
    library(pheatmap)
    library(fossil)
    library(RColorBrewer)
    library(stringi)
    library(dplyr)
    library(tm)
    library(ggnewscale)
  }
  
  # - R packages to upload
  if(1){

  library(ShortRead)
  library(Biostrings)
    library(Cairo)
    library(gplots) 
    library(tidyr)
    library(dplyr)
    library(ggplot2)
  }
  # - Global pathways upload
  if(1){
    
    Parent_dir <- 'YOUR_PARENT'
    # for example:  
    Parent_dir <-'/Volumes/My_Book_3/PIONEAR_seq_Project_Back_Up_2022/Manuscript/preparation/GEO_repository/MNase_seq_datasets_to_upload/' 
    
    if (! dir.exists( Parent_dir ))
    {system( paste( "mkdir ",    Parent_dir ) )}
    setwd(  Parent_dir )
    # Directory with the sequencing read files paired by PEAR 
    Paired_dir <-  paste( sep="", Parent_dir, "/MNase_seq_Paired_files/") 
    list.files(Paired_dir)
    # it should contain the following files :
    # "F3-203-1-MNase.assembled.fastq.gz"  NRCAM single template  Replicate 1 
    # "G3-203-2-Mnase.assembled.fastq.gz" NRCAM single template Replicate 2
    # "E2-337-2-Mnase.assembled.fastq.gz" NRCAM single template with TTGCGCAA insertion on left side
    # "F2-338-1-Mnase.assembled.fastq.gz" NRCAM single template with TTGCGCAA insertion on right side
    # "C1-WANC-1-Input.assembled.fastq.gz" Pionear-seq NCP input library before MNase Replicate 1 
    # "D1-WANC-1-Mnase.assembled.fastq.gz" Pionear-seq NCP input library after MNase Replicate 1 
    # "G1-WANC-2-Input.assembled.fastq.gz" Pionear-seq NCP input library before MNase Replicate 2
    # "H1-WANC-2-Mnase.assembled.fastq.gz" Pionear-seq NCP input library after MNase Replicate 2
     
    Graphs_dir <-paste( sep="", Parent_dir, "Graphs/" )
    if (! dir.exists(  Graphs_dir )) { system( paste( "mkdir ",   Graphs_dir) ) }
    
    Sequencing_Statistics_dir <- paste(sep="", Parent_dir,"Sequencing_Statistics/")
    if(!file.exists(  Sequencing_Statistics_dir))  {   system( paste( "mkdir ",    Sequencing_Statistics_dir ) )}
     
  }  
  
  # - Global variables 
  if(1){ 
     
    template_name_vec <- c(  "ALB1", "CX3" , "NRCAM" , "W601" )
    template_parallel_vec<-c(1:length(template_name_vec))
    
    ALB1_sequence= toupper("CTGCTCTGTCAGCAGGGCACTGTACTTGCTGATACCAGGGAATGTTTGTTCTTAAATACCATCATTCCGGACGTGTTTGCCTTGGCCAGTTTTCCATGTACATGCAGAAAGAAGTTTGGACTGATCCATACAGTCCTCTGCCTTTTC")
    CX3_sequence=toupper("GGGCCTCTCGGCTGCTGATCTTCAGCTGGTTGCTGAGAGTTGCAGCATTGCTGAGTCTTAGCAATGGATACTTgCCGATTCCCCTCACAAAAATAGGTCAGTCTGTCTGGCTAGTTCTGTACTTGCAGACACAGGGCATGTGGGGTT")
    NRCAM_sequence=toupper("ACTTCTGAAACAGATGACTCCCAGCAGCTGCTGCCTGTGGCgCACAGGGCTaCCTGCCCTGCATGACAGCTGCACATCACATCCTGTGGTCATACTACTTCAGCCGCTTCTACGGCCAGATACAAAAGTGGGTGGcGAACATAGGCA")
    W601_sequence <- "ATCGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACcTGTCAGATATATAGATCCGAT"
    sequence_template_vec <- c( ALB1_sequence , CX3_sequence , NRCAM_sequence , W601_sequence )
    
    
    sequence_337 <- sequence_338 <-sequence_203 <- DNAString( NRCAM_sequence )
    sequence_338[114:121] <- sequence_337[32:39] <- DNAString( "TTGCGCAA")
  
    sequence_vec <- c( "203"=sequence_203 , "337"=sequence_337 , "338"=sequence_338 )
    
    tiled_ALB1_start_vec <-  seq(21,102,3 )
    tiled_CX3_start_vec <- seq(18,108,3 )
    tiled_NRCAM_start_vec <-  seq(26,107,3 )
    tiled_W601_start_vec <- seq(22,100,3 )
    tile_sequence_start_vec_list <- list(tiled_ALB1_start_vec, tiled_CX3_start_vec, tiled_NRCAM_start_vec, tiled_W601_start_vec)
    names( sequence_template_vec) <- names( tile_sequence_start_vec_list) <-    template_name_vec
  
    mapping <- diag(4)
    dimnames(mapping) <- list(DNA_BASES, DNA_BASES)
    
  }
}
#
####################################################################################
################################ Analysis of single templates #########################
#####################################################################################
#
if(1)
{
  paired_file_vec <-   list.files( Paired_dir  )  
  
  single_sequence_filename_vec <- paired_file_vec[grepl( '-', paired_file_vec)&!grepl( 'WANC', paired_file_vec)]
  single_sequence_filename_vec <- single_sequence_filename_vec [file.size(paste( sep="", Paired_dir, single_sequence_filename_vec ) )>10^5]

  for( single_sequence_filename in single_sequence_filename_vec )
  {
    paired_file <- paste( sep="", Paired_dir,single_sequence_filename )
    dataframe_file_name <- sub( ".assembled.fastq.gz","", substring(paired_file, regexpr("Paired_files/", paired_file) + 13)) 
    sequence_number <- strsplit( single_sequence_filename, "-")[[1]][2]
    sequence_info <- sequence_vec[sequence_number]

    if(0) #computationally expensive: this part aligns the digested reads to the original template
    {
      paired_file_fastq <- readFastq(  paired_file )
      paired_file_fastq <- trimTailw( paired_file_fastq, 1, "?", halfwidth=147) 
      paired_file_fastq <- paired_file_fastq[width( paired_file_fastq@sread )<148  ]
      paired_file_fastq <- paired_file_fastq[width( paired_file_fastq@sread )>80  ]

      pairwiseAlignment_vec <-  FW_pairwiseAlignment_vec <- pairwiseAlignment (  pattern=paired_file_fastq@sread  ,subject= sequence_info, type='overlap', substitutionMatrix=mapping)
      RC_pairwiseAlignment_vec <- pairwiseAlignment (  pattern=reverseComplement(paired_file_fastq)@sread , subject=sequence_info , type='overlap', substitutionMatrix=mapping)
      RC_lab <-  which(FW_pairwiseAlignment_vec@score<  RC_pairwiseAlignment_vec@score)
      pairwiseAlignment_vec  [ RC_lab ]  <-  RC_pairwiseAlignment_vec[ RC_lab ]   
      pairwiseAlignment_vec <-  pairwiseAlignment_vec[ pairwiseAlignment_vec@pattern@range@width>80]
      pairwiseAlignment_vec <-  pairwiseAlignment_vec[ sapply( pairwiseAlignment_vec@pattern@mismatch, length)<6]
      pairwiseAlignment_vec <-  pairwiseAlignment_vec[  pairwiseAlignment_vec@subject@range@width==pairwiseAlignment_vec@pattern@range@width]
      save( pairwiseAlignment_vec, file=paste(sep="", Sequencing_Statistics_dir, dataframe_file_name, "_pairwiseAlignment.rdata" ) )
      }

  load( file=paste(sep="", Sequencing_Statistics_dir, dataframe_file_name, "_pairwiseAlignment.rdata" ) )
  
  uncut_fraction <- sum( pairwiseAlignment_vec@pattern@range@width==147)/length(pairwiseAlignment_vec)
  cut_percentage <- round(100*(1- uncut_fraction ), digits=1)
  left_cut_frequency_vec <- right_cut_frequency_vec <- rep(0, 72)
  names(left_cut_frequency_vec ) <-  1:67
  names( right_cut_frequency_vec ) <- 147:81
  left_frequency_vec <- table( pairwiseAlignment_vec@subject@range@start)
  left_cut_frequency_vec[names(left_frequency_vec)] <- left_frequency_vec
  right_frequency_vec <- table( pairwiseAlignment_vec@subject@range@start+ pairwiseAlignment_vec@subject@range@width-1)
  right_cut_frequency_vec[names(right_frequency_vec )] <- right_frequency_vec
  max_count_freq <- max( left_frequency_vec [-1] , right_cut_frequency_vec[-1])
  # Plot for Extended Data Panel S1d
  graph_filename=paste( sep="", Sequencing_Statistics_dir, dataframe_file_name , "_Cut_Histogram.png" )
  CairoPNG(width = 450, height = 300 , file=  graph_filename )
  plot( c(0,max_count_freq),ylab="Cut Counts", xlab="Cut Position (bp)", xlim=c(1,147), type="l", col="white")
  left_cut_percentage <- round(100*(1-left_cut_frequency_vec [1]/sum(left_cut_frequency_vec )), 1)
  right_cut_percentage <- round(100*(1-right_cut_frequency_vec [1]/sum(right_cut_frequency_vec )), 1)
  if( is.element(sequence_number, c("203","204", "335", "336", "337","338") )){
      segments (y0=0, y1 = 0,  x0=32, x1=39, col="grey", lwd=4)
      segments (y0=0, y1 = 0,  x0=114, x1=121, col="grey", lwd=4)
    }
  if( is.element(sequence_number, c("337") )){segments (y0=0, y1 = 0,  x0=32, x1=39, col="red", lwd=4)}
  if( is.element(sequence_number, c( "338") )){segments (y0=0, y1 = 0,  x0=114, x1=121, col="red", lwd=4)}
  title(paste( "Template", sequence_number, ": Left Cuts", left_cut_percentage , "% VS Right Cut", right_cut_percentage, "%"))
  points(x=names(left_cut_frequency_vec [-1]), y=left_cut_frequency_vec [-1] , type="l", col="black", lwd=1)
  points(x=names(right_cut_frequency_vec [-1]), y=right_cut_frequency_vec [-1] , type="l", col="black", lwd=1)
  dev.off()

  histogram_filename <- paste0(     Parent_dir,"Histograms_to_upload/", dataframe_file_name, "_read_ends.txt"  )
  write.table( c(left_cut_frequency_vec[!is.na(names(left_cut_frequency_vec))],
  rev(right_cut_frequency_vec[!is.na(names(right_cut_frequency_vec))])),
  row.names=T, col.names = T, quote=F,       sep = "\t",
  file=   histogram_filename
  )
  
  if( dataframe_file_name=='F3-203-1') { 
    system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_left_top.png") ) 
    system( paste0( "cp ",   histogram_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_left_top.txt") ) 
    }
  if( dataframe_file_name=='G3-203-2'){  
    system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_left_bottom.png") ) 
    system( paste0( "cp ", histogram_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_left_bottom.txt") ) 
  }
  if( dataframe_file_name=='F2-338-1') {
    system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_right_bottom.png") ) 
    system( paste0( "cp ", histogram_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_right_bottom.txt") ) 
    }
  if( dataframe_file_name=='E2-337-2') { 
    system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_right_top.png") ) 
    system( paste0( "cp ", histogram_filename, " " , Graphs_dir, "Ext_Data_Fig_1d_right_top.txt") ) 
    }
  }
}
#
####################################################################################
################################ Analysis of PIONEAR-seq libraries #########################
#####################################################################################

if(1)
  {
  
  paired_file_vec <-   list.files( Paired_dir  )
  WANC_filename_vec <- paired_file_vec[ grepl( 'WANC', paired_file_vec)]
  read_number_matrix <- matrix(data=NA, nrow=4, ncol=length(WANC_filename_vec), dimnames=list( template_name_vec,    WANC_filename_vec ) )
  
  for( WANC_filename in WANC_filename_vec )
  {
    paired_file <- paste( sep="", Paired_dir, WANC_filename )
    dataframe_file_name <- sub( ".assembled.fastq.gz","",  substring(paired_file, regexpr("Paired_files/", paired_file) + 13)) 
    paired_file_fastq <- readFastq(  paired_file )
    if( length(paired_file_fastq)>10^6 ){ paired_file_fastq <- paired_file_fastq[sample(length(paired_file_fastq),10^6)]}
    paired_file_fastq_unique_reads <- unique(paired_file_fastq@sread)
    paired_file_fastq <- paired_file_fastq[match( paired_file_fastq_unique_reads, paired_file_fastq@sread)]
    paired_file_fastq <- trimTailw( paired_file_fastq, 1, "?", halfwidth=147) 
    paired_file_fastq <- paired_file_fastq[width( paired_file_fastq@sread )<=147  ]
    paired_file_fastq   <-  paired_file_fastq[width( paired_file_fastq@sread )>80  ]
    write( paste(paired_file, " has ", length(paired_file_fastq), "unique reads between 81 and 147 bp \n" ), stdout()) 
    for(index_template in 1:length(template_name_vec))
    {
      template_name <- template_name_vec[index_template]
      tile_sequence_start_vec <- tile_sequence_start_vec_list[[template_name]]
      template_sequence <- DNAString(sequence_template_vec[template_name])
      sublibrary_dir <- paste(sep="", Sequencing_Statistics_dir, dataframe_file_name,"_" , template_name,"/")
      if( !dir.exists( sublibrary_dir) ) { system( paste( "mkdir ",   sublibrary_dir ) ) }
      if(0) # Computationally expensive: this part aligns the digested reads to the original template 
      {
        FW_pairwiseAlignment_vec <-   pairwiseAlignment(pattern=   paired_file_fastq@sread,
                                     subject=  template_sequence, type="overlap",  substitutionMatrix=mapping, gapOpening=1000, gapExtension=0)
        RC_pairwiseAlignment_vec  <-pairwiseAlignment(pattern=  reverseComplement( paired_file_fastq@sread),
                                     subject=  template_sequence, type="overlap",   substitutionMatrix=mapping, gapOpening=1000, gapExtension=0)
        pairwiseAlignment_vec <- FW_pairwiseAlignment_vec
        RC_lab <-  which(FW_pairwiseAlignment_vec@score<  RC_pairwiseAlignment_vec@score)
        pairwiseAlignment_vec  [ RC_lab ]  <-  RC_pairwiseAlignment_vec[ RC_lab ]   
        pairwiseAlignment_vec <- pairwiseAlignment_vec [pairwiseAlignment_vec@score>(pairwiseAlignment_vec@pattern@range@width-20)]
        pairwiseAlignment_vec <- pairwiseAlignment_vec[which(sapply(pairwiseAlignment_vec@pattern@mismatch,length)>0)]
        mismatched_strech_vec <- sapply(pairwiseAlignment_vec@subject@mismatch,max) - sapply(pairwiseAlignment_vec@subject@mismatch,min)
        pairwiseAlignment_vec <- pairwiseAlignment_vec[which(  mismatched_strech_vec <20&mismatched_strech_vec >3)]
        no_error_in_fw_primer_lab <- sapply(pairwiseAlignment_vec@subject@mismatch,min)>=tile_sequence_start_vec[1]
        no_error_in_rv_primer_lab <- sapply(pairwiseAlignment_vec@subject@mismatch,max)<tail(tile_sequence_start_vec,1)+20
        pairwiseAlignment_vec <- pairwiseAlignment_vec [which(no_error_in_fw_primer_lab&no_error_in_rv_primer_lab) ]
        save( pairwiseAlignment_vec, file=paste(sep="", sublibrary_dir, dataframe_file_name,"_" , template_name ,"_pairwiseAlignment.rdata" ) )
        }
      load( file=paste(sep="", sublibrary_dir, dataframe_file_name,"_" , template_name ,"_pairwiseAlignment.rdata" ) )
      read_number_matrix[ template_name, WANC_filename]  <- length(pairwiseAlignment_vec)
      min_mismatch_position_vec <- sapply(pairwiseAlignment_vec@subject@mismatch,min)
      max_mismatch_position_vec <- sapply(pairwiseAlignment_vec@subject@mismatch,max)
      start_position_vec <- pairwiseAlignment_vec@subject@range@start
      end_position_vec<- pairwiseAlignment_vec@subject@range@start +pairwiseAlignment_vec@subject@range@width-1
      left_end_cut_read_lab <- which(max_mismatch_position_vec-start_position_vec<20)
      right_end_cut_read_lab <-which(end_position_vec-min_mismatch_position_vec<20)
      assigned_tile_lab <- rep( F, length(pairwiseAlignment_vec)) 
      bp_frequency_matrix <-  matrix( data=0, ncol=147, nrow=length( tile_sequence_start_vec ) )
      rownames( bp_frequency_matrix) <- tile_sequence_start_vec
      colnames( bp_frequency_matrix) <- 1:147
      bp_count_matrix <-bp_frequency_matrix
      for( index_tile_start in 1:length( tile_sequence_start_vec ) )
      {
        tile_start <- tile_sequence_start_vec[ index_tile_start ]
        tile_end <- tile_start+19
        tile_lab <- (min_mismatch_position_vec >= tile_start) & (tile_end>=max_mismatch_position_vec)
        tile_lab[right_end_cut_read_lab] <- abs(min_mismatch_position_vec[right_end_cut_read_lab] - tile_start-1.5) < 2
        tile_lab[left_end_cut_read_lab] <- abs(tile_end -1 - max_mismatch_position_vec[left_end_cut_read_lab]-1.5) < 2
        sublibrary_pairwiseAlignment_vec <- pairwiseAlignment_vec[which(tile_lab)]
        assigned_tile_lab <- tile_lab|  assigned_tile_lab
        start_bp_count_vec <- rep( 0, 147)
        names(start_bp_count_vec) <- 1:147  
        start_bp_freq_vec <- start_bp_count_vec
        FW_start_bp_vec <- sublibrary_pairwiseAlignment_vec@subject@range@start
        FW_count_start_bp_table <- table(FW_start_bp_vec)
        start_bp_count_vec[names(FW_count_start_bp_table)] <-  FW_count_start_bp_table
        start_bp_freq_vec[names(FW_count_start_bp_table)] <-  FW_count_start_bp_table/sum(FW_count_start_bp_table)
        RC_start_bp_vec <- sublibrary_pairwiseAlignment_vec@subject@range@start+sublibrary_pairwiseAlignment_vec@subject@range@width-1
        RC_count_start_bp_table <- table(RC_start_bp_vec)
        start_bp_count_vec[names(RC_count_start_bp_table)] <-  RC_count_start_bp_table
        start_bp_freq_vec[names(RC_count_start_bp_table)] <-  RC_count_start_bp_table/sum(RC_count_start_bp_table)
        bp_count_matrix[ index_tile_start, ] <- start_bp_count_vec
        
        bp_frequency_matrix[ index_tile_start, ] <- start_bp_freq_vec
        uncut_fraction <- sum(sublibrary_pairwiseAlignment_vec@subject@range@width==147)/length(sublibrary_pairwiseAlignment_vec)
        cut_percentage <- round(100*(1- uncut_fraction ), digits=1)
        bp_count_matrix[index_tile_start,tile_start:tile_end] <-         bp_frequency_matrix[index_tile_start,tile_start:tile_end] <- NA
      }
      rownames(bp_count_matrix) <- paste0( "Tile_from_", tile_sequence_start_vec ,"_to_", tile_sequence_start_vec+19 , "_bp")
      count_matrix_file=paste(sep="", sublibrary_dir, dataframe_file_name,"_" , template_name, "_bp_count_matrix.txt")
      write.table(rbind( Cut_position_bp = 1:147, bp_count_matrix) ,row.names=T, col.names = F, quote=F,       sep = "\t",
                   file=   count_matrix_file)
      system( paste( "cp" ,  count_matrix_file, "Histograms_to_upload/."))
      
      save(bp_frequency_matrix, file=paste(sep="", sublibrary_dir, dataframe_file_name,"_" , template_name, "_bp_frequency_matrix.Rdata"))
      load( file=paste(sep="", sublibrary_dir, dataframe_file_name,"_" , template_name, "_bp_frequency_matrix.Rdata"))
      max_frequency=10
    
      # Plot for Extended Data Panels S1e-f
      dataframe_template_file_name <- paste( sep="",    dataframe_file_name,"_" , template_name)
      graph_filename = paste( sep="",  sublibrary_dir ,       dataframe_template_file_name , "_Heatmap.png") 
      CairoPNG(width = 800, height = 800 ,file=   graph_filename)
      heatmap.2( 100*bp_frequency_matrix[,2:146],symm = FALSE, dendrogram = "none",Rowv=F, Colv=F, 
              breaks=seq(0,max_frequency,max_frequency/100),
               col=c( grey.colors(100,start=1, end=0, rev=T) ) ,
              density.info=c("none"), na.color="blue", trace="none", main= template_name)
      dev.off()
      
      if( dataframe_template_file_name=='D1-WANC-1-MNase_NRCAM')  system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1e_bottom.png") )
      if( dataframe_template_file_name=='H1-WANC-2-MNase_NRCAM')  system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1e_top.png") )
      if( dataframe_template_file_name=='H1-WANC-2-MNase_CX3')  system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1f_left.png") )
      if( dataframe_template_file_name=='H1-WANC-2-MNase_ALB1')  system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1f_center.png") )
      if(  dataframe_template_file_name=='H1-WANC-2-MNase_W601')  system( paste0( "cp ", graph_filename, " " , Graphs_dir, "Ext_Data_Fig_1f_right.png") )
    }
    count_matrix_files=paste(sep="",  dataframe_file_name ,"_" , template_name_vec,"_bp_count_matrix.txt", collapse = " ")
    setwd("Histograms_to_upload/")
    system(paste0( "tar -cf ", dataframe_file_name, "_read_ends.tar ",   count_matrix_files) )
    system( paste( "rm" ,   count_matrix_files) )
    setwd("..")
    }

  save( read_number_matrix,    file=paste(sep="", Sequencing_Statistics_dir, "NCP_Input_Library_read_number_matrix.rdata" ) )
  load(     file=paste(sep="", Sequencing_Statistics_dir, "NCP_Input_Library_read_number_matrix.rdata" ) )
  
  percentage_read_matrix <- 100*t(read_number_matrix)/colSums(read_number_matrix)
  save( percentage_read_matrix,    file=paste(sep="", Sequencing_Statistics_dir, "NCP_Input_Library_percentage_read_matrix.rdata" ) )
      
  lib_percentage_mat=as.data.frame(percentage_read_matrix)
  lib_percentage_mat$exp = rownames(lib_percentage_mat)
  lib_percentage_mat_long =  pivot_longer(lib_percentage_mat, cols=-exp)
  lib_percentage_mat_long$group = grepl( "MNase",lib_percentage_mat_long$exp)
  avg_lib_1 = lib_percentage_mat_long %>% filter(grepl( "MNase",lib_percentage_mat_long$exp)) %>% group_by(name) %>% summarise(avg=mean(value)) 
  
 # Plot for panel 3b
  ggplot(avg_lib_1, aes(x=reorder(name,-avg), y=avg, fill=name))+geom_bar(stat="identity", show.legend = F)+theme_classic()+
    geom_point(data=lib_percentage_mat_long %>% filter(group), aes(x=name, y=value, shape=exp), color="black",show.legend = F, size=2)+
    labs(x="Library template sequence", y="Protected reads (%)")+
    scale_fill_brewer(palette = "Set2")+theme(text=element_text(size=7), title=element_text(size=7))
  ggsave(paste(sep="", Sequencing_Statistics_dir, "template_MNase_protection_barplots.png"), width=60, height=80, units = "mm")
  ggsave(paste(sep="", Graphs_dir, "Fig_3b.png"), width=60, height=80, units = "mm")
  dev.off()
  
  lib_percentage_mat_long$group = case_when(
      grepl("1-MNase", lib_percentage_mat_long$exp) ~ "After digestion replicate 1 (   )",
      grepl("2-MNase", lib_percentage_mat_long$exp) ~ "After digestion replicate 2 (   )",
      grepl("1-Input", lib_percentage_mat_long$exp) ~ "Before digestion replicate 1",
      grepl("2-Input", lib_percentage_mat_long$exp) ~ "Before digestion replicate 2")
   
   avg_lib_2 = lib_percentage_mat_long %>% group_by(name, group) %>% summarise(avg=mean(value))
   lib_percentage_mat_long$label = paste0(round(lib_percentage_mat_long$value, digits=2), "%")

   # Plot for panel S6a   
   ggplot(lib_percentage_mat_long,
          aes(x=factor(group, levels=rev(c("Before digestion replicate 1",
                                           "Before digestion replicate 2",
                                           "After digestion replicate 1 (   )",
                                           "After digestion replicate 2 (   )"))),
          y=value/100, fill=name, label=label))+ geom_bar(stat="identity")+
          geom_text(size = 3, position = position_stack(vjust = 0.5))+
          scale_fill_brewer(palette = "Set2")+ scale_y_continuous(labels = scales::percent)+
          theme_classic()+labs(x="", y="Percentage of reads in library", fill="Template")+
          theme(axis.text.x = element_text(angle=90, vjust=0.5))+coord_flip()
   ggsave(paste(sep="", Sequencing_Statistics_dir, "template_MNase_before_after.png"), width=180, height=70, units = "mm")
   ggsave(paste(sep="", Graphs_dir, "Ext_Data_Fig_6a.png"), width=180, height=70, units = "mm")
   dev.off()
   }   

############################ END #################################
    
    
    
    
    
    
    
         
        
        
           
      


