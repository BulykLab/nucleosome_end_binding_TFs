############################################################################
# To initialize the program for any of the subsequent parts #
############################################################################

if(TRUE)
  {  
  rm(list=ls())
  
  # - R Libraries upload
  if(1)
    {
    options(java.parameters="-Xmx16g")     # to try to avoid  "java.lang.OutOfMemoryError: Java heap space"
    library(qdap)
    library(ShortRead)
    library(rJava)
    library(SELEX)  
    library(stringr)
    library(reshape2)
    library(Cairo)  
    library(pheatmap)  
    library(gplots) 
    }
  
  # Directories
  if(1)
    { 
    
      Parent_dir <-  $YOUR_PARENTAL_DIR  
    #for example 
    # Parent_dir <- "/Volumes/My_Book_3/PIONEAR_seq_Project_Back_Up_2022/Manuscript/preparation/GEO_repository/PIONEAR_seq_datasets_to_upload/"
    setwd(Parent_dir)
    
    
    Raw_sequencing_data_dir <- paste( sep="" ,Parent_dir, "Raw_sequencing_files/" )

    Paired_sequencing_file_input_dir <-  paste(sep="" ,  Parent_dir, "Paired_files/")
    if (! dir.exists(Paired_sequencing_file_input_dir )) {system( paste( "mkdir ",  Paired_sequencing_file_input_dir  ) ) }
    
    Preprocessed_files_dir <- paste( sep="",  Parent_dir, "preprocessed_files/" )
    if (! dir.exists(  Preprocessed_files_dir )) { system( paste( "mkdir ",   Preprocessed_files_dir  ) ) }
       
    Verbose_dir <- paste(sep="",    Preprocessed_files_dir, "Verbose/")
    if( !dir.exists(Verbose_dir)){system( paste( "mkdir ",  Verbose_dir ) )    }      
   
    isotile_16bp_file_dir <- paste(sep="",    Preprocessed_files_dir, "isotile_16bp_files/")
    if( !dir.exists(isotile_16bp_file_dir)){system( paste( "mkdir ",  isotile_16bp_file_dir ) )    }      
    
    isotile_16bp_downsized_file_dir <- gsub( "16bp_files/", "16bp_downsized_files/", isotile_16bp_file_dir)
    if( !dir.exists(isotile_16bp_downsized_file_dir)){  system( paste( "mkdir", isotile_16bp_downsized_file_dir) )}
    
    isotile_16bp_positioned_file_dir <- paste(sep="",    Preprocessed_files_dir, "isotile_16bp_positioned_files/")
    if( !dir.exists(isotile_16bp_positioned_file_dir)){system( paste( "mkdir ",  isotile_16bp_positioned_file_dir ) )    }      
      
 
    Kmer_Enrichment_dir <- paste(sep="", Parent_dir,  "Kmer_Enrichment_for_isotile_files/")
    if (! dir.exists(  Kmer_Enrichment_dir  )) { system( paste( "mkdir ",   Kmer_Enrichment_dir ) )}
    
    Graphs_for_Kmer_enrichment_dir <- paste( sep="",   Parent_dir , "Graphs_for_Kmer_enrichment_for_isotile_files/" )
    if (! dir.exists(      Graphs_for_Kmer_enrichment_dir )) { system( paste( "mkdir ",    Graphs_for_Kmer_enrichment_dir ) )}
    
    Graphs_for_Kmer_profiles_dir <- paste( sep="",   Parent_dir , "Graphs_for_Kmer_positional_preference_for_isotile_files/" )
    if (! dir.exists( Graphs_for_Kmer_profiles_dir )) { system( paste( "mkdir ",  Graphs_for_Kmer_profiles_dir ) )}
    
    Graphs_dir <-paste( sep="", Parent_dir, "Graphs/" )
    if (! dir.exists(  Graphs_dir )) { system( paste( "mkdir ",   Graphs_dir) ) }
    
    }  
  
  # - Variables 
  if(1)
    {
    ALB1_sequence <-  toupper("CTGCTCTGTCAGCAGGGCACTGTACTTGCTGATACCAGGGAATGTTTGTTCTTAAATACCATCATTCCGGACGTGTTTGCCTTGGCCAGTTTTCCATGTACATGCAGAAAGAAGTTTGGACTGATCCATACAGTCCTCTGCCTTTTC")
    CX3_sequence <- toupper("GGGCCTCTCGGCTGCTGATCTTCAGCTGGTTGCTGAGAGTTGCAGCATTGCTGAGTCTTAGCAATGGATACTTgCCGATTCCCCTCACAAAAATAGGTCAGTCTGTCTGGCTAGTTCTGTACTTGCAGACACAGGGCATGTGGGGTT")
    NRCAM_sequence <- toupper("ACTTCTGAAACAGATGACTCCCAGCAGCTGCTGCCTGTGGCgCACAGGGCTaCCTGCCCTGCATGACAGCTGCACATCACATCCTGTGGTCATACTACTTCAGCCGCTTCTACGGCCAGATACAAAAGTGGGTGGcGAACATAGGCA")
    W601_sequence <- toupper("ATCGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACcTGTCAGATATATAGATCCGAT")
      
    tiled_ALB1_start_vec <-  seq(21,102,3 )
    tiled_CX3_start_vec <- seq(18,108,3 )
    tiled_NRCAM_start_vec <-  seq(26,107,3 )
    tiled_W601_start_vec <- seq(22,100,3 )
      
    sequence_template_vec <- c( ALB1_sequence , CX3_sequence , NRCAM_sequence , W601_sequence )
    tile_sequence_start_vec_list <- list(tiled_ALB1_start_vec, tiled_CX3_start_vec, tiled_NRCAM_start_vec, tiled_W601_start_vec)
    
    template_name_vec <- c(  "ALB1", "CX3" , "NRCAM" , "W601" )
    names( sequence_template_vec) <- names( tile_sequence_start_vec_list) <-    template_name_vec
      
    TF_kmer_assigned_length_vec <-  c( "ASCL1"=8, "CEBPB"=8, "FOXA1"=9, "GATA1"=8, "GATA4"=8, "GST"=8, "KLF4"=9, "OCT4"=10, "PAX7"=8, "PU1"  =9, "RUNX2"=8, "ZNF143"=9)
    TF_name_vec <- names(TF_kmer_assigned_length_vec)
    
    mapping <- diag(4)
    dimnames(mapping) <- list(DNA_BASES, DNA_BASES)
    
    
    corr_method <- "spearman"
    }
  
  }

############################################################################
######## To Merge Paired-end reads via PEAR # ###################
############################################################################

if(1)
  {       
  PEAR_executable <- $YOUR_PEAR_EXECUTABLE #For example 
  # for example
  # PEAR_executable <- "/Applications/pear/src/pear"
  R1_R2_fastq_file_vec <- list.files(path =     Raw_sequencing_data_dir  )
  R1_suffix <- "_R1_001.fastq.gz"
  R2_suffix <- "_R2_001.fastq.gz"
  R1_fastq_file_vec <-  R1_R2_fastq_file_vec[ grep(  R1_suffix , R1_R2_fastq_file_vec)]
  R2_fastq_file_vec <-  R1_R2_fastq_file_vec[ grep(  R2_suffix , R1_R2_fastq_file_vec)]
  #Check that all the files are there
  match( gsub( R1_suffix ,"", R1_fastq_file_vec), gsub( R2_suffix ,"",R2_fastq_file_vec))
  # Check no extra files are there (must be empty)
  which(!grepl(  R1_suffix , R1_R2_fastq_file_vec)&!grepl(  R2_suffix , R1_R2_fastq_file_vec))
  
  fastq_file_vec <- gsub( "_R1_001.fastq.gz", "", R1_fastq_file_vec)
  for( index_fastq_file in 1:length(fastq_file_vec))
  {  
    fastq_file <- fastq_file_vec[index_fastq_file]
     
    zipped_R1_fastq_file_address <- paste(sep="",  Raw_sequencing_data_dir, fastq_file, "_R1_001.fastq.gz")
    zipped_R2_fastq_file_address <- gsub("_R1_", "_R2_", zipped_R1_fastq_file_address )
    system( paste( "gunzip " ,zipped_R1_fastq_file_address ) )
    system( paste( "gunzip " , zipped_R2_fastq_file_address) )
    
    unzipped_R1_fastq_file_address <- gsub("_001.fastq.gz", "_001.fastq", zipped_R1_fastq_file_address )
    unzipped_R2_fastq_file_address <-  gsub("_R1_", "_R2_", unzipped_R1_fastq_file_address )
   
    Paired_fastq_file_address <- paste(sep="",  Paired_sequencing_file_input_dir, fastq_file)
    
    PEAR_command <-paste( PEAR_executable , " -f ",  unzipped_R1_fastq_file_address   ,
                           " -r ",   unzipped_R2_fastq_file_address , " -o ",  Paired_fastq_address )
    write(  PEAR_command , stdout())     
    system( PEAR_command )
    system( paste( "gzip " , unzipped_R1_fastq_file_address ) )
    system( paste( "gzip " , unzipped_R2_fastq_file_address ) )
    system( paste("gzip", unzipped_Paired_template_filename ))
  }
  
  # Remove extra files
  system( paste(sep="", "rm ", Paired_sequencing_file_input_dir,"*.discarded.fastq") )
  system( paste(sep="", "rm ", Paired_sequencing_file_input_dir,"*.unassembled*fastq") )
 
   # Optional: 1 - Use the cat command to merge the inputs 
  
   #           2 - Zip files if you want to save space ("  system( paste("gzip", unzipped_Paired_template_filename )))
  #                by change replcing if(0) to if(01) 
  if(1)
    {
    setwd( Paired_sequencing_file_input_dir)
    Paired_sequencing_file_input_vec <- list.files()
    for( index_file in 1:length(Paired_sequencing_file_input_vec)) 
      {
      Paired_template_filename <-  Paired_sequencing_file_input_vec [index_file]
      system( paste("gzip", Paired_template_filename ))
      }
    }
  write(  "PAIRING END" , stdout()) 
  }

##########################################################################
###  To align and organize assembled reads into tile files  ####
##############################################################################

if(1)
  {
  # To aliogn the paired files to templates and subdivide them according to the random tiles
  setwd( Paired_sequencing_file_input_dir)
  Paired_sequencing_file_input_vec <- list.files()
  Paired_sequencing_file_size_vec <- file.size( Paired_sequencing_file_input_vec)
  
  names(  Paired_sequencing_file_size_vec ) <- names( Paired_sequencing_file_input_vec ) <- gsub(".assembled.fastq.gz","",  Paired_sequencing_file_input_vec)

  Ordered_Paired_sequencing_file_input_vec <- Paired_sequencing_file_input_vec[order( Paired_sequencing_file_size_vec, decreasing = F)]
  Paired_sequencing_file_size_vec [ names(Ordered_Paired_sequencing_file_input_vec )]
  random_region_size <- 20
 
   read_limit <-  10^6 # this parameter can vary accordingly to the CPU to make the program to process reads in subsets sequentially
  
  for( index_file in  1:length( Ordered_Paired_sequencing_file_input_vec ))
    {
    Paired_sequencing_file <- Ordered_Paired_sequencing_file_input_vec[index_file ]
    template_name <-  strsplit(  names(Paired_sequencing_file) , "-")[[1]][2]
    tile_sequence_start_vec <- tile_sequence_start_vec_list[[template_name]]
    template_sequence <- sequence_template_vec[template_name]
      
    Single_file_dir <- paste(sep="",    Verbose_dir , names(   Paired_sequencing_file),"/")
      
    if( !dir.exists(Single_file_dir)){system( paste( "mkdir ",  Single_file_dir  ) )    } 
    Tile_file_dir <-  paste( sep="", Single_file_dir, "tile_files/")
    if( !dir.exists(Tile_file_dir)) { system( paste( "mkdir ",  Tile_file_dir ) ) } 
    Tile_suffix_vec <- paste(sep="", names(   Paired_sequencing_file), "_tile_from_",   tile_sequence_start_vec ,  "_to_",  tile_sequence_start_vec + random_region_size-1,  ".fastq.gz")
    Tile_file_outputnames_vec <-  paste( sep="", Tile_file_dir,   Tile_suffix_vec )
   
    read_file_fq <- readFastq(Paired_sequencing_file )
   
    #Make reads unique and trim the UMIs
    duplicated_read_file_freq <-  duplicated(sread(read_file_fq))
    trim_1_read_file_fq_full <-read_file_fq[!duplicated_read_file_freq] 
    initial_random_divergent_tag_size <- 6
    trim_1_read_file_fq_full <- narrow(  trim_1_read_file_fq_full, start=initial_random_divergent_tag_size+1 )
    remove( read_file_fq , duplicated_read_file_freq )

    #Filter the reads for sequencing quality and correct length
    rng <- trimEnds(sread(trim_1_read_file_fq_full ), "N", relation="==", ranges=TRUE)
    trim_1_read_file_fq <-  narrow(trim_1_read_file_fq_full, start(rng), end(rng))
    trim_1_read_file_fq <- trimTailw(  trim_1_read_file_fq, 1, "#", halfwidth=147) 
    trim_1_read_file_fq <- trim_1_read_file_fq[width(  trim_1_read_file_fq@sread) ==147  ]
    trim_1_read_file_fq <- trim_1_read_file_fq[ trim_1_read_file_fq@sread !=template_sequence  ]
    remove( trim_1_read_file_fq_full ,rng)
      
    read_limit_count <- ceiling(length(trim_1_read_file_fq)/read_limit)
    read_limit_lab <- ((read_limit_count-1)*read_limit+1):length(trim_1_read_file_fq)
    writing_mode <- "w"
        
    # Read alignment to the template
    for(   index_read_limit in read_limit_count:1)
      {
      
      pairwiseAlignment_vec <-   pairwiseAlignment(pattern=    trim_1_read_file_fq@sread[read_limit_lab],
                                                        subject=  template_sequence, type="overlap",  substitutionMatrix=mapping, gapOpening=1000, gapExtension=0)
       
      good_score_lab <- pairwiseAlignment_vec@score>(147-random_region_size)
      max_mismatch_position_vec <- sapply(pairwiseAlignment_vec@subject@mismatch,max)
      min_mismatch_position_vec <- sapply(pairwiseAlignment_vec@subject@mismatch,min)
      remove( pairwiseAlignment_vec )
      
      mismatched_strech_vec <- max_mismatch_position_vec - min_mismatch_position_vec
      long_stretch_lab <-   mismatched_strech_vec <= random_region_size & mismatched_strech_vec > (random_region_size-4)
        
      no_error_in_fw_primer_lab <- min_mismatch_position_vec >= tile_sequence_start_vec[1]
      no_error_in_rv_primer_lab <- max_mismatch_position_vec < tail(tile_sequence_start_vec,1)+random_region_size
      
      full_lab <- good_score_lab & long_stretch_lab & no_error_in_fw_primer_lab  & no_error_in_rv_primer_lab
      trim_2_read_file_fq <- trim_1_read_file_fq[read_limit_lab [full_lab]]
      min_mismatch_position_vec <- min_mismatch_position_vec [full_lab]
      max_mismatch_position_vec <- max_mismatch_position_vec [full_lab]
      assigned_tile_lab <- rep( F, sum(full_lab)) 
      
      # Position of the internal random tiles
      for( index_tile_start in 1:length( tile_sequence_start_vec ) )
        {
        Tile_file_outputname <- Tile_file_outputnames_vec[ index_tile_start]
        tile_start <- tile_sequence_start_vec[ index_tile_start ]
        tile_end <- tile_start+ random_region_size -1
        tile_lab <- (min_mismatch_position_vec >= tile_start) & (tile_end>=max_mismatch_position_vec)
        tile_lab <- tile_lab&!assigned_tile_lab
        subset_read_file_fq <- trim_2_read_file_fq[tile_lab]
        assigned_tile_lab <- tile_lab|  assigned_tile_lab
        tile_read_file_fq <- narrow( subset_read_file_fq, start = tile_start, end = tile_end )
        writeFastq( tile_read_file_fq, Tile_file_outputname, writing_mode, compress=T)
        remove( subset_read_file_fq , tile_read_file_fq)
          }
        remove(  trim_2_read_file_fq)
        writing_mode <- "a"
        read_limit_lab <-(1:read_limit)+read_limit*(index_read_limit-2)
      } 
        remove( trim_1_read_file_fq )
      
  }
  
}


################################################################
###  To check and remove elements excessively duplicated  ####
################################################################

if(1)
{
  setwd( Verbose_dir)
  Verbose_dir_vec <- list.files( Verbose_dir)
  GST_Verbose_dir_vec <- Verbose_dir_vec[ grepl("GST|INPUT",Verbose_dir_vec )]
  
  max_vec_array_list <- list()
  for (template_name in template_name_vec){

    template_Verbose_dir_vec <- GST_Verbose_dir_vec[grep( template_name,   GST_Verbose_dir_vec)]
    max_vec <- c()

    for( index_file in 1:length(  template_Verbose_dir_vec ))
    {
      Single_file_name <-   template_Verbose_dir_vec[index_file ]
      Single_file_dir <- paste(sep="",    Verbose_dir , Single_file_name,"/")
      Tile_file_dir <-  paste( sep="", Single_file_dir, "tile_files/")
    
     # to check tiles for excessive duplicates 
      if(1){
        setwd( Tile_file_dir)
        file_to_check_vec <- list.files()
        for( index_tile in 1:length(  file_to_check_vec  ) )
        {
          file_to_check <-      file_to_check_vec  [  index_tile]
          Tile_read_fq <-   readFastq(  file_to_check )
          sread_table <- table( sread( Tile_read_fq))
          max_vec[index_tile] <- max( sread_table )
          }
        }
      names( max_vec ) <- file_to_check_vec
      if(index_file==1)    {max_vec_array  <- max_vec  }
      if(index_file>1)   {max_vec_array <- cbind(max_vec_array  ,max_vec  )}
      }
  
    colnames(max_vec_array) <-   template_Verbose_dir_vec
  
    rownames(max_vec_array) <- gsub( template_Verbose_dir_vec[1], "max", gsub(".fastq.gz","", rownames(max_vec_array) ) )
  
    max_vec_array_list[[template_name]] <- max_vec_array
    }

  names(max_vec_array_list) <- template_name_vec
  CairoPNG(  1200,600,file=     paste0( Graphs_for_Kmer_enrichment_dir, "Boxplot_GST_duplicated_tile_read_statistics_between_templates.png")   , pointsize = 40/2 ,bg = "white" , res=72*2)
  par(mar = c(4,4, 1, 1) ) 
  boxplot(  max_vec_array_list,log='y',  Main="Max Duplicated Tile (#)",ylab="Counts")
  dev.off()

  boxplot(  max_vec_array_list,log='y',   Main="Duplicated tile read count",ylab="Counts (Log10)")
  # look at these variable: 
  max_vec_array_list$W601
  max_vec_array_list$ALB1
  max_vec_array_list$CX3
  max_vec_array_list$NRCAM

  # after manually checking, tiles with excessive duplicated reads are
  # ALB1 tile from 54bp to 73bp: "AAATACCATTATGCGCATAG"
  # ALB1 tile from 60bp to 79bp:       "CATTATGCGCATAGTGTTTG" ,
  # ALB1 tile from 63bp to 82b:          "TATGCGCATAGTGTTTGCCT"
  # W601 tile from 55bp to 74bp:  "GCTTAAACTATGCGCATAGC"
  # W601 tile from 61bp to 80bp:        "ACTATGCGCATAGCTGTCCC"

  ALB1_tile_to_remove_vec <-  c( tile_from_54_to_73    ="AAATACCATTATGCGCATAG", 
                               tile_from_60_to_79    =      "CATTATGCGCATAGTGTTTG" ,
                               tile_from_63_to_82    =         "TATGCGCATAGTGTTTGCCT" )
  W601_tile_to_remove_vec <-  c( tile_from_55_to_74="GCTTAAACTATGCGCATAGC", 
                               tile_from_61_to_80    =  "ACTATGCGCATAGCTGTCCC")

  # to clean tile files from tiles with excessive duplicated reads
  Clean_Verbose_dir_vec <- Verbose_dir_vec[grepl( "W601|ALB1", Verbose_dir_vec)]
  Clean_Verbose_dir_vec <- Clean_Verbose_dir_vec[!grepl( "INPUT", Clean_Verbose_dir_vec)]
  for( index_file in 1:length(   Clean_Verbose_dir_vec ))
  {
    Single_file_name <-   Clean_Verbose_dir_vec[index_file ]
    Single_file_dir <- paste(sep="",    Verbose_dir , Single_file_name,"/")
    Tile_file_dir <-  paste( sep="", Single_file_dir, "tile_files/")
    setwd( Tile_file_dir)
    file_to_isotile_vec <- list.files()
    preclean_tile_file_dir <-  "../precleaned_tile_files/"
    if(  grepl("W601",  Single_file_name))
      {  
      if( !dir.exists(preclean_Tile_file_dir)) {system( paste( "mkdir"  , preclean_Tile_file_dir  ))}
      tile_to_clean_vec <-  c("55_to_74", "61_to_80" )
      file_to_clean_vec <- file_to_isotile_vec[as.logical(rowSums(
      sapply( tile_to_clean_vec ,  grepl,file_to_isotile_vec  )))]
      for( file_to_clean in file_to_clean_vec)
        {
        if(! is.element( file_to_clean , list.files(preclean_tile_file_dir))){system(  paste("mv", file_to_clean,  preclean_tile_file_dir) ) }
        W601_tile_to_remove <-  W601_tile_to_remove_vec[sapply( tile_to_clean_vec ,grepl,  file_to_clean )]
        preclean_Tile_read_fq <-   readFastq(paste0( preclean_tile_file_dir,file_to_clean))
        pairwiseAlignment_vec <-   pairwiseAlignment(pattern=    sread(preclean_tile_read_fq),
                                                   subject=   W601_tile_to_remove , 
                                                   type="overlap",  substitutionMatrix=mapping, gapOpening=1000, gapExtension=0)
        Tile_read_fq <- preclean_Tile_read_fq[pairwiseAlignment_vec@score<15]
        if( file.exists(file_to_clean)) {     system( paste( "rm", file_to_clean))             }    
        writeFastq( Tile_read_fq,file_to_clean )
        } 
      }
    if(  grepl("ALB1",  Single_file_name)){  
    if( !dir.exists(preclean_tile_file_dir)){  system( paste( "mkdir"  , preclean_tile_file_dir  ))}
    
    tile_to_clean_vec <-  c("54_to_73", "60_to_79", "63_to_82" )
    file_to_clean_vec <- file_to_isotile_vec[as.logical(rowSums(
    sapply( tile_to_clean_vec ,  grepl,file_to_isotile_vec  )))]
    for( file_to_clean in file_to_clean_vec)
      {
      if(! is.element( file_to_clean , list.files(preclean_tile_file_dir))){system(  paste("mv", file_to_clean,  preclean_tile_file_dir) ) }
      ALB1_tile_to_remove <-  ALB1_tile_to_remove_vec[sapply( tile_to_clean_vec ,grepl,  file_to_clean )]
      preclean_Tile_read_fq <-   readFastq(paste0( preclean_tile_file_dir,file_to_clean))
      pairwiseAlignment_vec <-   pairwiseAlignment(pattern=    sread(preclean_Tile_read_fq),
                                                   subject=   ALB1_tile_to_remove , 
                                                   type="overlap",  substitutionMatrix=mapping, gapOpening=1000, gapExtension=0)
      
      Tile_read_fq <- preclean_Tile_read_fq[pairwiseAlignment_vec@score<15]
      if( file.exists(file_to_clean)) {     system( paste( "rm", file_to_clean))             }    
      writeFastq( Tile_read_fq, file_to_clean)       
      } 
    } 
  }

# to check that it worked, rerun this part of the program and look again at the variable in max_vec_array_list

}

##########################################################################
###  To create isotile files from tile files  ####
##############################################################################

if(1)
{
  setwd( Verbose_dir)
  Verbose_dir_vec <- list.files( Verbose_dir)
  
   for( index_file in 1:length(   Verbose_dir_vec ))
  {
    Single_file_name <-   Verbose_dir_vec[index_file ]
    Single_file_dir <- paste(sep="",    Verbose_dir , Single_file_name,"/")
    Tile_file_dir <-  paste( sep="", Single_file_dir, "tile_files/")
    Isotile_file_dir <-  gsub("tile_",  "isotile_",Tile_file_dir)
    
    # to create isotile fastq files from the single 20-bp tile fastq files 
    if(1){
      setwd( Tile_file_dir)
      file_to_isotile_vec <- list.files()
      
      isotile_size <- round(median(sapply(  file_to_isotile_vec , countLines)/4))
      Isotile_file_dir <-  gsub("tile_",  "isotile_",Tile_file_dir)
      if( dir.exists(   Isotile_file_dir )) { system( paste( "rm -r ",    Isotile_file_dir ) ) } 
      system( paste( "mkdir ",    Isotile_file_dir ) ) 
      for( index_tile in 1:length( file_to_isotile_vec ) ){
        file_to_isotile <-     file_to_isotile_vec  [  index_tile]
        Tile_read_fq <-   readFastq(  file_to_isotile )
        Tile_read_size <- length(Tile_read_fq)
        Isotile_read_fq <-   Tile_read_fq [ sample( Tile_read_size, isotile_size, replace = ( Tile_read_size <= isotile_size))]
        Isotile_file_outputname <- paste( sep="",  Isotile_file_dir, gsub( "_tile", "_isotile", file_to_isotile))
        writeFastq(  Isotile_read_fq , Isotile_file_outputname  , "w")
      }
    }
    
    Merge_file_dir <-  paste( sep="", Single_file_dir, "merge_files/" )
   
     # to construct merge 16bp files from single 20bp tile files, both Fastq and positioned fasta
    if(1){
       if( dir.exists(  Merge_file_dir )) { system( paste( "rm -r ",   Merge_file_dir ) ) } 
        system( paste( "mkdir ",   Merge_file_dir ) )  

        setwd(Isotile_file_dir)
        file_to_merge_vec <-  list.files( )  
        Merge_16bp_isotile_file_outputname  <-   paste(sep="", Merge_file_dir,    gsub("_files_","_16bp.fastq.gz" ,    gsub("/","_" , sub('.*Verbose/', '',       Isotile_file_dir ) ) ) )
        Merge_16bp_positioned_isotile_file_outputname   <-   gsub("_16bp.fastq","_16bp_positioned.fasta" ,   Merge_16bp_isotile_file_outputname )
        mode_type <- "w"
        for( index_tile in 1:length( file_to_merge_vec ) )
        {
          file_to_merge <-     file_to_merge_vec  [  index_tile]
          tile_limit_vec <-  as.integer(strsplit( gsub('.fastq.gz', '', sub('.*from_', '', file_to_merge)), "_to_")[[1]])
          file_to_merge_read_fq <-   readFastq(   file_to_merge )
          Intertile_read_fq <- narrow( file_to_merge_read_fq, start=3, end =18)
          
          writeFastq(  Intertile_read_fq ,   Merge_16bp_isotile_file_outputname, mode_type, compress=T)
          Intertile_read_fq@sread  <- DNAStringSet(paste( sep="" , paste( rep( "N",   tile_limit_vec[1]  +1) , collapse="") ,
                                                            sread( Intertile_read_fq)  ,    paste( rep( "N",  147- tile_limit_vec[2] + 2 ), collapse="")))
          writeFasta(  Intertile_read_fq ,   Merge_16bp_positioned_isotile_file_outputname, mode_type, compress=T)
          mode_type <- "a"
        }
        
    }
    
    # to copy merge files in the general directories for Kmer analyses
    if(1){
      setwd(Merge_file_dir)
      merge_file_vec <- list.files()
       file_to_cp_from_name <-  paste(sep="", Single_file_name, "_isotile_16bp.fastq.gz")
      system (paste( "cp", file_to_cp_from_name , isotile_16bp_file_dir) )
      file_to_cp_from_name <-  paste(sep="", Single_file_name, "_isotile_16bp_positioned.fasta.gz")
      system (paste( "cp", file_to_cp_from_name , isotile_16bp_positioned_file_dir) )
    }
  }
  }


#######################################################################
#####  To downsizing isotile files for efficient kmer enrichment ######
#######################################################################

if(1) 
  { 
  setwd( isotile_16bp_file_dir)
  to_downsize_file_vec <- list.files()
  max_read_to_count_TF <- 10^6
  max_read_to_count_BG <- max_read_to_count_TF*5 
  min_read_to_count_TF <-   5*10^4
  min_read_to_count_BG <-   min_read_to_count_TF*5
  BG_file_lab <- grepl("GST", to_downsize_file_vec)|grepl("INPUT", to_downsize_file_vec)
  
  # check
  if(1)
    { 
    countLines_to_downsize_file_vec <- countLines(to_downsize_file_vec)
    sample_fq_read_count_vec <- countLines_to_downsize_file_vec/4
    if( min( sample_fq_read_count_vec)<min_read_to_count_TF)
      { 
      too_small_file_lab <-  sample_fq_read_count_vec<min_read_to_count_TF
      too_small_sample_fq_read_count_vec <- sample_fq_read_count_vec[too_small_file_lab]
      write( paste("Files eliminated because too Small:", names( too_small_sample_fq_read_count_vec ) ,"\n Reads:",  too_small_sample_fq_read_count_vec) ,  stdout())
      to_downsize_file_vec <- to_downsize_file_vec[!too_small_file_lab]
    }
    min(sample_fq_read_count_vec)
  }
  
  max_read_to_count_vec <-ifelse(grepl("INPUT",to_downsize_file_vec) | grepl("GST",to_downsize_file_vec),
                                 max_read_to_count_BG, max_read_to_count_TF )
  names(max_read_to_count_vec) <- to_downsize_file_vec
  
  for( to_downsize_file in  to_downsize_file_vec )
    {
    sample_fq <- readFastq( to_downsize_file )
    sample_length <- length(sample_fq)
    downsize_read_to_count <- min( sample_length, max_read_to_count_vec[to_downsize_file])
    randomization_vec <- sample(  sample_length ,downsize_read_to_count,  replace = FALSE )    
    sample_fq <- sample_fq[ randomization_vec ]
    isotile_16bp_downsized_filename <- paste(sep="", isotile_16bp_downsized_file_dir, to_downsize_file)
    if(file.exists( isotile_16bp_downsized_filename) ) {system( paste( "rm" ,  isotile_16bp_downsized_filename)  )}
    writeFastq(   sample_fq, isotile_16bp_downsized_filename, "w")                       
    remove( sample_fq)
    }
  }

######################################################
#####  Run of Selex package for kmer enrichment ######
######################################################

if(1) 
  { 
  Selex_Analysis_dir <- paste(sep="", Parent_dir,  "Selex_analysis/")
  # Not needed at the beginning. Useful in case you want to rerun the enrichment files from scratch
  # if ( dir.exists(  Selex_Analysis_dir   )) { system( paste( "rm -r ",  Selex_Analysis_dir ) )}
  system( paste( "mkdir ",  Selex_Analysis_dir ))
  
  setwd( isotile_16bp_downsized_file_dir)
  downsized_file_vec <- list.files()
  sequence_size <- 16 
  
  for( template_name in template_name_vec[3])
    {
    template_downsized_file_vec <- downsized_file_vec[grep( template_name, downsized_file_vec )]
    sample_name_vec <-sub('_.*', '', template_downsized_file_vec)
    round_vec <-as.integer(str_match(sample_name_vec, "-Round\\s*(.*?)\\s*-Rep")[,2])
    round_vec[is.na(round_vec)] <- 0
    sample_DF = data.frame(  seqName =  sample_name_vec   , sampleName=  sample_name_vec  ,
                               round=  round_vec,stringsAsFactors=F,   seqFile=   template_downsized_file_vec ,          
                               leftBarcode= '',  rightBarcode= '', leftFlank=  '',  rightFlank= '' )
    Selex_Analysis_Template_dir <- paste(sep="", Selex_Analysis_dir, template_name,"/" ) 
    selex.setwd(path= Selex_Analysis_Template_dir)
    selex.config(workingDir =  Selex_Analysis_Template_dir  , verbose=F, maxThreadNumber=4) ## this can be verbose=FALSE
    sample_handle_list <- list()
    for( index_sample in  1:nrow(sample_DF) )
      {
      sample_vec <- sample_DF[index_sample ,]
      selex.defineSample(seqName =  sample_vec$seqName , sampleName =  sample_vec$sampleName  , 
                         leftBarcode=  sample_vec$leftBarcode,  rightBarcode=  sample_vec$rightBarcode, 
                         leftFlank=   sample_vec$leftFlank,  rightFlank=  sample_vec$rightFlank,
                         round   =  sample_vec$round   , varLength  =  sequence_size,    seqFile=  sample_vec$seqFile )
      sample_handle_list[[index_sample]] <-   selex.sample(seqName=sample_vec$seqName, sampleName=sample_vec$sampleName , round=sample_vec$round)
      }
 
    selex.sampleSummary()
    names( sample_handle_list ) <-  sample_name_vec
      
# To plit the Rd0 sample into testing and training sets ##
    Selex_foreground_file_vec  <- sample_name_vec[!grepl("INPUT", sample_name_vec)]
    Selex_background_file  <- sample_name_vec[grepl("INPUT", sample_name_vec)]
    r0.split = selex.split(sample=sample_handle_list[[    Selex_background_file ]], ratios=c(.5,.5) )

# To find kmax on test dataset and build Markov model on training dataset #####
    k = selex.kmax(sample=r0.split$test)
    mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test, Kmax=k)
    selex.mmSummary()
      
    for( index_comparison in  1:length(Selex_foreground_file_vec) )
      {
      foreground_file_name <-   Selex_foreground_file_vec [index_comparison] 
      Comparison_Kmer_Enrichment_dir <- paste(sep="", Kmer_Enrichment_dir, foreground_file_name,"/")
      system(paste( "mkdir", Comparison_Kmer_Enrichment_dir))
      sample_handle <- sample_handle_list[[ foreground_file_name  ]]
      selex.infogain(sample=   sample_handle ,   markovModel=mm) ## up to 16-mer evaluation
      infogainSummary_df <- selex.infogainSummary()  
      infogainSummary_df <-  infogainSummary_df[grep(foreground_file_name ,infogainSummary_df$Sample),]
      InformationGain_vec <-  infogainSummary_df$InformationGain
      names(InformationGain_vec) <- kmer_length_vec <- infogainSummary_df$K
      save(InformationGain_vec, file= paste(sep="", Comparison_Kmer_Enrichment_dir,"InformationGain_vec.rdata"))
      Comparison_Kmer_Enrichment_subdir <- paste( sep="", Comparison_Kmer_Enrichment_dir, "by_kmer_length/")
      system( paste( " mkdir",  Comparison_Kmer_Enrichment_subdir )) 
      for( kmer_length_index in   1:length(kmer_length_vec) )
        {
        skip_to_next <- FALSE
        kmer_length <- kmer_length_vec[ kmer_length_index ]
        tryCatch(  kmer_table <- selex.counts(sample=sample_handle , k=kmer_length, markovModel=mm), error = function(e) { skip_to_next <<- TRUE})
        if(skip_to_next) { next }  
        if( dim( kmer_table)[1]>0)
          { 
          kmer_table$Increase <- round( kmer_table$ObservedCount/kmer_table$ExpectedCount, digit=2)
          kmer_table <- kmer_table[order(kmer_table$Increase, decreasing = T),]
          kmer_table$ID <- 1:dim( kmer_table )[1]
          kmer_table <- kmer_table[,c("ID", "Kmer" , "ObservedCount", "ExpectedCount", "Increase")]  
          write.table( kmer_table  ,file=paste(sep="",   Comparison_Kmer_Enrichment_subdir,kmer_length,"mer_table.txt"), row.names=F, col.names = T, quote=F) 
          }
        }
      }
    }
  }
  
###########################################################################
######## Process Information Gain and Top Kmer Enrichment (Supp. Table S5) ########
###########################################################################

if(1)
  {
  setwd( Kmer_Enrichment_dir)
  Kmer_Enrichment_subdir_vec <- list.files()
  TF_name_with_7PFs_and_GST_vec <- c("GST", "OCT4","RUNX2",   "PAX7","FOXA1","GATA4","CEBPB","GATA1") 
  TF_name_with_5nonPFs_and_GST_vec <- c( "GST", "ASCL1",  "KLF4",      "PU1"  ,    "ZNF143")
  TF_set_list <- list(  'PFs'=TF_name_with_7PFs_and_GST_vec  ,   'nonPFs'=TF_name_with_5nonPFs_and_GST_vec )
  Kmer_Enrichment_db <- data.frame(t( matrix(c(unlist(strsplit( Kmer_Enrichment_subdir_vec, "-" ))), 5, length(Kmer_Enrichment_subdir_vec),   
  dimnames=list(  c("condition","template","TF", "Round","Replicate") , Kmer_Enrichment_subdir_vec))), stringsAsFactors=F)
  assigned_IG_vec <- c()
  kmer_length_assigned_IG_vec <- c()
  top_kmer_assigned_IG_vec <- c()
  top_kmer_increase_assigned_IG_vec <- c()
  
  for( index_subdir in 1:length(Kmer_Enrichment_subdir_vec))
  {
  Kmer_Enrichment_subdir <- Kmer_Enrichment_subdir_vec[  index_subdir ]
  load( paste( sep="", Kmer_Enrichment_subdir, "/InformationGain_vec.rdata"))

  assigned_kmer_length <- TF_kmer_assigned_length_vec[Kmer_Enrichment_db$TF[index_subdir]]
  assigned_IG_vec[ index_subdir] <- InformationGain_vec[as.character(assigned_kmer_length) ]
  kmer_length_assigned_IG_vec[ index_subdir] <- assigned_kmer_length
  kmer_enrichment_table <- read.table(file= paste( sep="", Kmer_Enrichment_subdir, "/by_kmer_length/", kmer_length_assigned_IG_vec[ index_subdir], "mer_table.txt") , header=T)
  top_kmer_assigned_IG_vec[ index_subdir ] <-  kmer_enrichment_table$Kmer[1]
  top_kmer_increase_assigned_IG_vec[ index_subdir ]  <-  kmer_enrichment_table$Increase[1]
  }
  
  
  Kmer_Enrichment_db$InfoGain <- assigned_IG_vec  
  Kmer_Enrichment_db$kmer_length <- kmer_length_assigned_IG_vec 
  Kmer_Enrichment_db$top_kmer <- top_kmer_assigned_IG_vec 
  Kmer_Enrichment_db$top_kmer_increase <- top_kmer_increase_assigned_IG_vec 
  
  Round4_Kmer_Enrichment_db <-   subset( Kmer_Enrichment_db, Round=="Round4")
  save(  Round4_Kmer_Enrichment_db, file=paste(sep="",  Graphs_for_Kmer_enrichment_dir, "DNA_binding_specificity_table.rdata") )
  
  Round4_Kmer_Enrichment_db_list <- sapply( TF_name_vec, names)
  
  for( TF_name in TF_name_vec)
    {
    subtable <- Round4_Kmer_Enrichment_db_list[[TF_name]] <- subset( Round4_Kmer_Enrichment_db, TF==  TF_name)
    if(TF_name== TF_name_vec[1]){  Round4_Kmer_Enrichment_table <- subtable }
    if(TF_name != TF_name_vec[1]){  Round4_Kmer_Enrichment_table <- rbind(Round4_Kmer_Enrichment_table, subtable ) }
    }

  write.table(  cbind( sample=rownames(Round4_Kmer_Enrichment_table),  Round4_Kmer_Enrichment_table),    row.names=F, col.names = T, quote=F,       sep = "\t",
                file=paste(sep="",  Graphs_dir, "Table_S5.txt")) 
  }

############################################################
#   Plots for Panels 1b, 1e and  S2b, S2c, S2d
#####################################################################################
if(1)
  {
  setwd( Graphs_for_Kmer_enrichment_dir)
  load(  file=paste(sep="",  "DNA_binding_specificity_table.rdata") )

  TF_name_with_7PFs_and_GST_vec <- c("GST", "OCT4","RUNX2",   "PAX7","FOXA1","GATA4","CEBPB","GATA1") 
  TF_name_with_5nonPFs_and_GST_vec <- c( "GST", "ASCL1",  "KLF4",      "PU1"  ,    "ZNF143")
  TF_set_list <- list(  'PFs'=TF_name_with_7PFs_and_GST_vec  ,   'nonPFs'=TF_name_with_5nonPFs_and_GST_vec )
  
  graphics.off()
  for(index_TF_set in 1:2)
  {
    TF_set <-  TF_set_list[[index_TF_set]]
    TF_set_suffix <- names(TF_set_list[index_TF_set])
    PNGlength <- ifelse( index_TF_set==1,   1000,    800)
    TF_set_Round4_Kmer_Enrichment_db <-  Round4_Kmer_Enrichment_db[is.element(  Round4_Kmer_Enrichment_db$TF, TF_set),]
    score_tag_vec <- c(  "top_kmer_increase" , "InfoGain")
    for( score_tag in  score_tag_vec)
      {
      TF_name_matrix <-    matrix( NA, nrow=4, ncol=length(  TF_set ), dimnames=list( template_name_vec,   TF_set ))
      Symbolplot_DNA_VS_NCP_for_TF_set_and_GST <-     paste( sep="",   Graphs_for_Kmer_enrichment_dir,  "Symbolplot_", score_tag,"_DNA_VS_NCP_for_", TF_set_suffix , "_and_GST.png")
      CairoPNG(   PNGlength,500,file= Symbolplot_DNA_VS_NCP_for_TF_set_and_GST   , pointsize = 40/2 ,bg = "white" , res=72*2)
      par(mar=c(4,4,1,1))
      main_text  <-   paste (      score_tag   )
      max_y <- max(TF_set_Round4_Kmer_Enrichment_db[,score_tag]  ,na.rm=T)*1.2
      min_y <- min(TF_set_Round4_Kmer_Enrichment_db[,score_tag]  ,na.rm=T)*0.5
      boxplot(  TF_name_matrix  , las=2 , ylab=score_tag,  col="white",  border="white",
              ylim=c(min_y, max_y) , log="y"  )
       
      label_WANC_list <-  enrichment_ACN_list <-enrichment_W601_list <- list()
      for(condition_ID in c("DNA","NCP"))
        {
        condition_TF_set_Round4_Kmer_Enrichment_db <-  subset( TF_set_Round4_Kmer_Enrichment_db, condition==condition_ID )
        label_WANC_vec <- enrichment_W601_vec <- enrichment_ACN_vec <- c()
        for(Rep_ID in c("Rep1","Rep2"))
          {
          Rep_condition_TF_set_Round4_Kmer_Enrichment_db <-  subset( condition_TF_set_Round4_Kmer_Enrichment_db, Replicate==Rep_ID)
          casted_matrix <- acast(  Rep_condition_TF_set_Round4_Kmer_Enrichment_db[, c("template","TF",score_tag )] , template~TF)
          condition_ID_TF_name_matrix <- TF_name_matrix
          condition_ID_TF_name_matrix[,colnames(casted_matrix)] <- casted_matrix
  
          if( condition_ID =="DNA"){ color_ID="#f4b553"; plot_space=-1}
          if( condition_ID =="NCP"){ color_ID="#7F3F98"; plot_space=0}
  
          x_coord <- 1:length(TF_set)+plot_space*0.2
          points(x_coord,  condition_ID_TF_name_matrix ["ALB1", ],  pch=15 , col= color_ID )
          points( x_coord, condition_ID_TF_name_matrix ["CX3", ], pch=16,col= color_ID)
          points( x_coord ,  condition_ID_TF_name_matrix ["NRCAM", ] , pch=17, col= color_ID )
          points( x_coord,  condition_ID_TF_name_matrix ["W601", ], pch=18, col= color_ID )
          enrichment_W601_vec <- c(enrichment_W601_vec, condition_ID_TF_name_matrix["W601",setdiff(TF_set, "GST")  ] )
          enrichment_ACN_vec <- c(enrichment_ACN_vec, colMeans(condition_ID_TF_name_matrix[c("ALB1", "CX3", "NRCAM"),setdiff(TF_set, "GST")], na.rm=T))
          label_WANC_vec <- c( label_WANC_vec, paste0(condition_ID,"_",Rep_ID,"_", setdiff(TF_set, "GST") ))
          }
        enrichment_ACN_list[[condition_ID]] <-    enrichment_ACN_vec
        enrichment_W601_list[[condition_ID]] <-    enrichment_W601_vec
        label_WANC_list[[condition_ID]] <- label_WANC_vec 
        }
    
      treshold_line <-  2*max(subset( TF_set_Round4_Kmer_Enrichment_db , TF=="GST" )[, score_tag])

      if(!grepl("top_kmer",    score_tag  ) ){  abline( h=   treshold_line  , col="black" , lty=5) }
      dev.off()

      # Plot for panel 1b 
      if(  grepl( "top_kmer_increase_DNA_VS_NCP_for_PFs_and_GST", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST))
      {system( paste0( "cp ", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST , " ",Graphs_dir,"Fig_1b.png" ))}
  
      # Plot for panel S2b 
      if(  grepl( "top_kmer_increase_DNA_VS_NCP_for_nonPFs_and_GST", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST))
      {system( paste0( "cp ", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST , " ",Graphs_dir,"Fig_S2b.png" ))}
    
      # Plot for panel S2c 
      if(  grepl( "InfoGain_DNA_VS_NCP_for_PFs_and_GST", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST))
      {system( paste0( "cp ", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST , " ",Graphs_dir,"Fig_S2c.png" ))}
    
      # Plot for panel S2d 
      if(  grepl( "InfoGain_DNA_VS_NCP_for_nonPFs_and_GST", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST))
      {system( paste0( "cp ", Symbolplot_DNA_VS_NCP_for_TF_set_and_GST , " ",Graphs_dir,"Fig_S2d.png" ))}
    
      # Fig 1e
      if( score_tag=="InfoGain" & TF_set_suffix =="PFs") 
      {
        Round4_score_connected_dots_W601_VS_ACN_with_PFs_and_GST <-   paste( sep="",   Graphs_for_Kmer_enrichment_dir,  
                                    "Round4_", score_tag,"_connected_dots_W601_VS_ACN_with_",  TF_set_suffix , "_and_GST.png")
        CairoPNG(  750,750,file=      Round4_score_connected_dots_W601_VS_ACN_with_PFs_and_GST   , pointsize = 40/2 ,bg = "white" , res=72*2)
      
        max_value <-     max(unlist(enrichment_W601_list, enrichment_ACN_list), na.rm=T)
        min_value <-     min(unlist(enrichment_W601_list, enrichment_ACN_list), na.rm=T)
      
        plot(c(),xlim=c(0,5), ylim=c(0.2* min_value,1.2*max_value), ylab=score_tag, log="y" )    
        color_ID <- "#f4b553"
        segments( 1,enrichment_W601_list$DNA, 2,enrichment_ACN_list$DNA , color_ID )
        points( rep(1,length(enrichment_W601_list$DNA)), enrichment_W601_list$DNA, pch=18, col= color_ID  )
        points( rep(2,length(enrichment_ACN_list$DNA)), enrichment_ACN_list$DNA, pch=8, col= color_ID  )
        T_test_DNA_P_value <- t.test(enrichment_W601_list$DNA, enrichment_ACN_list$DNA, paired = TRUE, alternative = "less") # alternative = "two.sided"
  
        color_ID <- "#7F3F98"
        segments( 3,enrichment_W601_list$NCP, 4,enrichment_ACN_list$NCP, col=color_ID )
        points( rep(3,length(enrichment_W601_list$NCP)), enrichment_W601_list$NCP, pch=18, col= color_ID )
        points( rep(4,length(enrichment_ACN_list$NCP)), enrichment_ACN_list$NCP, pch=8, col= color_ID )
        T_test_NCP_P_value <- t.test(enrichment_W601_list$NCP, enrichment_ACN_list$NCP, paired = TRUE, alternative = "less") # alternative = "two.sided"
        title( paste( "P-v: DNA =", round( T_test_DNA_P_value$p.value, digits=3),
                    ", NCP =" , round( T_test_NCP_P_value$p.value, digits=3) ) )
        dev.off()

        system( paste0( "cp ", Round4_score_connected_dots_W601_VS_ACN_with_PFs_and_GST , " ",Graphs_dir,"Fig_1e.png" ))
    
        Fig1e_table_DNA <- cbind( W601 = enrichment_W601_list$DNA, ALBN1_CX3CR1_NRCAM=enrichment_ACN_list$DNA)
        rownames(Fig1e_table_DNA) <- label_WANC_list$DNA
        Fig1e_table_DNA <- Fig1e_table_DNA[!is.na(Fig1e_table_DNA[,1]),]
        Fig1e_table_NCP <- cbind( W601 = enrichment_W601_list$NCP,  ALBN1_CX3CR1_NRCAM=enrichment_ACN_list$NCP)
        rownames(Fig1e_table_NCP) <- label_WANC_list$NCP
        Fig1e_table_NCP <- Fig1e_table_NCP[!is.na(Fig1e_table_NCP[,1]),]
        Fig1e_table <- rbind(    Fig1e_table_DNA ,    Fig1e_table_NCP)
        write.table( cbind( samples=rownames(Fig1e_table), Fig1e_table), row.names=F, col.names = T, quote=F,       sep = "\t",
                  file=paste(sep="",  Graphs_dir, "Fig_1e_datasets.txt")) 
        }      
      }
    }
  }

###########################################################################
######## Main Fig 1c-d  (Correlation Between Replicates) ########
###########################################################################

if(1) 
  {
  graphics.off() 
  setwd(Kmer_Enrichment_dir)
  Correlation_between_replicates_scatterplot_dir <-   paste( sep="",  Graphs_for_Kmer_enrichment_dir , "Correlation_between_replicates_scatterplots/")
  if (! dir.exists(   Correlation_between_replicates_scatterplot_dir  )) { system( paste( "mkdir ",    Correlation_between_replicates_scatterplot_dir ) )}
  
  PF_set_vec <- c("OCT4","RUNX2", "PAX7","FOXA1","GATA4","CEBPB","GATA1")
  Kmer_Enrichment_subdir_vec <- list.files(  )
  Round4_Kmer_Enrichment_subdir_vec <- Kmer_Enrichment_subdir_vec[ grepl("Round4", Kmer_Enrichment_subdir_vec)]
  Round4_Kmer_Enrichment_subdir_vec <- Round4_Kmer_Enrichment_subdir_vec[ grepl("NCP", Round4_Kmer_Enrichment_subdir_vec)]
  PF_Round4_Kmer_Enrichment_subdir_vec <- Round4_Kmer_Enrichment_subdir_vec[as.logical(rowSums( sapply( PF_set_vec, grepl, Round4_Kmer_Enrichment_subdir_vec)) ) ]
  correlation_var_vec <- c(corr_method ,"first","second")
  Kmer_Enrichment_correlation_df = data.frame(matrix(vector(), 0, length(correlation_var_vec),   dimnames=list(c(), correlation_var_vec)), stringsAsFactors=F)
  
  for( first_file_subdir in PF_Round4_Kmer_Enrichment_subdir_vec)
    {
    first_file_variable_vec <-    strsplit(first_file_subdir,"-")[[1]]
    PF_kmer_length <- TF_kmer_assigned_length_vec[first_file_variable_vec[3]]
    first_Rep_ID <- first_file_variable_vec[5]
    second_Rep_ID <-ifelse( first_Rep_ID=="Rep1", "Rep2", "Rep1")
    second_file_subdir <- gsub( first_Rep_ID,  second_Rep_ID, first_file_subdir)

    if(dir.exists(second_file_subdir) )
      {
      first_file_graph_dir <- paste(sep="", Correlation_between_replicates_scatterplot_dir, first_file_subdir, "/")
      if( !dir.exists(first_file_graph_dir)) { system( paste("mkdir",first_file_graph_dir) ) }
      first_file_address <-    paste(sep="",    first_file_subdir, "/by_kmer_length/",  PF_kmer_length ,"mer_table.txt")
      first_file_kmer_enrichment_table <-   read.table(     first_file_address  , header=T)
      second_file_address <-    paste(sep="",    second_file_subdir, "/by_kmer_length/",  PF_kmer_length ,"mer_table.txt")
      second_file_kmer_enrichment_table <-   read.table(         second_file_address  , header=T)
      kmer_enrichment_match_vec <- match(second_file_kmer_enrichment_table$Kmer, first_file_kmer_enrichment_table$Kmer )
      corr_table <-   cbind( first_file_kmer_enrichment_table [kmer_enrichment_match_vec,], second_file_kmer_enrichment_table )
      corr_table  <-   corr_table [which(!is.na(kmer_enrichment_match_vec)),] 
      rownames( corr_table ) <- corr_table$Kmer
      corr_table  <-   corr_table [,c(5,10)] 
      colnames( corr_table ) <- c( first_file_subdir, second_file_subdir  )
      corr_coeff_value <-  cor( corr_table, method = corr_method)[2,1]  
      Kmer_Enrichment_correlation_df[dim(Kmer_Enrichment_correlation_df)[1]+1, ] <-   c( corr_coeff_value , first_file_subdir, second_file_subdir )
      Between_replicates_Kmer_enrichment_plot <-  
      paste( sep="", first_file_graph_dir ,first_file_subdir, "_VS_", second_Rep_ID, ".png")
      CairoPNG(  1000,1000,file= Between_replicates_Kmer_enrichment_plot  , pointsize = 40/2 ,bg = "white" , res=72*2)
      plot(  corr_table, xlab=first_file_subdir, ylab=second_file_subdir, 
             pch=16, main= paste( corr_method , "=" ,   round( corr_coeff_value, digits=3) ))
      dev.off()
      
      write.table(  cbind( kmer=rownames(corr_table ),corr_table )  ,  row.names=F, quote=F,       sep = "\t",
                    file=gsub( ".png" , ".txt", Between_replicates_Kmer_enrichment_plot )) 
      
      if( grepl( "NCP-W601-PAX7-Round4-Rep1_VS_Rep2", Between_replicates_Kmer_enrichment_plot))
        {
        system( paste0( "cp " , Between_replicates_Kmer_enrichment_plot, " ", Graphs_dir,"Fig_1c.png" ))
        system( paste0( "cp " ,gsub( ".png" , ".txt", Between_replicates_Kmer_enrichment_plot ), " ", Graphs_dir,"Fig_1c_datasets.txt" ))
        }
      }
    }
  
  Kmer_Enrichment_correlation_df_filename <-  paste( sep="", Graphs_for_Kmer_enrichment_dir, "Kmer_Enrichment_correlation_between_replicates_df.rdata")
  save( Kmer_Enrichment_correlation_df,  file= Kmer_Enrichment_correlation_df_filename )
  write.table( Kmer_Enrichment_correlation_df [ c(T,F) , ],  file= "../Graphs/Fig1d_datasets.txt",
               row.names=F, col.names = c( "Corr_coeff", "first_sample", "second_sample"  ), quote=F,       sep = "\t",)
  
  for( condition in c("NCP"))
  {
    condition_correlation_df <- Kmer_Enrichment_correlation_df[grepl(condition , Kmer_Enrichment_correlation_df$first) , ]
    first_W601_lab <- grepl("W601" , condition_correlation_df$first)
    correlation_W601_df <- condition_correlation_df[first_W601_lab , ]
    correlation_W601_df <- correlation_W601_df[ !duplicated(correlation_W601_df[,corr_method]), ]
    correlation_nonW601_df <- condition_correlation_df[!first_W601_lab , ]
    correlation_nonW601_df <- correlation_nonW601_df[ !duplicated(correlation_nonW601_df[,corr_method]), ]
    corr_coeff_W601 <- as.double( correlation_W601_df[,corr_method])
    corr_coeff_nonW601 <- as.double( correlation_nonW601_df[,corr_method])
    corr_coeff_list <- list ( W601=corr_coeff_W601,  nonW601=corr_coeff_nonW601) 
    Between_replicates_Kmer_enrichment_corr_coeff_boxplot <-      paste( sep="", Graphs_for_Kmer_enrichment_dir ,  "Between_replicates_" ,  condition,"_Kmer_enrichment_corr_coeff_boxplot.png")

    CairoPNG(  400,600,file=      Between_replicates_Kmer_enrichment_corr_coeff_boxplot   , pointsize = 40/2 ,bg = "white" , res=72*2)
    par(mar = c(1,1, 1, 1) ) 
    max_value <-     max(unlist( corr_coeff_list), na.rm=T)
    min_value <-     min(unlist( corr_coeff_list), na.rm=T)
    boxplot( corr_coeff_list , ylim=c( min_value, 1), outline=FALSE)
    points(  rep(1,length(corr_coeff_list[[1]])), corr_coeff_list[[1]] , pch=18)
    points(  rep(2,length(corr_coeff_list[[2]])), corr_coeff_list[[2]] , pch=8)
    Pvalue <- round( wilcox.test( corr_coeff_list$W601, corr_coeff_list$nonW601 ,  alternative = c("less"))$p.value, 3)
    title( paste( condition, "Pvalue", Pvalue ))
    dev.off()
    system( paste0( "cp ",  Between_replicates_Kmer_enrichment_corr_coeff_boxplot ,
                    " ", Graphs_dir,"Fig_1d.png" ))
    }
  }

###########################################################################
######## Main Fig 1f-g  (Correlation Between Templates) ########
###########################################################################

if(1) 
  {
  graphics.off() 
  setwd(Kmer_Enrichment_dir)
  Correlation_between_templates_scatterplot_dir <-   paste( sep="",  Graphs_for_Kmer_enrichment_dir , "Correlation_between_templates_scatterplots/")
  if (! dir.exists(   Correlation_between_templates_scatterplot_dir  )) { system( paste( "mkdir ",    Correlation_between_templates_scatterplot_dir ) )}
  PF_set_vec <- c("OCT4","RUNX2",   "PAX7","FOXA1","GATA4","CEBPB","GATA1")
  Kmer_Enrichment_subdir_vec <- list.files(  )
  Round4_Kmer_Enrichment_subdir_vec <- Kmer_Enrichment_subdir_vec[ grepl("Round4", Kmer_Enrichment_subdir_vec)]
  Round4_Kmer_Enrichment_subdir_vec <- Round4_Kmer_Enrichment_subdir_vec[ grepl("NCP", Round4_Kmer_Enrichment_subdir_vec)]
  PF_Round4_Kmer_Enrichment_subdir_vec <- Round4_Kmer_Enrichment_subdir_vec[as.logical(rowSums( sapply( PF_set_vec, grepl, Round4_Kmer_Enrichment_subdir_vec)) ) ]
  correlation_var_vec <- c(  corr_method,"first","second")
  Kmer_Enrichment_correlation_df = data.frame(matrix(vector(), 0, length(correlation_var_vec),   dimnames=list(c(), correlation_var_vec)), stringsAsFactors=F)
  
  for( first_file_subdir in PF_Round4_Kmer_Enrichment_subdir_vec)
  {
    first_file_graph_dir <- paste(sep="", Correlation_between_templates_scatterplot_dir, first_file_subdir, "/")
    if( !dir.exists(first_file_graph_dir)) { system( paste("mkdir",first_file_graph_dir) ) }
    first_file_variable_vec <-    strsplit(first_file_subdir,"-")[[1]]
    PF_kmer_length <- TF_kmer_assigned_length_vec[first_file_variable_vec[3]]
    first_template <- first_file_variable_vec[2]
    
    first_template_file_address <-    paste(sep="",    first_file_subdir, "/by_kmer_length/",  PF_kmer_length ,"mer_table.txt")
    first_template_kmer_enrichment_table <-   read.table(     first_template_file_address  , header=T)
    
    second_template_vec <-   setdiff( template_name_vec, first_template)
    for( second_template in second_template_vec )
      {
      second_file_subdir <- gsub( first_template, second_template,   first_file_subdir)
      if(is.element( second_file_subdir, PF_Round4_Kmer_Enrichment_subdir_vec) )
        {
        second_template_file_address <-    paste(sep="",    second_file_subdir, "/by_kmer_length/",  PF_kmer_length ,"mer_table.txt")
        second_template_kmer_enrichment_table <-   read.table(         second_template_file_address  , header=T)
        kmer_enrichment_match_vec <- match(second_template_kmer_enrichment_table$Kmer, first_template_kmer_enrichment_table$Kmer )
        corr_table <-   cbind( first_template_kmer_enrichment_table [kmer_enrichment_match_vec,], second_template_kmer_enrichment_table )
        corr_table  <-   corr_table [which(!is.na(kmer_enrichment_match_vec)),] 
        rownames( corr_table ) <- corr_table$Kmer
        corr_table  <-   corr_table [,c(5,10)] 
        colnames( corr_table ) <- c( first_file_subdir, second_file_subdir  )
        corr_coeff_value <-  cor( corr_table, method = corr_method)[2,1]  
        Kmer_Enrichment_correlation_df[dim(Kmer_Enrichment_correlation_df)[1]+1, ] <-   c( corr_coeff_value , first_file_subdir, second_file_subdir )
        Between_templates_Kmer_enrichment_plot <-  
        paste( sep="", first_file_graph_dir ,first_file_subdir, "_with_", second_template, ".png")
        CairoPNG(  1000,1000,file= Between_templates_Kmer_enrichment_plot  , pointsize = 40/2 ,bg = "white" , res=72*2)
        plot(  corr_table, xlab=first_file_subdir, ylab=second_file_subdir, 
               pch=16, main= paste( corr_method , "=" ,   round( corr_coeff_value, digits=3) ))
        dev.off()
        
        write.table(  cbind( kmer=rownames(corr_table ),corr_table )  ,  row.names=T, col.names = T, quote=F,       sep = "\t",
                      file=gsub( ".png" , ".txt", Between_templates_Kmer_enrichment_plot )) 
        if( grepl( "NCP-ALB1-PAX7-Round4-Rep2_with_W601", Between_templates_Kmer_enrichment_plot))
          {
          system( paste0( "cp " , Between_templates_Kmer_enrichment_plot , " ", Graphs_dir,"Fig_1f.png" ))
          system( paste0( "cp " , gsub( ".png" , ".txt", Between_templates_Kmer_enrichment_plot ), " ", Graphs_dir,"Fig_1f_datasets.txt" ))
          }
        }
      }      
    }
  
  Kmer_Enrichment_correlation_df_filename <-  paste( sep="", Graphs_for_Kmer_enrichment_dir, "Kmer_Enrichment_correlation_between_templates_df.rdata")
  save( Kmer_Enrichment_correlation_df,  file= Kmer_Enrichment_correlation_df_filename )
  write.table( Kmer_Enrichment_correlation_df ,  file= "../Graphs/Fig1g_datasets.txt",
               row.names=F, col.names = c( "Corr_coeff", "first_sample", "second_sample"  ), quote=F,       sep = "\t",)
  
  for( condition in c( "NCP"))
  {
    condition_correlation_df <- Kmer_Enrichment_correlation_df[grepl(condition , Kmer_Enrichment_correlation_df$first) , ]
    first_W601_lab <- grepl("W601" , condition_correlation_df$first)
    correlation_W601_df <- condition_correlation_df[first_W601_lab , ]
    second_W601_lab <- grepl("W601" , condition_correlation_df$second)
    correlation_nonW601_df <- condition_correlation_df[!first_W601_lab & !second_W601_lab , ]
    correlation_nonW601_df <- correlation_nonW601_df[ !duplicated(correlation_nonW601_df[,corr_method]), ]
    corr_coeff_W601 <- as.double( correlation_W601_df[,corr_method])
    corr_coeff_nonW601 <- as.double( correlation_nonW601_df[,corr_method])
    corr_coeff_list <- list ( W601=corr_coeff_W601,  nonW601=corr_coeff_nonW601) 
    Between_templates_Kmer_enrichment_corr_coeff_boxplot <-      paste( sep="", Graphs_for_Kmer_enrichment_dir ,  "Between_templates_" ,  condition,"_Kmer_enrichment_corr_coeff_boxplot.png")
    CairoPNG(  400,600,file=      Between_templates_Kmer_enrichment_corr_coeff_boxplot   , pointsize = 40/2 ,bg = "white" , res=72*2)
    par(mar = c(1,1, 1, 1) ) 
    max_value <-     max(unlist( corr_coeff_list), na.rm=T)
    min_value <-     min(unlist( corr_coeff_list), na.rm=T)
    
    boxplot( corr_coeff_list , ylim=c( min_value, 1), outline=FALSE)
    
    points(  rep(1,length(corr_coeff_list[[1]])), corr_coeff_list[[1]] , pch=18)
    points(  rep(2,length(corr_coeff_list[[2]])), corr_coeff_list[[2]] , pch=8)
    Pvalue <- round( wilcox.test( corr_coeff_list$W601, corr_coeff_list$nonW601 ,  alternative = c("less"))$p.value, 3)
    title( paste( condition, "Pvalue", Pvalue ))
    dev.off()
    
    system( paste0( "cp ",       Between_templates_Kmer_enrichment_corr_coeff_boxplot ,
                    " ", Graphs_dir,"Fig_1g.png" ))
    }
  }

#######################################################################################
## Main Fig 2a and Ext. Data Fig. S3a, 3b and 3c (Correlation DNA VS NCP)  ####
#######################################################################################

if(1) 
  {
  graphics.off() 
  setwd(Kmer_Enrichment_dir)
  Correlation_DNA_VS_NCP_scatterplot_dir <-   paste( sep="",  Graphs_for_Kmer_enrichment_dir , "Correlation_DNA_VS_NCP_scatterplots/")
  if (! dir.exists(   Correlation_DNA_VS_NCP_scatterplot_dir  )) { system( paste( "mkdir ",    Correlation_DNA_VS_NCP_scatterplot_dir ) )}
  
  PF_set_vec <- c("OCT4","RUNX2", "PAX7","FOXA1","GATA4","CEBPB","GATA1")
  Kmer_Enrichment_subdir_vec <- list.files(  )
  NCP_Kmer_Enrichment_subdir_vec <- Kmer_Enrichment_subdir_vec[ grepl("NCP", Kmer_Enrichment_subdir_vec)]
  PF_NCP_Kmer_Enrichment_subdir_vec <- NCP_Kmer_Enrichment_subdir_vec[as.logical(rowSums( sapply( PF_set_vec, grepl, NCP_Kmer_Enrichment_subdir_vec)) ) ]
  DNA_Kmer_Enrichment_subdir_vec <- Kmer_Enrichment_subdir_vec[ grepl("DNA", Kmer_Enrichment_subdir_vec)]
  
  correlation_var_vec <- c(corr_method,"first","second")
  Kmer_Enrichment_correlation_df = data.frame(matrix(vector(), 0, length(correlation_var_vec),   dimnames=list(c(), correlation_var_vec)), stringsAsFactors=F)
  
  NCP_VS_DNA_sample_vec <- c( "NCP-W601-CEBPB-Round4-Rep1_VS_DNA-W601-CEBPB-Round1-Rep1",
                              "NCP-W601-FOXA1-Round4-Rep1_VS_DNA-W601-FOXA1-Round2-Rep2",
                              "NCP-W601-GATA1-Round4-Rep2_VS_DNA-W601-GATA1-Round2-Rep1",
                              "NCP-W601-GATA4-Round4-Rep2_VS_DNA-W601-GATA4-Round1-Rep1",
                              "NCP-W601-OCT4-Round4-Rep2_VS_DNA-W601-OCT4-Round2-Rep1",
                              "NCP-W601-PAX7-Round4-Rep2_VS_DNA-W601-PAX7-Round1-Rep1",
                              "NCP-W601-RUNX2-Round4-Rep1_VS_DNA-W601-RUNX2-Round2-Rep1",
                              "NCP-NRCAM-PAX7-Round4-Rep1_VS_DNA-NRCAM-PAX7-Round1-Rep1",
                              "NCP-CX3-RUNX2-Round4-Rep2_VS_DNA-CX3-RUNX2-Round1-Rep1")
  names( NCP_VS_DNA_sample_vec ) <- c("Fig_2a", "Fig_S3a_FOXA1","Fig_S3a_GATA1","Fig_S3a_GATA4","Fig_S3a_OCT4" ,"Fig_S3a_PAX7", "Fig_S3a_RUNX2",
                                       "Fig_S3b", "Fig_S3c" )
  
  for( first_file_subdir in PF_NCP_Kmer_Enrichment_subdir_vec)
    {
    first_file_variable_vec <-    strsplit(first_file_subdir,"-")[[1]]
    PF_kmer_length <- TF_kmer_assigned_length_vec[first_file_variable_vec[3]]
    tag_to_grep <-str_match(first_file_subdir, "NCP-\\s*(.*?)\\s*-Round")[,2]
    second_file_subdir_vec <- DNA_Kmer_Enrichment_subdir_vec[ grepl(  tag_to_grep, DNA_Kmer_Enrichment_subdir_vec) ]
      
    if(length(second_file_subdir_vec)>0 )
      {
      first_file_graph_dir <- paste(sep="", Correlation_DNA_VS_NCP_scatterplot_dir
                                    , first_file_subdir, "/")
      if( !dir.exists(first_file_graph_dir)) { system( paste("mkdir",first_file_graph_dir) ) }
      first_file_address <-    paste(sep="",    first_file_subdir, "/by_kmer_length/",  PF_kmer_length ,"mer_table.txt")
      first_file_kmer_enrichment_table <-   read.table(     first_file_address  , header=T)
      
      for(second_file_subdir in second_file_subdir_vec )
        {
        second_file_address <-    paste(sep="",    second_file_subdir, "/by_kmer_length/",  PF_kmer_length ,"mer_table.txt")
        second_file_kmer_enrichment_table <-   read.table(         second_file_address  , header=T)
        kmer_enrichment_match_vec <- match(second_file_kmer_enrichment_table$Kmer, first_file_kmer_enrichment_table$Kmer )
        corr_table <-   cbind( first_file_kmer_enrichment_table [kmer_enrichment_match_vec,], second_file_kmer_enrichment_table )
        corr_table  <-   corr_table [which(!is.na(kmer_enrichment_match_vec)),] 
        rownames( corr_table ) <- corr_table$Kmer
        corr_table  <-   corr_table [,c(5,10)] 
        colnames( corr_table ) <- c( first_file_subdir, second_file_subdir  )
        corr_coeff_value <-  cor( corr_table, method = corr_method)[2,1]  
        Kmer_Enrichment_correlation_df[dim(Kmer_Enrichment_correlation_df)[1]+1, ] <-   c( corr_coeff_value , first_file_subdir, second_file_subdir )
      
        NCP_VS_DNA_Kmer_enrichment_plot <-  
        paste( sep="", first_file_graph_dir ,first_file_subdir, "_VS_", second_file_subdir , ".png")
        CairoPNG(  1000,1000,file= NCP_VS_DNA_Kmer_enrichment_plot  , pointsize = 40/2 ,bg = "white" , res=72*2)
        plot(  corr_table, xlab=first_file_subdir, ylab=second_file_subdir, 
               pch=16, main= paste( corr_method , "=" ,   round( corr_coeff_value, digits=3) ))
        dev.off()
      
        DNA_VS_NCP_Kmer_enrichment_plot <-  paste( sep="", first_file_graph_dir ,second_file_subdir, "_VS_", first_file_subdir, ".png")
        CairoPNG(  1000,1000,file= DNA_VS_NCP_Kmer_enrichment_plot  , pointsize = 40/2 ,bg = "white" , res=72*2)
        plot(  corr_table[,2:1], xlab=second_file_subdir, ylab=first_file_subdir, 
               pch=16, main= paste( corr_method , "=" ,   round( corr_coeff_value, digits=3) ))
        dev.off()

        write.table(  cbind( kmer=rownames(corr_table ),corr_table )  ,  row.names=F, col.names = T, quote=F,       sep = "\t",
                    file=gsub( ".png" , ".txt", NCP_VS_DNA_Kmer_enrichment_plot )) 
        for( Fig_name in names( NCP_VS_DNA_sample_vec ))
          {
          if( grepl( NCP_VS_DNA_sample_vec[Fig_name], NCP_VS_DNA_Kmer_enrichment_plot ))
            {
            system( paste0( "cp " , DNA_VS_NCP_Kmer_enrichment_plot , " ", Graphs_dir, Fig_name,".png" ))
            system( paste0( "cp " , gsub( ".png" , ".txt", NCP_VS_DNA_Kmer_enrichment_plot ), " ", 
                        Graphs_dir,Fig_name, "_datasets.txt" ))
            }
          }
        }
      }
    }
  }

######################################################################################
############# Main Fig. 2b, Fig. 3a and Ext. Data Fig. S7   #################################################
#######################################################################################

if(1)
  {
  # - Function to call kmer profiles in files
  to_get_kmer_match_frequency <- function(  single_kmer, many_sequences){
    kmer_length <- nchar(single_kmer)
    end_count <- 147-kmer_length+1
    merged_single_sequence_from_all <- paste0( many_sequences, collapse="")
    matches_kmer <-   matchPattern( single_kmer ,  merged_single_sequence_from_all   ) 
    corrected_start_sites <- start(matches_kmer) %% nchar(many_sequences)[1]
    start_site_predistribution <- table(corrected_start_sites)
    start_site_distribution <- setNames(rep(0,times= end_count), 1:end_count)
    start_site_distribution[names(start_site_predistribution)] <- start_site_predistribution
    output_df <- start_site_distribution
    return( output_df    )
  }
  
  setwd(isotile_16bp_positioned_file_dir)
  graph_subir_suffix <- "isotile_16bp_positioned"

  Graphs_subdir <- paste( sep="",    Graphs_for_Kmer_profiles_dir,graph_subir_suffix ,"/" )
  if( !dir.exists(Graphs_subdir)){ system( paste( "mkdir", Graphs_subdir ))}
  dir_file_suffix <- paste( sep="", "_", graph_subir_suffix ,".fasta.gz" )
  sequence_size <- 16 
  
  ### most enriched kmers to scan (both FW and RV orientations)
  CEBPB_kmer_vec <-     c( "TTGCGCAA","TTGCGTAA", "TTACGCAA","TTGCGAAA","TTTCGCAA")
  FOXA1_kmer_vec <- c( "TGTTTACTT",  "AAGTAAACA" )
  GATA_kmer_vec <- c(  "GATAAGAT", "ATCTTATC",    "ATCTTGAT", "ATCAAGAT"  )
  OCT4_kmer_vec <- c(  "AATTAGCATA",  "TATGCTAATT",  "AATTTGCATA",  "TATGCAAATT")# "AATTAACATA","TATGTTAATT")
  PAX7_kmer_vec <-  c( "AATCGATT",  "AATTGATT", "AATCAATT" )
  RUNX2_kmer_vec <- c( "TTGCGGTT",  "AACCGCAA" )
  
  kmer_positional_perference_graph_list <-    list(
      Fig_S3d_DNA  =list( 'CEBPB-Round1-Rep1'=  CEBPB_kmer_vec  ,#problematic on CX3
                           'PAX7-Round1-Rep1' = PAX7_kmer_vec,      
                           'RUNX2-Round1-Rep1' = RUNX2_kmer_vec),
     
      Fig_2b_NCP=list( 'CEBPB-Round4-Rep1'=  CEBPB_kmer_vec ),
         
      Fig_3a_NCP =list( 'CEBPB-Round4-Rep1'=  CEBPB_kmer_vec  ,
                          'PAX7-Round4-Rep2' = PAX7_kmer_vec,      
                          'RUNX2-Round4-Rep1' = RUNX2_kmer_vec),
     
      Fig_S7a_NCP=list( 'OCT4-Round4-Rep2'=  OCT4_kmer_vec),
       
      Fig_S7a_DNA=list( 'OCT4-Round2-Rep1'=  OCT4_kmer_vec),
        
      Fig_S7b_NCP=list( 'FOXA1-Round4-Rep1'=  FOXA1_kmer_vec),
      
      Fig_S7b_DNA=list( 'FOXA1-Round2-Rep2'=  FOXA1_kmer_vec),
      
      Fig_S7c_NCP=list( 'GATA1-Round4-Rep2'=  GATA_kmer_vec),
        
      Fig_S7c_DNA=list( 'GATA1-Round2-Rep1'=  GATA_kmer_vec),
      
      Fig_S7d_NCP=list( 'GATA4-Round4-Rep2'=  GATA_kmer_vec),

      Fig_S7d_DNA=list( 'GATA4-Round1-Rep1'=  GATA_kmer_vec))

  max_frequency <- 0.05

  for( index_graph in  1:length(kmer_positional_perference_graph_list) ) 
    { 
    prefix_graph_name <- names(kmer_positional_perference_graph_list)[index_graph]
    prefix_graph_dir <- paste( sep="", Graphs_subdir, prefix_graph_name,"/" )
    if( !dir.exists(prefix_graph_dir)){ system( paste( "mkdir", prefix_graph_dir))}
    kmer_positional_perference_file_list <- kmer_positional_perference_graph_list [[prefix_graph_name]]
    all_kmers <-     unlist(  kmer_positional_perference_file_list)    
    condition <- ifelse( grepl("DNA", prefix_graph_name),"DNA", "NCP" )
    
    for( template_name in template_name_vec)
      {
      prefix_file <- paste(sep="-", condition, template_name )
      kmer_length <- nchar(all_kmers[1]) 
      end_window <- 147- kmer_length+1
      all_kmer_distribution_matrix <-  array(data=NA, dim=c(length( all_kmers), end_window), dimnames=list( all_kmers, 1:end_window ) )
      for( index_file in 1:length(  kmer_positional_perference_file_list) )
        { 
        suffix_file_name <- names(  kmer_positional_perference_file_list)[index_file]
        file_kmer_to_scan_vec <-   kmer_positional_perference_file_list [[ suffix_file_name]]
        file_name <- paste(sep="-", prefix_file, suffix_file_name)
        Fasta_file_address <-   paste( sep="", file_name, dir_file_suffix)
        if( file.exists(Fasta_file_address))
          {
          Fasta_file <- readFasta(    Fasta_file_address  )
          if( length(Fasta_file)>10^6)
            { 
            resampling_lab <- sample( length(Fasta_file), 10^6 )
            Fasta_file <- Fasta_file[ resampling_lab ] 
            }
          TF_kmer_match_info_list_df <-    sapply( file_kmer_to_scan_vec ,  to_get_kmer_match_frequency , sread( Fasta_file ))
          start_point <-  tile_sequence_start_vec_list[[template_name]][1] + 2
          end_point  <-  tail( tile_sequence_start_vec_list[[template_name]], n=1)+sequence_size+2-kmer_length
          all_kmer_distribution_matrix[file_kmer_to_scan_vec,start_point:end_point] <- t(TF_kmer_match_info_list_df[start_point:end_point,file_kmer_to_scan_vec])
          }                
        }
      read_frequency_vec <- rowSums( all_kmer_distribution_matrix ,na.rm=T)
      norm_all_kmer_distribution_matrix  <- all_kmer_distribution_matrix/read_frequency_vec

      max_freq_value <- round( max(norm_all_kmer_distribution_matrix, na.rm =T), digits=3)
      graph_name_Heatmap <- paste( sep="" , prefix_graph_dir,"Heatmap_" ,  prefix_file ,    "_Kmer_position_frequence.png")
      CairoPNG( 600,50*(2+dim(norm_all_kmer_distribution_matrix)[1]) ,file=    graph_name_Heatmap   , pointsize = 40/2 ,bg = "white" , res=72*2)
      pheatmap(    norm_all_kmer_distribution_matrix  , scale = "none",  cluster_cols = F, cluster_rows = F,
               main= paste(template_name, condition, ": Max Freq. =", max_freq_value),
               breaks=seq(0,  max_frequency,max_frequency/100),
               col=c( grey.colors(50,start=0, end=1, rev=T),redgreen(100)[50:1] ) ,
               key  = T , keysize=1, key.title= NA,  key.xlab="Rel.Freq.",
               symbreaks=T,   margins=c(8,20),   cexCol=1, annotation_legend=T,  border_color=NA, 
               labels_col = ifelse( is.element( 1:end_window, seq(10, end_window, 10)), 1:end_window, ""),
               cexRow=1 ,trace = "none",   density.info = "none", na_col="white" ) 
      dev.off()
      output_graph_name <- paste0( Graphs_dir, prefix_graph_name, "_",template_name,".png")
      if( prefix_graph_name !="Fig_2b_NCP"| template_name=="W601")
        {
        system( paste( "cp" , graph_name_Heatmap, output_graph_name ))
        write.table(  cbind( "start_bp"=rownames(norm_all_kmer_distribution_matrix ),norm_all_kmer_distribution_matrix), row.names=T, col.names = T, quote=F,       sep = "\t",
                     file=gsub( ".png" , "_datasets.txt",        output_graph_name  ))    
        }
      }
    }
  }

#### End ###