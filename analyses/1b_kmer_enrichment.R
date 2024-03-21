library(SELEX)
library(stringr)
library(dplyr)

source("utils.R")
source("config.R")

tf_list = read.table(SELEX_tf_list_file)$V1
var_length <- 101
exprs = c("NCAP-SELEX_cycle-4")

input_name <- "cycle-0" ### dir of input library
input_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/", input_name, "/")

info_gain_df = data.frame(matrix(ncol=3,nrow=0))
info_gain_path <- "../info_gain.txt"


for (tf_name in tf_list){
  # Define variables
  ### Modify paths, ligands, random region lengths
  tf_data_dir <- paste0(SELEX_data_dir, ligand, "/", tf_name, "/")
  ef_threshold = 1 ### Enrichment threshold to eliminate non-enriched kmers
  
  for (expr in exprs){
    if (!dir.exists(paste0(tf_data_dir, expr))){
      next
    }
    
    # Set up selex config
    ### Modify path for cache
    selex_wkdir <- "/n/data2/bch/medicine/bulyk/Katrina/capstone/cache/"
    selex.config(workingDir = selex_wkdir, verbose = F, maxThreadNumber = 4)
    
    # Retrieve cleaned merged reads
    get_clean_file <- function(ddir) {
      all_files <- list.files(ddir)
      clean_files <- all_files[grepl("clean", all_files)]
      return(paste0(ddir, clean_files[1]))
    }
    
    # Load input sample
    input_clean_file <- get_clean_file(input_dir)
    input_seq_name <- paste(ligand, input_name, sep = "_")
    selex.defineSample(
      seqName = input_seq_name, seqFile = input_clean_file,
      sampleName = input_name, round = 0, rightBarcode = "",
      leftBarcode = "", varLength = 101 
    )
    r0 <- selex.sample(seqName = input_seq_name, sampleName = input_name, round = 0)
    
    # Build Markov model
    r0.split <- selex.split(r0)
    k <- 6
    mm <- selex.mm(
      sample = r0.split$train, order = NA,
      crossValidationSample = r0.split$test, Kmax = k
    )
    
    
    
    # Load TF sample directories
    ### Modify experiment directory
    tf_sample_dirs <- c(expr)
    
    get_samples <- function(sample_dir) {
      # Get sample information
      seq_file <- get_clean_file(paste0(tf_data_dir, sample_dir, "/"))
      cycle <- str_split(sample_dir, "_")[[1]][2]
      round <- strtoi(substr(cycle, nchar(cycle), nchar(cycle)))
      sample_seq_name <- paste(ligand, tf_name, sample_dir, sep = "_")
      # Load sample
      selex.defineSample(
        seqName = sample_seq_name, seqFile = seq_file, sampleName = sample_dir,
        round = round, varLength = var_length, leftBarcode = "",
        rightBarcode = ""
      )
      sample <- selex.sample(
        seqName = sample_seq_name, sampleName = sample_dir,
        round = round
      )
      return(sample)
    }
    
    tf_selex_samples <- sapply(tf_sample_dirs, get_samples)
    
    # Get information gain
    tf_sample_info_gains <- as.data.frame(sapply(
      tf_selex_samples,
      function(s) {
        ig = selex.infogain(s, k = c(k), markovModel = mm)
        print(ig)
        return (ig)
      }
    ))
    
    info_gain_df = rbind(info_gain_df, data.frame(TF=tf_name, Experiment=expr, InfoGain=tf_sample_info_gains[[1]]))
    
    # Get affinities
    get_affinites <- function(sample, sample_dir) {
      for (k in c(7,8,9)){
        ig = selex.infogain(sample, k=c(k), markovModel = mm)
        print(ig)
        
        affinities <- selex.affinities(sample, k = k, markovModel = mm)
        affinities <- affinities %>% mutate(enrichmentFold = ObservedCount / ExpectedCount)
        sorted_affinities = affinities[order(affinities$Affinity,
                                             decreasing = TRUE), ]
        potential_affinities = sorted_affinities[sorted_affinities$enrichmentFold>ef_threshold,]
        print(head(sorted_affinities))
        write.csv(potential_affinities,
                  file = paste0(tf_data_dir, sample_dir, "/aff_mat_potential_", k,".txt"), quote = F
        )
      }
      return(affinities)
    }
    aff_mats <- mapply(get_affinites, tf_selex_samples, tf_sample_dirs)
  }
}

write.csv(info_gain_df, info_gain_path)