library(SELEX)

source("utils.R")
source("config.R")

tf_list = c("SOX11")#read.table(SELEX_tf_list_file)$V1
var_length <- 101
exprs = c("NCAP-SELEX_cycle-4", "HT-SELEX_cycle-4")

input_name <- "cycle-0" ### dir of input library
input_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/", input_name, "/")

# Store library information gain
info_gain_df = data.frame(matrix(ncol=3,nrow=0))



for (tf_name in tf_list){
  # Define variables
  tf_data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/", tf_name, "/")
  ef_threshold = 1 ### Enrichment threshold to eliminate non-enriched kmers
  
  for (expr in exprs){
    if (!dir.exists(paste0(tf_data_dir, expr))){
      next
    }
    
    # Set up selex config
    ### Modify path for cache
    selex.config(workingDir = SELEX_cache_dir, verbose = F, maxThreadNumber = 4)
    
    # Retrieve cleaned merged reads
    get_clean_file <- function(ddir) {
      all_files <- list.files(ddir)
      clean_files <- all_files[grepl("clean", all_files)]
      return(paste0(ddir, clean_files[1]))
    }
    
    # Load input sample
    input_clean_file <- get_clean_file(input_dir)
    input_seq_name <- paste(SELEX_ligand, input_name, sep = "_")
    selex.defineSample(
      seqName = input_seq_name, seqFile = input_clean_file,
      sampleName = input_name, round = 0, rightBarcode = "",
      leftBarcode = "", varLength = var_length
    )
    r0 <- selex.sample(seqName = input_seq_name, sampleName = input_name, 
                       round = 0)
    
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
      sample_seq_name <- paste(SELEX_ligand, tf_name, sample_dir, sep = "_")
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
    
    k = 7
    # Get information gain
    tf_sample_info_gains <- as.data.frame(sapply(
      tf_selex_samples,
      function(s) {
        ig = selex.infogain(s, k = c(k), markovModel = mm)
        print(ig)
        print(k)
        return (ig)
      }
    ))
    info_gain_df = rbind(info_gain_df, data.frame(TF=tf_name, 
                                                  Experiment=expr, 
                                                  InfoGain=tf_sample_info_gains[[1]]))
    
    # Get affinities
    get_affinites <- function(sample, sample_dir) {
      for (k in c(7,8,9)){
        ig = selex.infogain(sample, k=c(k), markovModel = mm)
        print(ig)
        
        affinities <- selex.affinities(sample, k = k, markovModel = mm)
        # Calculate K-mer enrichment fold and filter for enriched k-mers
        affinities <- affinities %>% 
          mutate(enrichmentFold = ObservedCount / ExpectedCount)
        sorted_affinities = affinities[order(affinities$Affinity,
                                             decreasing = TRUE), ]
        potential_affinities = 
          sorted_affinities[sorted_affinities$enrichmentFold>ef_threshold,]
        write.csv(potential_affinities,
                  file = paste0(tf_data_dir, sample_dir, 
                                "/aff_mat_potential_", k,".txt"), quote = F
        )
      }
      return(affinities)
    }
    aff_mats <- mapply(get_affinites, tf_selex_samples, tf_sample_dirs)
  }
}

# Save information gain to file
write.csv(info_gain_df, info_gain_path)