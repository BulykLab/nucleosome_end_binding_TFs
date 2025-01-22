source("utils.R")
source("config.R")
options(timeout = 500)

tf_list = c("SOX11")#read.table(SELEX_library_file)$V1


### Modify path here
data_dir <- paste0(SELEX_data_dir, SELEX_ligand, "/")

# Create the data folder if it does not exist
if (!file.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}



for (tf_name in tf_list){
  # Check if downloading input library
  input <- (tf_name == "cycle-0")
  
  tf_data_dir <- paste0(data_dir,tf_name,"/")
  
  # Create the folder for TF if it does not exist
  if (!file.exists(tf_data_dir)) {
    dir.create(tf_data_dir)
  }
  
  # Read report
  report_tbl <- read.csv2(SELEX_data_report_tsv,
                          sep = "\t")
  
  # Find all urls for tf within ligand library
  if (input){
    tf_tbl <- report_tbl %>%
      filter(grepl(SELEX_ligand, submitted_ftp)) %>%
      filter(grepl(paste0(tf_name,"_"), submitted_ftp)) 
  } else {
    tf_tbl <- report_tbl %>%
      filter(grepl(SELEX_ligand, submitted_ftp)) %>%
      filter(grepl(paste0(tf_name,"_"), submitted_ftp)) %>%
      filter(grepl("cycle-4", submitted_ftp))
  }
  
 
  
  # Find destinations of the downloaded fastq
  get_dest <- function(sftp) {
    sftps <- str_split(sftp, ";")
    dests <- sapply(sftps, function(x) {
      dest <- gsub("^.*/", "", x)
      return(paste0(tf_data_dir, dest))
    })
    return(paste(dests, collapse = ";"))
  }
  
  destinations <- sapply(tf_tbl$submitted_ftp, get_dest)
  
  # Download two fastq_ftps to the destinations
  download_urls <- function(fftp, dest) {
    fftps <- str_split(fftp, ";")
    dests <- str_split(dest, ";")
    if (!file.exists(dests[[1]][1])) {
      download.file(fftps[[1]][1], dests[[1]][1])
    }
    if (!file.exists(dests[[1]][2])) {
      download.file(fftps[[1]][2], dests[[1]][2])
    }
  }
  
  mapply(download_urls, tf_tbl$fastq_ftp, destinations)
  write.table(destinations, row.names = F, paste0(tf_data_dir, "dest.txt"))
  
  # Call PEAR to merge the downloaded fastqs
  merge_pe_reads <- function(dests) {
    read_files <- str_split(dests, ";")
    if (input) {
      outputPrefix <- paste0(tf_data_dir, glue("input_{SELEX_ligand}_cycle-0"))
      print(outputPrefix)
      
      ### Modify path for PEAR
      system(glue("{pear_path} -f {read_files[[1]][1]} -r {read_files[[1]][2]} -v 5 -o {outputPrefix}"))
      return (paste0(outputPrefix,".assembled.fastq"))
    } else {
      infos <- str_split(gsub("^.*/", "", read_files[[1]][1]), "__")
      selex_method <- infos[[1]][1]
      cycle <- infos[[1]][4]
      if (cycle == "cycle-5") {
        output_dir <- paste0(tf_data_dir, selex_method, "_", cycle, "_",
                             infos[[1]][5], "/")
        outputPrefix <- paste0(infos[[1]][1:5], collapse = "_")
      } else {
        output_dir <- paste0(tf_data_dir, selex_method, "_", cycle, "/")
        outputPrefix <- paste0(infos[[1]][1:4], collapse = "_")
      }
      outputPrefix <- paste0(output_dir, outputPrefix)
      if (file.exists(paste0(outputPrefix,".assembled.fastq"))) {
        return(paste0(outputPrefix,".assembled.fastq"))
      }
      dir.create(output_dir)
      ### Modify path for PEAR
      system(glue("{pear_path} -f {read_files[[1]][1]} -r {read_files[[1]][2]} -v 5 -o {outputPrefix}"))
      return(paste0(outputPrefix,".assembled.fastq"))
    }
  }
  
  merged_dests = sapply(destinations, merge_pe_reads)
  
  # Create clean destinations
  clean_dests = sapply(merged_dests, function(x)gsub(".assembled",".assembled.clean",x))
  
  # Clean invalid and duplicated reads
  clean_dedup = function(origin_fastq, output){
    if (file.exists(output)){
      return (0)
      file.remove(output)
    }
    data_vec <-readFastq(origin_fastq)
    print(length(data_vec))
    
    N101_data = data_vec[width(data_vec)==101]
    print(length(N101_data))
    cleaned_data <- clean(N101_data)
    # How many reads remain after getting rid of "N" reads?
    print(length(cleaned_data))
    v <- unique(sread(cleaned_data))
    print(length(v))
    # Now writing up Clean version###
    unique_reads_fastq <- N101_data[match(v, N101_data@sread)]
    print(length(unique_reads_fastq))
    writeFastq(unique_reads_fastq, output, mode="w",
               compress=FALSE, format="fastq", qualityType = "FastqQuality")
  }
  
  mapply(clean_dedup,merged_dests,clean_dests)
}
