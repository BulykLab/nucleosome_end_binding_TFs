# DNA flexibility modules transcription factor binding to nucleosomes

## 0. Set up working directories and software paths
Open `config.R` file and update the paths to the PEAR software and other path-related variables.

## 1. Filtering, classification, and analysis of NCAP-SELEX data

### 1a. Downloading the SELEX data

### 1b. Find k-mer enrichment

### 1c. Determine k-mers used as TF binding motif

### 1d. Filter for end binding reads in each SELEX library

## 1.5 Predicting cyclizability of processed SELEX data 

Go to the `../DNAcycP/` directory for predicting sequence cyclizability. A bash script in the directory is prepared for submitting slurm jobs to predict all end binding sequence cyclizability in cloud computing.

## 2. Cyclizability analyses

### 2a. Load and summarize cyclizability prediction results

### 2b. Calculate the cyclizability slope in the internal sequence region

### 2c. Filter for TFs with high information gain in HT-SELEX library

### 2d. Filter for TFs with end binding enrichment and low end binding enrichment
