# DNA flexibility modules transcription factor binding to nucleosomes

## 0. Set up working directories and software paths
Open `config.R` file and update the paths to the PEAR software and other path-related variables.

## 1. Filtering, classification, and analysis of NCAP-SELEX data

### 1a. Downloading the SELEX data

### 1b. Find k-mer enrichment

### 1c. Determine k-mers used as TF binding motif

### 1d. Filter for end binding reads in each SELEX library

## 1.5 Predicting cyclizability of processed SELEX data 

Go to the directory `../Cyclizability-Prediction-Website/` for predicting sequence cyclizability of end binding reads in SELEX libraries. It is recommended to submit parallel jobs on a cloud computing platform. Example bash scripts are prepared for submitting slurm jobs  in the directory to predict all end binding sequence cyclizability in cloud computing, including `pred_fasta.sh` and `pred_fasta_cpu.sh`. To run the script with CPU, under the directory `../Cyclizability-Prediction-Website/`,

```{bash}
chmod 777 selex_cyc_script.sh # This is a script that will be generated in 1d_find_end_binding_reads.R
./selex_cyc_script.sh
```

## 2. Cyclizability analyses

### 2a. Load and summarize cyclizability prediction results

### 2b. Calculate the cyclizability slope in the internal sequence region

### 2c. Filter for TFs with high information gain in HT-SELEX library

### 2d. Filter for TFs with end binding enrichment and low end binding enrichment
