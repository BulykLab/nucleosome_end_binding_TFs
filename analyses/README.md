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
chmod 777 selex_internal_cyc_script.sh # This is a script that will be generated in 1e_find_internal_end_binding_reads.R
./selec_internal_cyc_script.sh
```

## 2. Cyclizability analyses of SELEX data

### 2a. Load and summarize cyclizability prediction results

### 2b. Calculate the cyclizability slope in the internal sequence region

### 2c. Filter for TFs with high information gain in HT-SELEX library

### 2d. Filter for TFs with end binding enrichment and low end binding enrichment

## 3. Cyclizability analysis of K562 nucleosomes bound\not bound by CEBPB

- Download MNase-seq data from ENCODE

```{bash}
wget https://www.encodeproject.org/files/ENCFF000VMJ/@@download/ENCFF000VMJ.bam
mv ENCFF000VMJ.bam ../data/
```

- Call nucleosome peaks using DANPOS3

DANPOS3 can be installed here: [https://github.com/sklasfeld/DANPOS3](https://github.com/sklasfeld/DANPOS3).
Go to the DANPOS directory, then run the following commands:

```{bash}
python danpos.py dpos <path/to/data/ENCFF000VMJ.bam> -o <path/to/data/k562_nucleosomes/>
```
Convert DANPOS output `.xls` file to bed format.

```{bash}
awk -v OFS='\t' '{if (NR>1) {print $1, $2, $3, $4, $5}} \
  <path/to/data/k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.xls> > \
  <path/to/data/k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.bed>'
```

- Convert hg19 nucleosomes to hg38 coordinates

Download liftOver executable here: [https://genome-store.ucsc.edu/](https://genome-store.ucsc.edu/) and the hg19 to hg38 chain file here: [https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz).

```{bash}
<liftover_executable> <path/to/data/k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.bed> \
  hg19ToHg38.over.chain.gz \
  <path/to/data/k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.hg38.bed> \
  unMapped
```

- Download ChIP-seq data

```{bash}
wget https://www.encodeproject.org/files/ENCFF712ZNR/@@download/ENCFF712ZNR.bed.gz
mv ENCFF712ZNR.bed.gz ../data/
gunzip ../data/ENCFF712ZNR.bed.gz
```

- Use Homer to find CEBPB binding motifs and their positions.
Homer can be downloaded and installed following instructions here: [http://homer.ucsd.edu/homer/](http://homer.ucsd.edu/homer/)

```{bash}
findMotifsGenome.pl ../data/ENCFF712ZNR.bed hg38 ../data/CEBPB_homer_out/ -size given \
  -preparsedDir preparsed
scanMotifGenomeWide.pl ../data/CEBPB_homer_out/homerResults/motif1.motif hg38 \
  -bed > ../data/CEBPB_homer_motif_hg38.bed
```

- Separate bound/unbound nucleosomes

```{bash}
./3a_get_binding_nuc.sh  ../data/k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.bed> \
  ../data/CEBPB_homer_motif_hg38.bed 8 ../data/
```

- Extract nucleosome sequence and predict for cyclizability