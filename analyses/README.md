# DNA flexibility modules transcription factor binding to nucleosomes

## 0. Set up working directories and software paths
Open `config.R` file and update the paths to the PEAR software and other path-related variables.

## 1. Filtering, classification, and analysis of NCAP-SELEX data

- a. Downloading the SELEX data with `1a_download_SELEX.R`

- b. Find k-mer enrichment with `1b_kmer_enrichment.R`

- c. Determine k-mers used as TF binding motif with `1c_find_binding_kmers.R`

- d. Filter for end binding reads in each SELEX library with `1d_find_end_binding_reads.R`

- e. Filter for internal binding reads in each NCAP-SELEX library with `1e_find_internal_end_finding_reads.R`
- f. Compute nucleotide spatial periodicity with FFT with `1f_FFT_end_binding_reads.R`

## 1.5 Predicting cyclizability of processed SELEX data 

Go to the directory `../Cyclizability-Prediction-Website/` for predicting sequence cyclizability of end binding reads in SELEX libraries. It is recommended to submit parallel jobs on a cloud computing platform. Example bash scripts are prepared for submitting slurm jobs  in the directory to predict all end binding sequence cyclizability in cloud computing, including `pred_fasta.sh` and `pred_fasta_cpu.sh`. To run the script with CPU, under the directory `../Cyclizability-Prediction-Website/`,

```{bash}
chmod 777 selex_cyc_script.sh # This is a script that will be generated in 1d_find_end_binding_reads.R
./selex_cyc_script.sh
chmod 777 selex_internal_cyc_script.sh # This is a script that will be generated in 1e_find_internal_end_binding_reads.R
./selex_internal_cyc_script.sh
```

## 2. Cyclizability analyses of SELEX data

- a. Load and summarize cyclizability prediction results with `2a_load_SELEX_cyc.R`

- b. Calculate the cyclizability slope in the internal sequence region with `2b_calc_cyc_slope_avg.R`

- c. Filter for TFs with high information gain in HT-SELEX library with `2c_info_gain_filter.R`

- d. Filter for TFs with end binding enrichment and low end binding enrichment with `2d_end_binding_filter.R`

## 3. Cyclizability analysis of K562 nucleosomes bound\not bound by CEBPB

- Download MNase-seq data from ENCODE
In the `data` directory at the top level of this repo, 
```{bash}
wget https://www.encodeproject.org/files/ENCFF000VMJ/@@download/ENCFF000VMJ.bam
```
- Index the BAM file with samtools

```{bash}
samtools index ENCFF000VMJ.bam
```

- Call nucleosome peaks using DANPOS3

DANPOS3 can be installed here: [https://github.com/sklasfeld/DANPOS3](https://github.com/sklasfeld/DANPOS3).
Go to the DANPOS directory, then run the following commands:

```{bash}
python danpos.py dpos <path/to/data/ENCFF000VMJ.bam> -o <path/to/data/k562_nucleosomes/>
```
In the `data` directory at the top level of this repo, convert DANPOS output `.xls` file to bed format.

```{bash}
awk -v OFS='\t' '{if (NR>1) {print $1, $2, $3, $4, $5}} k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.xls > k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.bed
```

- Verify nucleosome positions are not biased

Required packages: bedtools, bedops, seqoutbias

Seqoutbias can be downloaded here: [https://guertinlab.github.io/seqOutBias/](https://guertinlab.github.io/seqOutBias/).


```{bash}
 ./seqOutBias hg19.fa <path/to/data/ENCFF000VMJ.bam> --no-scale --skip-bw --stranded --shift-counts --pdist=100:400

# Postprocessing reads to signal
awk -v OFS='\t' '{if ($6=="+") {print $1, $2, $2+35, $4, $5, $6} else {print $1, $3-35,$3, $4, -$5, $6}}' <path/to/data/ENCFF000VMJ_scaled.bed> > <path/to/data/ENCFF000VMM_scaled_pp.bed>
sort-bed <path/to/data/ENCFF000VMM_scaled_pp.bed> > <path/to/data/ENCFF000VMM_scaled_pp_sorted.bed>
bedops --partition <path/to/data/ENCFF000VMM_scaled_pp_sorted.bed> | bedmap --echo --sum --delim '\t'  - <path/to/data/ENCFF000VMM_scaled_pp_sorted.bed> > <path/to/data/ENCFF000VMM_scaled_pp_sorted_merged.bedgraph>
```

Positions can be compared using a genome browser with <path/to/data/ENCFF000VMM_scaled_pp_sorted_merged.bedgraph> and the original bam file.

- Convert hg19 nucleosomes to hg38 coordinates

Download liftOver executable here: [https://genome-store.ucsc.edu/](https://genome-store.ucsc.edu/) and the hg19 to hg38 chain file here: [https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz).

In the `data` directory at the top level of this repo,
```{bash}
<liftover_executable> data/k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.bed \
  <path/to/hg19ToHg38.over.chain.gz> \
  data/k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.hg38.bed \
  unMapped
```
- Filtering for nucleosome peaks with 140-bp
In the `data` directory at the top level of this repo,
```{bash}
awk -v OFS='\t' '{if($3-$2==140){print $1,$2,$3}}' k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.hg38.bed > k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.hg38.L140.bed
```

- Download ChIP-seq data
In the `data` directory at the top level of this repo, 
```{bash}
wget https://www.encodeproject.org/files/ENCFF712ZNR/@@download/ENCFF712ZNR.bed.gz
gunzip ENCFF712ZNR.bed.gz
mv ENCFF712ZNR.bed IDR_thres_peaks_ENCFF712ZNR.bed
wget https://www.encodeproject.org/files/ENCFF309ZJM/@@download/ENCFF309ZJM.bed.gz
gunzip ENCFF309ZJM.bed.gz
mv ENCFF309ZJM.bed all_peaks_ENCFF309ZJM.bed
```

- Use Homer to find CEBPB binding motifs and their positions.
Homer can be downloaded and installed following instructions here: [http://homer.ucsd.edu/homer/](http://homer.ucsd.edu/homer/)
In the `data` directory at the top level of this repo, 
```{bash}
findMotifsGenome.pl IDR_thres_peaks_ENCFF712ZNR.bed hg38 ../data/CEBPB_homer_out/ -size given \
  -preparsedDir preparsed
scanMotifGenomeWide.pl CEBPB_homer_out/homerResults/motif1.motif hg38 \
  -bed > CEBPB_homer_out/CEBPB_motif1_hg38.bed
```

- Separate bound/unbound motifs
In the `data` directory at the top level of this repo, 
```{bash}
bedtools intersect -wa -a CEBPB_motif1_hg38.bed -b IDR_thres_peaks_ENCFF712ZNR.bed -f 1> CEBPB_motif1_with_IDR_thres_peaks_hg38.bed

bedtools intersect -v -wa -a CEBPB_motif1_hg38.bed -b all_peaks_ENCFF309ZJM.bed > CEBPB_motif1_no_peak_hg38.bed
```

- Separate bound/unbound nucleosomes
In the `data` directory at the top level of this repo, 
```{bash}
../analysis/3a_get_binding_nuc.sh  k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.hg38.L140.bed CEBPB_motif1_with_IDR_thres_peaks_hg38.bed sample
../analysis/3a_get_binding_nuc.sh  k562_nucleosomes/pooled/ENCFF000VMJ.smooth.positions.hg38.L140.bed CEBPB_motif1_no_peak_hg38.bed control
```

- Extract nucleosome sequence and predict for cyclizability
In the `data` directory, the bound nucleosomes are stored in `head_sample.bed` and `tail_sample.bed`, while the unbound nucleosomes are stored in `head_control.bed` and `tail_control.bed`.
The following command is used to extract bed file sequences. 
```{bash}
homerTools extract <bed-file> hg38 -fa > <fasta-file> 
```
Cyclizability is predicted by using `pred_fasta_cpu.sh` script in the directory `Cyclizability-Prediction-Website/`.
