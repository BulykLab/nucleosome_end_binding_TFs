#!/bin/bash
echo $1 # Nucleosome Positions
echo $2 # Motif Positions
echo $3 # Sample Name

module load gcc bedtools homer

bedtools intersect -wa -a $1 -b $2 > with_motif.bed
awk -v OFS='\t' -v motif_len=$3 '{print $1, $2, $2+30}' with_motif.bed > with_chip_head.bed
awk -v OFS='\t' -v motif_len=$3 '{print $1, $3-30, $3}' with_motif.bed > with_chip_tail.bed
bedtools intersect -wa -a with_chip_head.bed -b $2 -F 1.0 | uniq > head_chip.bed
bedtools intersect -wa -a with_chip_tail.bed -b $2 -F 1.0 | uniq > tail_chip.bed
awk -v OFS='\t' -v motif_len=$3 '{print $1, $2, $2+140}' head_chip.bed > "head_${3}.bed"
awk -v OFS='\t' -v motif_len=$3 '{print $1, $3-140, $3}' tail_chip.bed > "tail_${3}.bed"
homerTools extract "head_${3}.bed" hg38 -fa > "head_${3}.fasta"
homerTools extract "tail_${3}.bed" hg38 -fa > "tail_${3}.fasta"

rm with_motif.bed
rm with_chip_head.bed
rm with_chip_tail.bed
rm head_chip.bed
rm tail_chip.bed
