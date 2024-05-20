#!/bin/bash
echo $1 # Nucleosome positions
echo $2 # Motif positions
echo $3 # Motif length
echo $4 # Output dir

module load gcc bedtools homer

bedtools intersect -wa -a $1 -b $2 > with_motif.bed
awk -v OFS='\t' -v motif_len=$3 '{print $1, $2+motif_len-1, $2+30-motif_len+1}' with_motif.bed > with_chip_head.bed
awk -v OFS='\t' -v motif_len=$3 '{print $1, $3-30+motif_len-1, $3-motif_len+1}' with_motif.bed > with_chip_tail.bed
bedtools intersect -wa -a with_chip_head.bed -b $2 > head_chip.bed
bedtools intersect -wa -a with_chip_tail.bed -b $2 > tail_chip.bed
awk -v OFS='\t' -v motif_len=$3 '{print $1, $2-motif_len, $2+140-motif_len}' head_chip.bed > head_chip_nuc.bed
awk -v OFS='\t' -v motif_len=$3 '{print $1, $3+motif_len-140, $3+motif_len}' tail_chip.bed > tail_chip_nuc.bed
bedtools intersect -v -wa -a head_chip_nuc.bed -b tail_chip_nuc.bed | uniq -u> head_chip_nuc_exc.bed
bedtools intersect -v -wa -a tail_chip_nuc.bed -b head_chip_nuc.bed | uniq -u > tail_chip_nuc_exc.bed
awk -v OFS='\t' '{if ($3-$2==140){print $1, $2, $3}}' head_chip_nuc_exc.bed > head_motif_nuc.bed
awk -v OFS='\t' '{if ($3-$2==140){print $1, $2, $3}}' tail_chip_nuc_exc.bed > tail_motif_nuc.bed
mkdir $4
homerTools extract head_motif_nuc.bed hg38 -fa > "${4}/head_nucs.fasta"
homerTools extract tail_motif_nuc.bed hg38 -fa > "${4}/tail_nucs.fasta"

rm with_chip_head.bed
rm with_chip_tail.bed
rm head_chip.bed
rm tail_chip.bed
rm head_chip_nuc.bed
rm tail_chip_nuc.bed 
rm head_chip_nuc_exc.bed
rm tail_chip_nuc_exc.bed
