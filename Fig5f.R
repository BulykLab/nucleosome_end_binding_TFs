library(ShortRead)
library(Cairo)

head_sample_bed=read.delim(  "../data/head_sample.bed", header=F)
head_sample = readFasta(paste0(  "../data/head_sample.fasta"))
tail_sample_bed=read.delim( "../data/tail_sample.bed", header=F)
tail_sample_rc= tail_sample = readFasta(  "../data/tail_sample.fasta")
tail_sample_rc@sread =reverseComplement(tail_sample@sread)
sample_fasta = DNAStringSet(c(sread(head_sample),sread(tail_sample_rc)))
sample_bed <-  rbind( head_sample_bed, tail_sample_bed) 

head_control_bed=read.delim(  "../data/head_control.bed", header=F)
head_control = readFasta(paste0(  "../data/head_control.fasta"))
tail_control_bed=read.delim( "../data/tail_control.bed", header=F)
tail_control_rc= tail_control = readFasta(  "../data/tail_control.fasta")
tail_control_rc@sread =reverseComplement(tail_control@sread)
control_fasta = DNAStringSet(c(sread(head_control),sread(tail_control_rc)))
control_bed <-  rbind( head_control_bed, tail_control_bed) 


sample_dens = sapply(1:136, function(x) colMeans(oligonucleotideFrequency( narrow(  sample_fasta, x, x+4) , 5 )))
control_dens = sapply(1:136, function (x) colMeans(oligonucleotideFrequency(narrow( control_fasta, x, x+4 ), 5 )))

sample_combined_dens = colSums(sample_dens[c("TAAAA", "TTAAA" ,"TTTAA", "TTTTA"  ) ,] )
control_combined_dens = colSums( control_dens[c("TAAAA", "TTAAA" ,"TTTAA" , "TTTTA"  ) ,]   )

kmer_sorted_sample_freq = mean(sample_combined_dens[30:110])
kmer_sorted_control_freq =  mean(control_combined_dens[30:110])

plot(NA, xlab="TAAAA & TTAAA start position (bp)",ylab="TAAAA & TTAAA frequency",xlim=c(1,136),main="TA-rich 5mers", ylim=c(
  0, max(sample_combined_dens,control_combined_dens )))
points(sample_combined_dens, type="line", col="red")
points(control_combined_dens ,type="line", col="blue")
dev.copy(png,paste(sep="","Fig_5f.png"),300, 300)
dev.off()

CairoPDF( 3,3,file=paste(sep="","Fig_5f.pdf"))

plot(NA, xlab="TAAAA & TTAAA start bp",ylab="TAAAA & TTAAA frequency",xlim=c(1,136),main="TA-rich 5mers", ylim=c(
  0, max(sample_combined_dens,control_combined_dens )))
points(sample_combined_dens, type="line", col="red")
points(control_combined_dens ,type="line", col="blue")
dev.off()
