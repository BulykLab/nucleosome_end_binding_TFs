library(ggplot2)
library(readxl)
library(Cairo)

Parent_dir <-  $YOUR_PARENTAL_DIR  
#for example 
#Parent_dir <- "/Volumes/My_Book_3/PIONEAR_seq_Project_Back_Up_2022/Manuscript/preparation/Github/EMI_NCAP_SELEX/"
setwd(Parent_dir)
system( "mkdir MS_files/") 
system( "mkdir data/") 

# download the supplementary data from PMID: 30250250 
# and save the section a, b and c of Table S3 (EMI penetration) 
# in the "data" subdirectory with the following names:
# Zhu_et_al_Nature_2018_TABLE_S3_section_a.xlsx
# Zhu_et_al_Nature_2018_TABLE_S3_section_b.xlsx
# Zhu_et_al_Nature_2018_TABLE_S3_section_c.xlsx

#### Fig. S2e-h (Lig200)#####

EMI_lig200_df <- read_xlsx("data/Zhu_et_al_Nature_2018_TABLE_S3_section_b.xlsx", col_names = T)

EMI_lig200_matrix <- as.matrix(EMI_lig200_df[,-1])
row.names( EMI_lig200_matrix) <- EMI_lig200_df$TF
normalized_EMI_lig200_matrix  <-   EMI_lig200_matrix/rowSums(EMI_lig200_matrix )

TF_vec <- c("VSX1","SREBF2",    "CEBPB",  "FOXA2","GATA2" , 
            "GATA4","PAX7", "RUNX3", 
             "MYOG", "KLF12" ,   "KLF13")

protection_vec <- c()
thresh <- 0.0025
write.table( normalized_EMI_lig200_matrix[TF_vec,], file="MS_files/Fig_S2e_f_g_datasets.txt" ,row.names=T, col.names = T, quote=F,       sep = "\t")

for( index_TF in 1:length(TF_vec))
{

  TF_to_plot <-TF_vec[index_TF]
  smoothed_EMI_lig200 <- loess.smooth( 1:149, normalized_EMI_lig200_matrix[TF_to_plot ,] ,   span=0.45, degree=1)
  min_v <- min(smoothed_EMI_lig200$y)
  
  protected_region_lab <- normalized_EMI_lig200_matrix[TF_to_plot ,]>thresh+min_v
  protection <- 0
  temp_protection <- 0
  x0_temp <- 149
  for(index_protection in 1:length( protected_region_lab))
  { 
    if(! protected_region_lab[index_protection]){
      temp_protection = temp_protection + 1
      x0_temp <- min(x0_temp, index_protection)
      }
    if(protected_region_lab[index_protection]){

      if( temp_protection > protection ){
        protection <- temp_protection
        x0 <- x0_temp
        x1 <- x0_temp+  protection
        }
      x0_temp <- 149
      temp_protection = 0
      }
  }
  if(!protected_region_lab[index_protection]&temp_protection>protection){
    protection <- temp_protection
    x0 <- x0_temp
    x1 <- 149
    }
  
  figure_address <- paste0( "MS_files/Fig_S2e_f_g_" , TF_to_plot,".png")
  CairoPNG( file=  figure_address ,  pointsize=12, 
            width = 320, height = 240, bg = "white")
  plot (x=1:149, xlab="position on lig200 (bp)", ylab="Normalized E-MI diagonals [AU]", 
        ylim=c(0,max(normalized_EMI_lig200_matrix[TF_to_plot ,] ) ), main=paste( TF_to_plot, "; protected bases",protection )) 
  lines(normalized_EMI_lig200_matrix[TF_to_plot ,])
  segments( x0, thresh+min_v, x1 ,  col="red")
  dev.off() 
  
  protection_vec[index_TF] <-   protection
}

names(  protection_vec ) <- TF_vec
protection_vec <- sort(  protection_vec)

CairoPNG( width = 320,  height = 240, 
          file=paste( sep="", "MS_files/Fig_S2h.png") ,  pointsize=12, 
          bg = "white")
barplot(   protection_vec, las=2)
dev.off()
write.table( protection_vec, file="MS_files/Fig_S2h_datasets.txt" ,col.names=T,  quote=F,       sep = "\t")


##### Fig S8b, Lig147 ####

for( SELEX_data in c( "NCAP", "HT") )
  {
  if(SELEX_data=="NCAP" )
  { 
    section="a"
    TF_vec <- c(   "CEBPB",  "FOXA2", "GATA2" , "GATA4" , "PAX7", "RUNX3"      )
  }
  if(SELEX_data=="HT" )
  { 
    section="c"
    TF_vec <- c(   "CEBPB",  "FOXA2",  "GATA4" , "PAX7", "RUNX3"      )
  }
  file_address <- paste0(   "data/Zhu_et_al_Nature_2018_TABLE_S3_section_", section,".xlsx")
  EMI_lig147_df <- read_xlsx( file_address , col_names = T)

  EMI_lig147_matrix <- as.matrix(EMI_lig147_df[,-1])
  row.names( EMI_lig147_matrix) <- EMI_lig147_df$TF
  normalized_EMI_lig147_matrix  <-   EMI_lig147_matrix/rowSums(EMI_lig147_matrix )
  normalized_EMI_lig147_matrix[TF_vec,]
  write.table( normalized_EMI_lig147_matrix[TF_vec,], file=paste0("MS_files/Fig_S8b_", SELEX_data,"-SELEX_datasets.txt") ,row.names=T, col.names = T, quote=F,       sep = "\t")
  
  for( index_TF in 1:length(TF_vec))
  {
    TF_to_plot <-TF_vec[index_TF]
  file_address <- paste0( "MS_files/Fig_S8b_" , TF_to_plot, "_", SELEX_data,"-SELEX.png")
    CairoPNG( file=   file_address,  pointsize=12, 
            width = 320, height = 240, bg = "white")
  plot (x=1:96, xlab="position on lig147 (bp)", ylab="Normalized E-MI diagonals [AU]", 
        ylim=c(0,max(normalized_EMI_lig147_matrix[TF_to_plot ,] ) ), main=paste( TF_to_plot, SELEX_data) )
  lines(normalized_EMI_lig147_matrix[TF_to_plot ,])
   dev.off() 

}
}
