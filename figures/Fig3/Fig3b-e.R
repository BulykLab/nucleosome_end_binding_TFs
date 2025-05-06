library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)


# Panel b
load("../data/MNase_VS_Input_Read_percentage.rdata") ##### NEED TO CHANGE #####
lib_percentage_mat=as.data.frame(percentage_read_matrix)
lib_percentage_mat$exp = rownames(lib_percentage_mat)
lib_percentage_mat_long =  pivot_longer(lib_percentage_mat, cols=-exp)

lib_percentage_mat_long$group = grepl( "Mnase",lib_percentage_mat_long$exp)
avg_lib = lib_percentage_mat_long %>% filter(grepl( "Mnase",lib_percentage_mat_long$exp)) %>% group_by(name) %>% summarise(avg=mean(value)) 

panel_b=ggplot(avg_lib, aes(x=reorder(name,-avg), y=avg, fill=name))+geom_bar(stat="identity", show.legend = F)+theme_classic()+
  geom_point(data=lib_percentage_mat_long %>% filter(group), aes(x=name, y=value, shape=exp), color="black",show.legend = F, size=2)+
  labs(x="Library template sequence", y="Protected reads (%)")+
  scale_fill_brewer(palette = "Set2")+theme(text=element_text(size=7), title=element_text(size=7))


# Panel d
library(ggplot2)
library(dplyr)

avg_plot_df = read.csv("../data/pionear_seq_library_cyc_summary.csv")
panel_d=ggplot(avg_plot_df,aes(x=factor(template, levels=c("W601", "NRCAM", "ALB1", "CX3")), fill=template))+
  geom_boxplot(aes(lower=mean-sd,upper=mean+sd,middle=mean,ymin=min,ymax=max),stat="identity", show.legend = F)+
  labs(x="Library template sequence", y="Average sequence cyclizability")+theme_classic()+scale_fill_brewer(palette="Set2")+
  theme(text=element_text(size=7), title=element_text(size=7))



# Panel e


ALB1 = read.table("../data/ALB1_prediction.txt")
CX3 = read.table("../data/CX3_prediction.txt")
NRCAM = read.table("../data/NRCAM_prediction.txt")
W601 = read.table("../data/W601_prediction.txt")

plot_df = rbind(
  data.frame(x=25:122, y=ALB1$V2, Template="ALB1"),
  data.frame(x=25:122, y=CX3$V2, Template="CX3"),
  data.frame(x=25:122, y=NRCAM$V2, Template="NRCAM"),
  data.frame(x=25:122, y=W601$V2, Template="W601")
)

panel_e=ggplot(plot_df, aes(x=x , y=y,color=Template))+geom_line(show.legend = F, linewidth=1)+
  facet_grid(cols = vars(factor(Template, levels=c("W601", "NRCAM", "ALB1", "CX3"))))+
  theme_bw()+labs(x="50-bp center position (bp)", y="Sequence cyclizability")+scale_color_brewer(palette="Set2")+
  theme(strip.background =element_rect(fill="white"))+
  theme(text=element_text(size=7), title=element_text(size=7))


design <- "
  555
  123
  444
  666
"

panel_b+plot_spacer()+panel_d+panel_e+plot_spacer()+plot_spacer()+plot_layout(design = design, widths = c(1,2,1), heights=c(2,1.5,1,1.5))
ggsave("Fig3.png", units="mm", width = 180, heigh=170)
