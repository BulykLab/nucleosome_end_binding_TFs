library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)


# Panel a
load("../data/MNase_VS_Input_Read_percentage.rdata") ##### NEED TO CHANGE #####
lib_percentage_mat=as.data.frame(percentage_read_matrix)
lib_percentage_mat$exp = rownames(lib_percentage_mat)
lib_percentage_mat_long =  pivot_longer(lib_percentage_mat, cols=-exp)

lib_percentage_mat_long$group = case_when(
  grepl("1-Mnase", lib_percentage_mat_long$exp) ~ "After digestion replicate 1 (   )",
  grepl("2-Mnase", lib_percentage_mat_long$exp) ~ "After digestion replicate 2 (   )",
  grepl("1-Input", lib_percentage_mat_long$exp) ~ "Before digestion replicate 1",
  grepl("2-Input", lib_percentage_mat_long$exp) ~ "Before digestion replicate 2"
)


avg_lib = lib_percentage_mat_long %>% group_by(name, group) %>% summarise(avg=mean(value)) 
lib_percentage_mat_long$label = paste0(round(lib_percentage_mat_long$value, digits=2), "%")

panel_a = ggplot(lib_percentage_mat_long, 
       aes(x=factor(group, levels=rev(c("Before digestion replicate 1", 
                                        "Before digestion replicate 2", 
                                        "After digestion replicate 1 (   )",
                                        "After digestion replicate 2 (   )"))), 
           y=value/100, fill=name, label=label))+
  geom_bar(stat="identity")+
  geom_text(size = 2, position = position_stack(vjust = 0.5))+
  scale_fill_brewer(palette = "Set2")+
  scale_y_continuous(labels = scales::percent)+
  theme_classic()+labs(x="", y="Percentage of reads in library", fill="Template")+
  theme(axis.text.x = element_text(angle=90, vjust=0.5))+coord_flip()+
  theme(text=element_text(size=7), title=element_text(size=7))



# Panel b
template_cyc_dir = "../data/pionear_seq_templates/_cycle_"
ALB1 = read.csv(paste0(template_cyc_dir, "ALB1.txt"))
CX3 = read.csv(paste0(template_cyc_dir, "CX3.txt"))
NRCAM = read.csv(paste0(template_cyc_dir, "NRCAM.txt"))
W601 = read.csv(paste0(template_cyc_dir, "W601.txt"))

ALB1$template="ALB1"
CX3$template="CX3"
NRCAM$template="NRCAM"
W601$template="W601"

plot_df = rbind(ALB1, CX3, NRCAM, W601)

panel_b = ggplot(plot_df, aes(x=posision, y=c_score_norm, color=template))+geom_line(show.legend = F, linewidth=1)+
  facet_grid(cols=vars(factor(template, levels=c("W601", "NRCAM", "ALB1", "CX3"))))+theme_bw()+
  labs(x="50-bp center position (bp)", y="Sequence cyclizability")+
  scale_color_brewer(palette="Set2")+theme(strip.background =element_rect(fill="white"))+
  theme(text=element_text(size=7), title=element_text(size=7))


# Panel c
ALB1 = read.csv("../data/ALB1_library_cyc.csv")
CX3 = read.csv("../data/CX3_library_cyc.csv")
NRCAM = read.csv("../data/NRCAM_library_cyc.csv")
W601 = read.csv("../data/W601_library_cyc.csv")

ALB1$template="ALB1"
CX3$template="CX3"
NRCAM$template="NRCAM"
W601$template="W601"
plot_df = rbind(ALB1, CX3, NRCAM, W601)

panel_c = ggplot(plot_df, aes(x=dist+24, y=avg, group=group, color=template))+geom_line(show.legend = F, alpha=0.8)+
  facet_grid(cols=vars(factor(template, levels=c("W601", "NRCAM", "ALB1", "CX3"))))+
  scale_color_brewer(palette="Set2")+
  labs(x="50-bp center position (bp)", y="Sequence cyclizability")+
  theme_bw()+theme(strip.background =element_rect(fill="white"))+
  theme(text=element_text(size=7), title=element_text(size=7))


free(panel_a)+panel_b+panel_c+plot_layout(nrow=3, heights=c(2,1,1))

ggsave("SuppFig7.png", width=180, height=120, units = "mm")

