library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggrepel)

# Panel a-b

tf_list = read.csv('../data/tf_DBD_class.csv')
ht_thres_tfs = read.csv('../data/ht_thres_TF.csv')
end_binders = read.csv('../data/end_binding.csv')
non_binders = read.csv('../data/non_end_binding.csv')

end_binders_pie=tf_list %>% filter(TF %in% ht_thres_tfs$x) %>% filter(TF %in% end_binders$tf) %>% group_by(DBD.Class) %>% summarise(count=n())
non_binders_pie=tf_list %>% filter(TF %in% ht_thres_tfs$x) %>% filter(TF %in% non_binders$tf) %>% group_by(DBD.Class) %>% summarise(count=n())


end_binders_pie$percent = end_binders_pie$count/sum(end_binders_pie$count)
non_binders_pie$percent = non_binders_pie$count/sum(non_binders_pie$count)

end_binders_pie$class = "End binders"
non_binders_pie$class = "Non-end binders"


panel_ab=ggplot(rbind(end_binders_pie, non_binders_pie),aes(x = factor(1),y=percent, fill=DBD.Class)) + 
  facet_wrap(~class) + 
  geom_bar(width = 1, stat = "identity") + 
  coord_polar(theta="y")+
  labs(fill="TF DBD class")+
  geom_text(aes(label = count), size=2,
            position = position_stack(vjust = 0.5),
            fontface="bold")+theme_void()+
  theme(text=element_text(size=7), title=element_text(size=7))+
  guides(fill=guide_legend(ncol=2))

# Panel c-d
slope_avg_plot_df = read.csv("../data/SELEX_cyc_slope_avg.csv")
labeling_TF = c("CEBPB")
slope_avg_plot_df$label = ifelse(slope_avg_plot_df$tf %in% labeling_TF,slope_avg_plot_df$tf,"" )

slope_avg_plot_df = slope_avg_plot_df %>% filter(tf %in% ht_thres_tfs$x)%>% filter(tf %in% non_binders$tf | tf %in% end_binders$tf)
slope_avg_plot_df$end = ifelse(slope_avg_plot_df$tf %in% end_binders$tf, "End binders (70 TFs)", "Non-end binders (29 TFs)")
panel_cd=ggplot(slope_avg_plot_df, aes(x=avg, y=slope, group=group, color=group, shape=group, label=label))+
  geom_point(size=1.5)+theme_bw()+
  scale_color_manual(values=c("#00A651", "#7F3F98", "#f4b553"),
                     breaks=c('Nucleosome-SELEX', 
                              'NCAP-SELEX', 
                              'HT-SELEX'))+
  scale_shape_manual(values=c(16,17,18),
                     breaks=c('Nucleosome-SELEX', 
                              'NCAP-SELEX', 
                              'HT-SELEX'))+
  labs(x = "Average internal sequence cyclizability", y="Internal sequence cyclizability slope",
       color="Experiment", shape="Experiment")+
  geom_text_repel(size=2,box.padding = 0.5, max.overlaps = Inf, show.legend = F, color="black")+
  facet_grid(cols=vars(end))+
  theme(text=element_text(size=7), title=element_text(size=7))


# Panel e
all_left_fft = read.csv("../data/SELEX_fft_A_35_74.csv")
all_right_fft = read.csv("../data/SELEX_fft_A_76_115.csv")

all_left_fft = all_left_fft %>% filter(TF %in% ht_thres_tfs$x) %>% filter(TF %in% end_binders$tf)
all_right_fft = all_right_fft %>% filter(TF %in% ht_thres_tfs$x) %>% filter(TF %in% end_binders$tf)

average_left = all_left_fft %>% group_by(x, library) %>% summarize(avg=mean(spectrum.density))
average_left$TF = "Mean"
average_left$pos = "Region 35-74"

average_right = all_right_fft %>% group_by(x, library) %>% summarize(avg=mean(spectrum.density))
average_right$TF = "Mean"
average_right$pos = "Region 76-115"

rect_df = rbind(all_left_fft, all_right_fft) %>% group_by(library, pos) %>% summarise(spectrum.density=mean(spectrum.density))
rect_df$TF = "rect"
rect_df$x=0

panel_e=ggplot(rbind(all_left_fft, all_right_fft) , aes(x=x, y=spectrum.density, group=TF, color=library))+
  facet_grid(cols=vars(library), rows=vars(pos),switch="y")+
  geom_line(show.legend = F)+
  geom_rect(data=rect_df, aes(xmin=0, xmax=0.6, ymin=0.15, ymax=0.75), alpha=0.8, color=NA, fill="white")+
  geom_line(data=rbind(average_left, average_right), aes(x=x, y=avg), linewidth=0.8,alpha=1,show.legend = F)+
  labs(x="Sequence spatial frequency of nucleotide A (bp-1)",
       y="Periodicity density")+
  scale_color_manual(values=c("#00A651", "#7F3F98", "#f4b553"),
                     breaks=c('Nucleosome-SELEX', 
                              'NCAP-SELEX', 
                              'HT-SELEX'))+theme_classic()+
  geom_vline(xintercept = 0.1, linewidth=0.2, linetype="dashed")+
  theme(text=element_text(size=7), title=element_text(size=7))

# Panel f
all_left_fft = read.csv("../data/SELEX_fft_A_35_74.csv")
all_right_fft = read.csv("../data/SELEX_fft_A_76_115.csv")

all_left_fft = all_left_fft %>% filter(TF %in% ht_thres_tfs$x) %>% filter(TF %in% non_binders$tf)
all_right_fft = all_right_fft %>% filter(TF %in% ht_thres_tfs$x) %>% filter(TF %in% non_binders$tf)

average_left = all_left_fft %>% group_by(x, library) %>% summarize(avg=mean(spectrum.density))
average_left$TF = "Mean"
average_left$pos = "Region 35-74"

average_right = all_right_fft %>% group_by(x, library) %>% summarize(avg=mean(spectrum.density))
average_right$TF = "Mean"
average_right$pos = "Region 76-115"

rect_df = rbind(all_left_fft, all_right_fft) %>% group_by(library, pos) %>% summarise(spectrum.density=mean(spectrum.density))
rect_df$TF = "rect"
rect_df$x=0

panel_f=ggplot(rbind(all_left_fft, all_right_fft), aes(x=x, y=spectrum.density, group=TF, color=library))+
  facet_grid(cols=vars(library), rows=vars(pos),switch="y")+
  geom_line(show.legend = F)+
  geom_rect(data=rect_df, aes(xmin=0, xmax=0.6, ymin=0.15, ymax=0.75), alpha=0.8, color=NA, fill="white")+
  geom_line(data=rbind(average_left, average_right), aes(x=x, y=avg), linewidth=0.8,alpha=1,show.legend = F)+
  labs(x="Sequence spatial frequency of nucleotide A (bp-1)",
       y="Periodicity density")+
  scale_color_manual(values=c("#00A651", "#7F3F98", "#f4b553"),
                     breaks=c('Nucleosome-SELEX', 
                              'NCAP-SELEX', 
                              'HT-SELEX'))+theme_classic()+
  geom_vline(xintercept = 0.1, linewidth=0.2, linetype="dashed")+
  theme(text=element_text(size=7), title=element_text(size=7))

design <- "
  11
  22
  34
"
free(panel_ab)+free(panel_cd)+panel_e+panel_f+plot_layout(design=design, heights = c(1.5,1.5,1))
ggsave("SuppFig11.png", width=180, height=170, units = "mm")
