library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(viridis)
library(rstatix)
# Panel a place holder

panel_a = plot_spacer()


# Panel b
CEBPB_cyc_results = read.csv("../data/CEBPB_plot_df.csv")
panel_b = ggplot(CEBPB_cyc_results,aes(x=dist+25, y=avg, color=group))+geom_line(linewidth=0.8)+
  theme_classic()+geom_ribbon(aes(ymin = low, ymax = high,fill=group), alpha = 0.2,  color=NA)+
  labs(x = "50-bp center position (bp)", 
       y="Sequence cyclizability (AU)", 
       fill="", 
       color="")+
  scale_color_manual(values=c("#00A651", "#7F3F98", "#f4b553"),
                     breaks=c('Nucleosome-SELEX', 
                              'NCAP-SELEX', 
                              'HT-SELEX'))+
  scale_fill_manual(values=c("#00A651", "#7F3F98", "#f4b553"),
                    breaks=c('Nucleosome-SELEX', 
                             'NCAP-SELEX', 
                             'HT-SELEX'))+
  theme(text=element_text(size=7), 
        legend.title = element_text( size=7),
        legend.text = element_text(size=7), 
        legend.spacing.y = unit(1.0, 'cm'),
        legend.margin=margin(t = 0, unit='cm'))+
  guides(fill = guide_legend(byrow = TRUE), color = guide_legend(byrow = TRUE))+
  geom_rect(data=CEBPB_cyc_results[1,], aes(xmin=25, xmax=25+23, ymin=-Inf, ymax=Inf),
            fill="gray", alpha=0.3, color=NA)+
  geom_rect(data=CEBPB_cyc_results[1,], aes(xmin=122-22, xmax=122, ymin=-Inf, ymax=Inf), 
            fill="gray", alpha=0.3, color=NA)+
  geom_rect(data=CEBPB_cyc_results[1,], aes(xmin=25, xmax=49+10, ymin=-Inf, ymax=Inf), 
            fill="#C1272D", alpha=0.2, color=NA)+
  geom_vline(xintercept = 74, linetype="dotted")



# Panel c
left_fft = read.csv("../data/SELEX_fft_A_35_74.csv")
right_fft = read.csv("../data/SELEX_fft_A_76_115.csv")
left_CEBPB_fft = left_fft %>% filter(TF=="CEBPB")
right_CEBPB_fft = right_fft %>% filter(TF=="CEBPB")
panel_c = ggplot(rbind(left_CEBPB_fft, right_CEBPB_fft), aes(x=x, y=spectrum.density, color=library))+
  geom_line(show.legend = F, linewidth=1)+facet_grid(cols=vars(pos))+theme_classic()+
  labs(x=bquote("Sequence spatial frequency of nucleotide A "(bp^-1)),
       y="Periodicity density (AU)")+
  geom_rect(data=rbind(left_fft[1,], right_fft[1,]), aes(xmin=0, xmax=0.5, ymin=0.2, ymax=0.5), color=NA, fill="white", alpha=0.5)+
  geom_line(data=rbind(left_CEBPB_fft, right_CEBPB_fft) %>% filter(library!="Nucleosome-SELEX"), aes(x=x, y=spectrum.density), show.legend = F,linewidth=1)+
  scale_color_manual(values=c("#00A651", "#7F3F98", "#f4b553"),
                     breaks=c('Nucleosome-SELEX', 
                              'NCAP-SELEX', 
                              'HT-SELEX'))+
  scale_fill_manual(values=c("#00A651", "#7F3F98", "#f4b553"),
                    breaks=c('Nucleosome-SELEX', 
                             'NCAP-SELEX', 
                             'HT-SELEX'))+
  geom_vline(xintercept = 0.1, linewidth=0.5, linetype="dashed")+
  theme(text=element_text(size=7), title=element_text(size=7))


# Panel d

cyc_results = read.csv("../data/SELEX_cyc_results.csv")
end_binder = read.csv("../data/end_binding.csv")$tf
non_end_binder = read.csv("../data/non_end_binding.csv")$tf
ht_thres_tf = read.csv("../data/ht_thres_TF.csv")$x

cyc_plot_df = cyc_results %>% filter(tf %in% ht_thres_tf) %>% filter((tf %in% end_binder) | (tf %in% non_end_binder))
cyc_plot_df$end_binding = ifelse(cyc_plot_df$tf %in% end_binder, "End binders", "Non-end binders")
avg_cycs = cyc_plot_df %>% group_by(end_binding, dist) %>% summarise(avg = mean(avg))

panel_d = ggplot(cyc_plot_df, aes(x=dist+24, y=avg, group=tf))+
  geom_line(color="#7F3F98")+
  geom_rect(data=cyc_plot_df[1,], aes(xmin=24, xmax=123, ymin=-Inf, ymax=Inf), fill="white", alpha=0.7)+
  geom_rect(data=cyc_plot_df[1,], aes(xmin=24, xmax=49, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3)+
  geom_rect(data=cyc_plot_df[1,], aes(xmin=122-22, xmax=123, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3)+
  geom_rect(data=cyc_plot_df[1,], aes(xmin=24, xmax=49+10, ymin=-Inf, ymax=Inf), fill="red", alpha=0.2)+
  geom_rect(data=cyc_plot_df[which(cyc_plot_df$end_binding=="End binders"),][1,], 
            aes(xmin=24, xmax=123, ymin=-Inf, ymax=Inf), fill="white", alpha=0.7)+  
  geom_rect(data=cyc_plot_df[which(cyc_plot_df$end_binding=="End binders"),][1,], 
            aes(xmin=24, xmax=49, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3)+
  geom_rect(data=cyc_plot_df[which(cyc_plot_df$end_binding=="End binders"),][1,], 
            aes(xmin=122-22, xmax=123, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3)+
  geom_rect(data=cyc_plot_df[which(cyc_plot_df$end_binding=="End binders"),][1,], 
            aes(xmin=24, xmax=49+10, ymin=-Inf, ymax=Inf), fill="red", alpha=0.2)+
  geom_vline(xintercept = 74,linetype="dotted")+
  theme_classic()+
  geom_line(data=avg_cycs, aes(x=dist+24, y=avg, group=end_binding), color="black")+
  facet_grid(cols=vars(end_binding))+
  labs(x="50-bp center position (bp)", y="Sequence cyclizability")+
  theme(text=element_text(size=7), title=element_text(size=7))


# Panel e
slope_avg = read.csv("../data/SELEX_cyc_slope_avg.csv") %>% filter(grepl("NCAP", group))
slope_avg  = slope_avg %>% filter(tf %in% ht_thres_tf) %>% filter((tf %in% end_binder) | (tf %in% non_end_binder))
slope_avg$end_binding = ifelse(slope_avg$tf %in% end_binder, "End\nbinders", "Non-end\nbinders")
panel_e = ggplot(slope_avg[order(slope_avg$end_binding),], aes(x=end_binding, 
                    y=slope))+geom_boxplot(color="#7F3F98")+ geom_jitter(color="#7F3F98",size=0.5)+
  geom_bracket(
    xmin = "End\nbinders", xmax = "Non-end\nbinders", y.position = 0.004,
    label = "***"
  )+
  labs(
    x="",
    y="Internal cyclizability slope (AU)"
  )+ylim(c(-0.0005, 0.0045))+theme_classic()+
  theme(title=element_text(size=7), text = element_text(size=7))


# Panel f
internal_cyc = read.csv("../data/SELEX_internal_cyc_results.csv")
end_binder_internal = internal_cyc %>% filter(tf %in% ht_thres_tf) %>% filter(tf%in%end_binder)
avg_height = end_binder_internal %>% filter(dist > 49) %>% group_by(tf, group) %>% summarize(height=mean(avg)) %>% ungroup()
avg_height$group = factor(avg_height$group)
panel_f=ggplot(avg_height, aes(x=group, y=height,  group=group, color=group,fill=group))+geom_boxplot(alpha=0.5)+
  scale_color_viridis(discrete = TRUE, limits = paste0(seq(25,44, 5)), begin=0.25, end=0.75,option="magma")+
  scale_fill_viridis(discrete = TRUE, limits = paste0(seq(25,44, 5)), begin=0.25, end=0.75,option="magma")+
  labs(x="TF binding position (bp)", y="Average internal\nsequence cyclizability (AU)", 
       color="", 
       fill="")+
  guides(fill = "none")+theme(text=element_text(size=9), title=element_text(size=9))+
  theme_classic()+ theme(text=element_text(size=7), title=element_text(size=7),
                         legend.key.height = unit(2.0, "line"),
                         legend.spacing.y = unit(0.5, 'cm'))+
  annotate("text", x=3, y=0.005, label=bquote("p = 1.08 x 10"^-27),size=2)

res.aov <- anova_test(data = avg_height, dv = height, wid = tf, within = group)
get_anova_table(res.aov)

panel_f_schema_spacer = plot_spacer()

design <- "
  111112229
  444433333
  555567778
"



panel_a + panel_b + panel_c + plot_spacer() + 
  panel_d + panel_e + panel_f+plot_spacer()+plot_spacer()+ plot_layout(design = design,
                                                                       widths = c(1,1,1,1,2,0.5,1,1,0.1))

ggsave("Fig4.png", width=180, height=150, units = "mm")
