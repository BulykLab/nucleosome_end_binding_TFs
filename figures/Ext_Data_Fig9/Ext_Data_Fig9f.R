library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggrepel)
library(viridis)

plot_df = read.csv("../data/internal_cyc.csv")

end_binder = read.csv("../data/end_binding.csv")
ht_thres = read.csv("../data/ht_thres_TF.csv")

end_plot_df = plot_df %>% filter(tf %in% ht_thres$x) %>% filter(tf%in%end_binder$tf)

avg_end_plot_df = end_plot_df %>% group_by(group, dist) %>% summarise(avg=mean(avg)) %>% ungroup()

CEBPB_plot_df = plot_df %>% filter(tf == "CEBPB")

g1=ggplot(CEBPB_plot_df, aes(x=dist+24, y=avg, group=as.factor(group), color=as.factor(group)))+geom_line(linewidth=1)+
  theme_classic()+xlim(c(25, 125))+
  labs(x="50-bp center position (bp)", 
       y="Sequence cyclizability",
       color="TF binding motif\n position (bp)")+
  geom_rect(data=CEBPB_plot_df[1,], aes(xmin=25, xmax=49, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3, color=NA)+
  geom_rect(data=CEBPB_plot_df[1,], aes(xmin=122-22, xmax=123, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3, color=NA)+
  geom_rect(data=CEBPB_plot_df[1,], aes(xmin=25, xmax=74, ymin=-Inf, ymax=Inf), fill="red", alpha=0.2, color=NA)+
  geom_vline(xintercept = 74,linetype="dotted")+
  scale_color_viridis(discrete=T,begin=0.25, end=0.75,option="magma")+
  theme(text=element_text(size=7), title=element_text(size=7),
        legend.key.height = unit(2.5, "line"),
        legend.spacing.y = unit(0.5, 'cm'))+
  ggtitle("CEBPB")


g2=ggplot(avg_end_plot_df, aes(x=dist+24, y=avg, group=as.factor(group), color=as.factor(group)))+geom_line(linewidth=1)+
  theme_classic()+xlim(c(25, 125))+
  labs(x="50-bp center position (bp)", 
       y="Sequence cyclizability",
       color="TF binding motif\n position (bp)")+
  geom_rect(data=avg_end_plot_df[1,], aes(xmin=25, xmax=49, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3, color=NA)+
  geom_rect(data=avg_end_plot_df[1,], aes(xmin=122-22, xmax=123, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.3, color=NA)+
  geom_rect(data=avg_end_plot_df[1,], aes(xmin=25, xmax=74, ymin=-Inf, ymax=Inf), fill="red", alpha=0.2, color=NA)+
  geom_vline(xintercept = 74,linetype="dotted")+
  scale_color_viridis(discrete=T,begin=0.25, end=0.75,option="magma")+
  theme(text=element_text(size=7), title=element_text(size=7),
        legend.key.height = unit(2.5, "line"),
        legend.spacing.y = unit(0.5, 'cm'))+
  ggtitle("End Binders")

g1+g2+plot_layout(nrow=1, guides = "collect")
ggsave("SuppFig12.png", units = "mm", width=180, height=80)
