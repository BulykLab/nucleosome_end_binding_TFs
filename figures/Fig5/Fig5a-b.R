library(ggplot2)

ex_chr2_cyc = read.table("../data/C0free_prediction_chr2_locus.txt")
ex_chr2_cyc$dist = 1:nrow(ex_chr2_cyc)

ggplot(ex_chr2_cyc, aes(x=dist, y=V2))+geom_line(color="#0000B2")+theme_classic()+labs(y="Cyclizability", x="")+guides(x = "none")

ggsave("5a_example_locus_chr2_cyc.png", width=170, height=50, units="mm")


sample = get_cycs("../data/CEBPB_k562_data/head_sample.csv",
         "../data/CEBPB_k562_data/tail_sample.csv")
control = get_cycs("../data/CEBPB_k562_data/head_control_clean.csv",
                   "../data/CEBPB_k562_data/tail_control_clean.csv")

sample_plot_df = build_plot_df(t(sample),"CEBPB", paste0("CEBPB-bound (n=",nrow(sample),")"))
control_plot_df = build_plot_df(t(control), "CEBPB", paste0("CEBPB-unbound (n=",nrow(control),")"))

library(ggplot2)
ggplot(rbind(sample_plot_df, control_plot_df), aes(x=dist+24, y=avg, ymin=low, ymax=high, color=group))+
  geom_line()+geom_ribbon(aes(fill=group), color=NA, alpha=0.1)+theme_classic()+
  scale_color_manual(values=c("#FF938B","#0000B2"))+
  scale_fill_manual(values=c("#FF938B","#0000B2"))+
  geom_rect(data=sample_plot_df[1,], aes(xmin=25, xmax=54, ymin=-Inf, ymax=Inf), fill="#C1272D", alpha=0.2, color=NA)+
  geom_vline(xintercept = 70, linetype="dotted")+
  labs(x="50-bp center position (bp)", y="Sequence cyclizability", color="K562 nucleosomes", fill="K562 nucleosomes")#+
  #theme(legend.position = c(0.75,0.3))

ggplot(sample_plot_df, aes(x=dist, y=avg, ymin=low, ymax=high, color=group))+geom_line()+geom_ribbon(aes(fill=group), color=NA, alpha=0.1)
ggsave("5b.png", width=100, height=100, units="mm")
