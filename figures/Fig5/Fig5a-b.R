library(ggplot2)

ex_chr2_cyc = read.table("../data/C0free_prediction_chr2_locus.txt")
ex_chr2_cyc$dist = 1:nrow(ex_chr2_cyc)

ggplot(ex_chr2_cyc, aes(x=dist, y=V2))+geom_line(color="#0000B2")+theme_classic()+labs(y="Cyclizability", x="")+guides(x = "none")

ggsave("example_locus_chr2_cyc.png", width=200, height=50, units="mm")
