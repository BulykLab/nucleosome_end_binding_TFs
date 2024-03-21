end_perc = read.csv("../data/end_binding_percent.csv")

strong_end_binders = end_perc[end_perc$end_bind_perc>0.08,2:3]
non_end_binders = end_perc[end_perc$end_bind_perc<0.06,2:3]

write.csv(strong_end_binders,"../data/strong_end_binding.csv")
write.csv(non_end_binders,"../data/non_end_binding.csv")