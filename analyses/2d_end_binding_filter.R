source("utils.R")
source("config.R")

end_perc = read.csv(SELEX_end_binding_perc)

strong_end_binders = end_perc[end_perc$end_bind_perc>0.08,2:3]
non_end_binders = end_perc[end_perc$end_bind_perc<0.06,2:3]

write.csv(strong_end_binders, SELEX_end_binder_list)
write.csv(non_end_binders, SELEX_non_end_binder_list)