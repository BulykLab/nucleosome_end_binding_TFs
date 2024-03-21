library(dplyr)

complete_df = read.csv("../data/cyc_results.csv", row.names = 1)
inner_region = 12:41

filtered_df = complete_df %>% filter(dist %in% inner_region)

fitted_model = filtered_df %>% group_by(tf, group) %>% do(slope = lm(avg~dist, data = .)$coefficients[[2]])

avgs = complete_df %>% filter(dist %in% inner_region) %>% 
  group_by(tf, group) %>% summarise(avg=mean(avg))

counts = complete_df %>% select(tf, group, n) %>% unique()

merged = inner_join(fitted_models, avgs, by=c("tf", "group"))
merged = inner_join(merged, counts, by=c("tf", "group"))
merged$slope = unlist(merged$slope)

write.csv(merged, "data/all_tfs_slope_avg.csv")