source("utils.R")
source("config.R")

complete_df = read.csv(SELEX_cyc_result_csv, row.names = 1)
inner_region = 35:66 # region without overlapping with motifs

filtered_df = complete_df %>% filter(dist %in% inner_region)
# Calculate cyclizability slope in the internal region
fitted_models = filtered_df %>% group_by(tf, group) %>%
  do(slope = lm(avg~dist, data = .)$coefficients[[2]])

# Calculate average
avgs = complete_df %>% filter(dist %in% inner_region) %>% 
  group_by(tf, group) %>% summarise(avg=mean(avg))

counts = complete_df %>% select(tf, group, n) %>% unique()

merged = inner_join(fitted_models, avgs, by=c("tf", "group"))
merged = inner_join(merged, counts, by=c("tf", "group"))
merged$slope = unlist(merged$slope)

write.csv(merged, SELEX_cyc_slope_avg)