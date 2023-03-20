
source("00_data_structuring.R")


# Filter data to keep only ASV_id and sort by abundance (count per sample)
# Add 'count_overall' = features count across the dataset (count_sample is by ... sample)

find_contaminats <- all_data %>% 
  select(sample_id, state, colony_ms, levels_colony_ms, state_colony_frgmt_ms, ASV_id, count_sample) %>% 
  group_by(ASV_id) %>% 
  arrange(desc(count_sample)) %>% 
  mutate(count_overall = sum(count_sample)) %>% 
  ungroup() %>% 
  relocate(count_overall, .after = count_sample) 


# Add 'rank_overall' = ranking by relative abundance across the whole dataset 
find_contaminats <- find_contaminats %>% 
  select(ASV_id, count_overall) %>% 
  unique() %>% 
  arrange(desc(count_overall)) %>% 
  mutate(rank_overall = row_number()) %>% 
  right_join(., find_contaminats, by = c("ASV_id", "count_overall")) %>% 
  relocate(count_sample, .after = rank_overall) # %>%   view()


# identify the features found in C1 (control) with > 100 reads (should be 9) ............................
C1_100 <- find_contaminats %>% 
  filter(sample_id == "C1-16S" & count_sample > 100) %>% 
  pull(ASV_id)

# Use C1_100 to filter sequences to plot 
# (to viz assess potential contaminates based on presence/absence between C1 and coral samples)
plot_find_contaminats <- find_contaminats %>% 
  filter(ASV_id %in% C1_100) # %>%  filter(count > 0) 


# PLOT  --------------------------------------------
ggplot(data = plot_find_contaminats, 
            aes(x = forcats::fct_reorder(state_colony_frgmt_ms, levels_colony_ms ), 
                y = as_factor(rank_overall),
                size = count_sample)) +
  geom_point(aes(color = colony_ms),
             shape = 16, 
             stroke = 0) + 
  scale_size_area(name = "seqs count \nper sample",
                  max_size = 10) +
  ggtitle("All features found in the PCR control with >100 reads") + 
  ylab("Feature ranking by overall frequency\n(1 = most frequent across the whole data set)") +
  theme_bw() +
  theme(
    text = element_text(colour = "white"),
    axis.text = element_text(color = "white"),
    plot.background = element_rect(fill = "black", color = "black"),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.key = element_rect(fill = "grey21"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) + 
  ggforce::facet_row(vars(forcats::fct_rev(factor(state))), 
                       scales = 'free_x', space = 'free')

ggsave("./out/Gfas_16S/find_contaminants/find_contaminants.png", 
       dpi = 330, 
       units = "cm", width = 25, height = 13)



# Make table to export  --------------------------------------------
tab1 <- find_contaminats %>% 
  filter(sample_id == "C1-16S") %>%
  arrange(desc(count_sample)) %>% 
  mutate(
    rank_ctrl = row_number(),
    count_ctrl = count_sample,
  ) %>% 
  # rename(ASV_id = taxon_name) %>% 
  select(ASV_id, count_overall, rank_overall, count_ctrl, rank_ctrl) #%>% view()

tab2 <- all_data %>% 
  filter(ASV_id %in% C1_100) %>% 
  select(ASV_id, c(Kingdom:Species)) %>% 
  unique()

tab3 <- asv_count_taxa %>% 
  select(ASV_id, Sequence) # to add the sequences!


potential_contaminants <- inner_join(tab1, tab2, by = "ASV_id") %>% 
  inner_join(., tab3, by = "ASV_id")

# Export table
# This will be used to BLASTn search and identify likely contaminants
write_csv(potential_contaminants, 
          "./out/Gfas_16S/find_contaminants/potential_contaminants.csv")


# Clean up 
rm(list = setdiff(ls(), c("all_data", "asv_count_taxa")))


# Now I will open this in excel and include the info that was manually searched (BLASTn)
# then flag as to keep or remove and import again and remove contaminant sequences

