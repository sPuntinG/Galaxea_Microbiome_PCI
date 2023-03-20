
source("02_remove_contaminants.R")

metadata <- read_csv("./in/metadata_ms.csv")

# Check how many reads in C1 -----------------------------
decont_long %>% 
  filter(taxon_level == "ASV_id") %>% 
  filter(count_sample != 0) %>% 
  group_by(sample_id) %>% 
  summarise(
    nr_reads = sum(count_sample),
    nr_features = length(unique(taxon_name)) 
  ) %>% #view()
  summarise(
    nr_samples = length(sample_id),
    tot_nr_reads = sum(nr_reads),
    tot_nr_features = sum(unique(nr_features))
  )

# nr_samples | tot_nr_reads | tot_nr_features
#    29      |    114902    |     1136


# Make summary table for SM: group, reads, features ---------------------
decont_long %>% 
  filter(taxon_level == "ASV_id") %>% 
  filter(count_sample != 0) %>% 
  group_by(state, colony_ms) %>% 
  summarise(
    nr_reads = sum(count_sample),
    nr_features = length(unique(taxon_name)),
    n = length(unique(sample_id)) ) %>%  #,
  arrange(desc(state), colony_ms) %>% 
  drop_na()  %>% #view()
  write_csv(., "./out/Gfas_16S/useful_tables/group_reads_features.csv")


# Make summary table for SM: sample, reads, features ---------------------
decont_long %>% 
  filter(taxon_level == "ASV_id") %>% 
  filter(count_sample != 0) %>% 
  group_by(sample_id) %>%
  summarise(
    nr_reads = sum(count_sample),
    nr_features = length(unique(taxon_name)) ) %>%  
  left_join(., metadata, by = "sample_id") %>% 
  select(sample_id, state_colony_frgmt_ms, state, colony_ms, nr_reads, nr_features) %>% 
  arrange(desc(state), colony_ms) %>%
  drop_na() %>%  
  write_csv(., "./out/Gfas_16S/useful_tables/sample_reads_features.csv")



# Remove control sample (C1) --------------------------------------------
decont_long <- decont_long %>% 
  filter(sample_id != "C1-16S")


# # Export RDS (keep 0s in count for presence/absence)  -------------------------------
# saveRDS(decont_long, "./out/RDS_files/decont_long_0s.rds")


# Remove zeroes -------------------------------------
decont_long <- decont_long %>% 
  filter(count_sample != 0)


# Add grouping (state x origin x colony) ----------------------------------------------
decont_long <- decont_long %>% 
  mutate(group = paste(state, colony_ms, sep = "_")) %>% 
  relocate(group, .after = colony_ms)



# Export RDS -------------------------------
saveRDS(decont_long, "./out/RDS_files/decont_long.rds")

