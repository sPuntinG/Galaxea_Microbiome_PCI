# ---------------------------------------------------------------------------------
# 
# In this script: I REMOVE the CONTAMINANTS
# which are 5 features (ASVs) that we have identified as 
# likely to be contaminants based on presence/absence and 
#  BLSTn matches (info about source or environment)
# 
# ---------------------------------------------------------------------------------

source("01_find_contaminants.R")


# Import the table that I have BLASTn-ed and annotated ------------------
contaminants <- read_csv("./out/Gfas_16S/find_contaminants/potential_contaminants_GP_MZ_GP_simplified_for_R.csv")

# Clean it up to have a list of the sequences to remove
contaminants <- contaminants %>% 
  select(Feature, Sequence, Frequency, Final_decision) %>% 
  drop_na() %>% 
  filter(Final_decision == "remove") %>% 
  select(-Final_decision) %>% 
  rename(ASV_id = Feature,
         count_overall = Frequency)


# Write .csv file for import in Qiime2
write_csv(contaminants, 
          "./out/Gfas_16S/find_contaminants/contaminants_features_sequences.csv")



# Remove CONTAMINANTS ----------------------------------------

decont <- contaminants %>% 
  select(ASV_id) %>% 
  anti_join(all_data, ., by = "ASV_id")


# Check how many features have been removed (ASV_id) - should be 5  .....................
all_data$ASV_id %>% unique() %>% length()     # 526 -
contaminants$ASV_id %>% unique() %>% length() #   5 =
decont$ASV_id %>% unique() %>% length()       # 521     :D

# Check how many reads have been removed ...............................................
all_data$count_sample %>% sum()      # 133892 -
contaminants$count_overall %>% sum() #  18990 = 
decont$count_sample %>% sum()        # 114902

sum(all_data$count_sample) - sum(contaminants$count_overall) == sum(decont$count_sample)
# # TRUE   all good :D

# So now we have removed the contaminants and we can proceed

# Let's clean up a bit of memory ... keep only 'decont'
rm(list = setdiff(ls(), c("decont", "metadata")))  # 


# Add RELATIVE ABUNDANCE -------------------------------------------
decont <- decont %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = count_sample / sum(count_sample)) %>% 
  relocate(rel_abund, .after = count_sample) %>% 
  ungroup()

# Check that rel_abund sums up to a tot of 1
decont %>% 
  group_by(sample_id) %>% 
  summarise(tot_abund = sum(rel_abund)) #%>% view()
# all good :D


# Make LONG -------------------------------------------------
decont_long <- decont %>% 
  pivot_longer(
    cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "uncX_FAM",
             "Genus", "Species", "ASV_id"),
    names_to = "taxon_level",
    values_to = "taxon_name"
  )

rm(decont)

# Ok, so now the data is ready to be plotted! 
# just source this code in the next scripts ;)