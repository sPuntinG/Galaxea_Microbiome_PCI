
source("03_check_data_structuring2.R")

############################################################
# Identify the most abundant FAMILIES across the data set ##
############################################################


# To pick the most frequent FAMILIES (not ASVs)
# across the whole data set (not sample-wise)

# Calculate the TOTAL NR of READS in the data set ...................................
tot_reads <- decont_long %>% 
  filter(taxon_level == "ASV_id") %>% 
  group_by(taxon_name) %>% 
  pull(count_sample) %>% sum()
# 112789

# 112'789 total nr of reads in the data set as of now 
# (= after removing contaminants and control C1)
# which is = 114902 - 2113 = tot - C1

# Calculate relative abundance of each FAM (across the data set) ......................
# = tot nr of counts for each family (how many reads are assigned to each family)
# divided by tot nr of reads in the data set
tot_abundances_byFAMILY <- decont_long %>% 
  # filter(count != 0) %>%  # remove zeroes
  filter(taxon_level == "uncX_FAM") %>% # keep only formatted family class
  group_by(taxon_name) %>%              # group by unique Family
  summarise(
    tot_counts = sum(count_sample),
    tot_rel_abund = sum(count_sample) / tot_reads) %>% 
  ungroup() %>% 
  arrange(desc(tot_counts))  #%>% #view()
# pull(tot_rel_abund) %>% sum() # 1 = sanity check passed ;)
# pull(tot_counts) %>% sum()      # 112789 = sanity check passed ;)
# view()

# This above is a reference table that summarized the number of reads by Family 
# across the whole data set (ALL samples pooled together)
# Note that this is NOT the count/abundance of the most common ASV's, BUT
# the count of the ASVs that are classified under the same Family
# (= each Family contains >=1 unique ASV_id!)
# (= different ASVs are pooled together by Family)

write_csv(tot_abundances_byFAMILY, "./out/Gfas_16S/useful_tables/tot_abundances_byFAMILY.csv")


# Identify the most representative families = tot_rel_abund >= 0.01 (1%) .............
top1perc <- tot_abundances_byFAMILY %>% 
  filter(tot_rel_abund >= 0.01)

top1perc$tot_rel_abund %>% sum() # 0.8369522 = ~ 83.7 %

# This results in 17 Families
# which account for a total of 83.7 % of all reads in the data set

# Extract their names
top_FAMs <- top1perc$taxon_name
length(top_FAMs) # 17


# Keep only FAMILY data ---------------------------------------------------
fam <- decont_long %>% 
  filter(taxon_level == "uncX_FAM") 

# Recalculate relative abundance by Family  ------------------------------------
# ... because the current counts and rel abundances are by ASV.
# But since we pooled diff ASV together by family, we need to update the counts
# => group by sample_id and family name ("taxon_name") and sum counts and rel_abund

# by sample
fam <- fam %>% 
  group_by(sample_id, taxon_name) %>%              # group by unique Family
  mutate(
    count_bysample_decont = sum(count_sample),
    rel_abund_bysample_decont = sum(rel_abund) ) %>% # rel_abund_decont = sum(rel_abund) )
  ungroup() %>% 
  select(-count_sample, -rel_abund, -Sequence) %>% # necessary to be able to remove duplicated rows
  unique()                       # remove duplicated rows

# Sanity check: check that 
# - tot counts (tot nr of reads) is the same as 'tot_reads' = 112'789
# - tot_abund_bysample = 1 (100%)
fam %>% 
  group_by(sample_id) %>%
  summarise(
    tot_counts = sum(count_bysample_decont),
    tot_abund_bysample = sum(rel_abund_bysample_decont)
  )  %>% 
 view()  # see that each sample tot_abund_bysample = 1 :D
# pull(tot_counts) %>% sum() # 112789  correct!           :D


# Pool "Others" together -----------------------------------------

# Current number of FAMs
fam$taxon_name %>% unique() %>% length() # 138

# Pool all low-abundance Families under "Others"
fam <- fam %>% 
  mutate(family = if_else(taxon_name %in% top_FAMs, taxon_name, "Others")) 


# Recalculate rel_abundance - needed?
fam <- fam %>% 
  group_by(sample_id, family) %>% 
  mutate(
    count_bysample = sum(count_bysample_decont), 
    rel_abund_bysample = sum(rel_abund_bysample_decont)
  ) %>% 
  ungroup() %>% 
  # view() %>% 
  select(-count_bysample_decont, -rel_abund_bysample_decont, -taxon_name) %>% 
  unique() # %>% view()

# Sanity check again .....
fam %>% 
  # summarise(tot_count = sum(count)) # 112789
  group_by(sample_id) %>% 
  summarise(rel_abund = sum(rel_abund_bysample)) # %>% view()
# passed :D



# Sanity check: How much "Others" accounts for? 
fam %>% filter(family == "Others") %>%  # view()
  pull(rel_abund_bysample) %>% mean()  # 0.1673156 makes sense (1 - 0.8369 ...)


# Clean up 
rm(list = setdiff(ls(), c("fam", "metadata")))
