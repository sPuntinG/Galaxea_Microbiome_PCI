# ---------------------------- What's going on here ----------------------------------
#
# This script is the start of data structuring. Meaning:
#  - gathering all the necessary info in one place (tibble)
#  - shaping it in a manner that is handy to use later on
#    e.g., for data exploration and plotting and other transformations.
# 
#  Structuring continues in 02_remove_contaminants.R where the table is 
#  reshaped to long format. 
# 
# ---------------------------------------------------------------------------------


# Load packages --------------------------------------------

library(here)
library(tidyverse)
library(janitor)
library(ggforce)


# Import data and formatting --------------------------------------------

metadata <- read_csv("./in/metadata_ms.csv")

features_freq <- read_csv("./in/q2_outputs/table-filtered_feature-freq.csv")

# This needs reshaping
asv_count_taxa <- read_tsv("./in/q2_outputs/table_filtered_count_wtaxa.tsv")

# Drop the fist row that is the data type from python
asv_count_taxa <- asv_count_taxa %>% 
  slice(-1) %>% 
  rename(ASV_id = id)
  

# Extract ASV counts (similar to "OTU_counts" ) ....................
ASV_counts <- asv_count_taxa %>% 
  select(ASV_id, contains("-16S")) # select all cols with sample names/id

ASV_counts <- ASV_counts %>% 
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname, values_from = value)

# Make first row into columns names
ASV_counts <- ASV_counts %>% 
  janitor::row_to_names(row_number = 1) %>% 
  rename(sample_id = ASV_id)

# Make counts numeric!
ASV_counts <- ASV_counts %>% 
  mutate_at(vars(-("sample_id")), as.numeric)

# Make long
ASV_counts <- ASV_counts %>% 
  pivot_longer(
    -sample_id,
    names_to = "ASV_id",
    values_to = "count_sample" ) # *** former "count" *** 


# Extract taxonomy = ASV_id + taxonomy .............
taxonomy <- asv_count_taxa %>% 
  select(ASV_id, Sequence, Taxon)

# taxa: Split 'Taxon' into 7 columns (one for each taxonomic level) -----------------------
taxonomy <- taxonomy %>% 
  mutate(Taxon = stringr::str_replace_all(Taxon, "[:alpha:]{1}__", "")) %>% 
  separate(Taxon, sep = "; ", 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# Example/guide:
# d__Bacteria; 
# p__Proteobacteria; 
# c__Gammaproteobacteria; 
# o__Alteromonadales; 
# f__Alteromonadaceae; 
# g__Alteromonas; 
# s__Alteromonas_sp.

# taxa: replace NA with "Unclassified" 
taxonomy <- taxonomy %>% 
  mutate(across(('Kingdom':'Species'), ~replace_na(., "Unclassified")))


# RENAME "Unclassified" and "uncultured" based on lowest known taxonomic level ..............

# "Kingdom" > "Phylum" > "Class" > "Order" > "Family"
taxonomy <- taxonomy %>% 
  mutate(
    unclassif_FAM = if_else(Family != "Unclassified", Family, 
                            if_else(Order != "Unclassified", paste0("unclassif_", Order), 
                                    if_else(Class != "Unclassified", paste0("unclassif_", Class),
                                            if_else(Phylum != "Unclassified", paste0("unclassif_", Phylum), "unclassif_Bacteria"))))
  ) %>% 
  mutate(
    uncult_FAM = if_else(Family != "uncultured", Family, 
                         if_else(Order != "uncultured", paste0("uncult_", Order), 
                                 if_else(Class != "uncultured", paste0("uncult_", Class),
                                         if_else(Phylum != "uncultured", paste0("uncult_", Phylum), "uncult_Bacteria"))))
    
  ) %>% 
  mutate(
    uncX_FAM = if_else(unclassif_FAM == "uncultured", uncult_FAM, unclassif_FAM)
  ) %>% 
  relocate(uncX_FAM, .after = Family) %>% 
  select(-unclassif_FAM, -uncult_FAM)


# Rename "PeM15" coz is not a real name for a Family 
taxonomy <- taxonomy %>% 
  mutate(uncX_FAM = if_else(
    Family == "PeM15", "unclassif_Actinobacteria", uncX_FAM)) 


# write csv - to use later with {phyloseq} ------------------------------

# write_csv(taxonomy, "./out/Gfas_16S/taxonomy.csv")
write_csv(taxonomy, "./out/Gfas_16S/useful_tables/taxonomy_all_nonraref.csv")


# Join all data together ----------------------------------------------
all_data <- metadata %>% 
  select(-c(new_name, sample_type, origin)) %>% 
  inner_join(., ASV_counts, by = "sample_id") %>% #view()
  inner_join(., taxonomy, by = "ASV_id") # %>% view()



# Clean up env
rm(list = setdiff(ls(), c("all_data", "asv_count_taxa")))

# Next script: identify potential contaminants

