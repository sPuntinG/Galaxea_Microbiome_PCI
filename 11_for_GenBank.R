
# Prep fasta file for repersentative sequences (ASV) submission to NCBI GenBank

library(tidyverse)
library(here)

data <- read_csv("./out/Gfas_16S/useful_tables/ASV_numbered_forGenBank.csv")

lowest_taxlevel <- data %>% 
  mutate(
    unclassif_Genus = if_else(Genus != "Unclassified", Genus, 
                            if_else(Family != "Unclassified", paste0("unclassif_", Family), 
                                    if_else(Order != "Unclassified", paste0("unclassif_", Order), 
                                            if_else(Class != "Unclassified", paste0("unclassif_", Class),
                                                    if_else(Phylum != "Unclassified", paste0("unclassif_", Phylum), "unclassif_Bacteria")))))
  ) %>% 
  mutate(
    uncult_Genus = if_else(Genus != "uncultured", Genus,
                         if_else(Family != "uncultured", paste0("uncult_", Family), 
                                 if_else(Order != "uncultured", paste0("uncult_", Order), 
                                        if_else(Class != "uncultured", paste0("uncult_", Class),
                                                if_else(Phylum != "uncultured", paste0("uncult_", Phylum), "uncult_Bacteria")))))
    
  ) %>% 
  mutate(
    lowest = if_else(unclassif_Genus == "uncultured", uncult_Genus, unclassif_Genus)
  ) %>% 
  relocate(lowest, .after = Genus) #%>% view()


lowest_taxlevel %>% 
  select(ASV_nr, lowest, Sequence) %>% view()


genbank <- lowest_taxlevel %>% 
  select(ASV_nr, Sequence) %>% view()


# Make a fasta file (.txt) in tidyverse ---------------------------------
# (faster than figureing out other packages)

genbank %>% 
  mutate(forfasta = paste0(">", ASV_nr, "\n", Sequence, "\n\n")) %>% 
  # view() 
  # select(forfasta) %>% 
  pull(forfasta) %>% 
  str_flatten() %>% 
  write_file(., "./out/Gfas_16S/useful_tables/forGenBank_Gfas16S.txt")

