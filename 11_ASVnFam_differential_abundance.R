######################################################
#
#   Differential abundance of ASV between symbiotic states 
#   (symb vs bleached)
#   on RAREFIED data
# 
#   Spoiler alert: nothing is significant here.
#   No clear "pattern" nor signature of bleached
#
######################################################

library(tidyverse)
library(here)


ps_rarefied.1000avg.rds <- readRDS("./out/RDS_files/ps_rarefied.1000avg.rds")
taxonomy_seqs <- read_csv("./out/Gfas_16S/useful_tables/taxonomy_all_nonraref.csv") 
metadata_ms <- read_csv("./in/metadata_ms.csv")

source("colors_GfasMS.R")


# Prep data - ASV level ------------------------------------------------------------------

# Get ASV table out of ps_rarefied  
ASV_table <- as(otu_table(ps_rarefied.1000avg.rds), "matrix") %>%  
  t() %>%
  as_tibble(., rownames = "sample_id")

# Rename samples to match ms format/naming
comp <- metadata_ms %>% 
  select(sample_id, state, origin_ms, colony_ms, state_colony_frgmt_ms, levels_colony_ms) %>% 
  inner_join(., ASV_table, by = "sample_id") %>% 
  select(-sample_id) %>% 
  mutate(levels_colony_ms = factor(levels_colony_ms))

# Make long 
comp <- comp %>% pivot_longer(., cols = where(is.numeric), 
                               names_to = "ASV", 
                               values_to = "abundance") 

# Calculate relative abundance of ASV by sample ('rel_abund')
comp <- comp %>% 
  group_by(state_colony_frgmt_ms) %>% 
  mutate(rel_abund = abundance / sum(abundance)) %>% 
  ungroup() 



# Pairwise comparisons of ASV rel.ab. - OVERALL symb vs bleached ----------------

# For each ASV, run pairwise comparison (Wilcoxon text) of 
# its rel.ab btw symbiotic and bleached samples
signif_ASV <- comp %>% 
  nest(data = -ASV) %>% 
  mutate(test = map(.x = data, 
                    ~wilcox.test(rel_abund~state, data = .x, exact = F) )) %>% 
  mutate(test = map(.x = test,
                     ~broom::tidy(.x))) %>% 
  unnest(test) %>% 
  # filter(p.value < 0.05) # only 14 ASVs ... # toogle this line
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p_adj < 0.05)

# Results: 
# - 14 ASVs signif. diff. between states with uncorrected Wilcoxon test
# - 0 after correction for Benjamini-Hochman ...

# Next steps: Repeat this but for Red Sea and Hong Kong separately



# Pairwise comparisons of ASV rel.ab. - RED SEA symb vs bleached ----------------
signif_ASV_RedSea <- comp %>% 
  filter(origin_ms == "Red Sea") %>% 
  nest(data = -ASV) %>% 
  mutate(test = map(.x = data, 
                    ~wilcox.test(rel_abund~state, data = .x, exact = F) )) %>% 
  mutate(test = map(.x = test,
                    ~broom::tidy(.x))) %>% 
  unnest(test) %>% 
  # filter(p.value < 0.05) # only 2 ASVs ... # toogle this line
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p_adj < 0.05)

# Again: NOTHING SIGNIFICANT after p.adjust() ...


# Pairwise comparisons of ASV rel.ab. - HONG KONG symb vs bleached ----------------
signif_ASV_HongKong <- comp %>% 
  filter(origin_ms == "Hong Kong") %>% 
  nest(data = -ASV) %>% 
  mutate(test = map(.x = data, 
                    ~wilcox.test(rel_abund~state, data = .x, exact = F) )) %>% 
  mutate(test = map(.x = test,
                    ~broom::tidy(.x))) %>% 
  unnest(test) %>% 
  # filter(p.value < 0.05) #  15 ASVs ...
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p_adj < 0.05)

# Again again: NOTHING SIGNIFICANT after p.adjust() ...

# CONCLUSION: there is not a single ASV that has a significantly different 
# relative abundance between symbiotic and bleached samples (even when grouped by origin)

# Next: try the same but at FAMILY level



# Prep data - FAMILY level ------------------------------------------------------------------

# Add family info to ASV
fam <- taxonomy_seqs %>% 
  rename(ASV = ASV_id) %>% 
  select(ASV, uncX_FAM) %>% 
  inner_join(., comp, by = "ASV") %>%
  select(-levels_colony_ms, -colony_ms) %>% 
  arrange(state_colony_frgmt_ms, uncX_FAM) #%>% select(-ASV)

# Recalculate abs and rel abundance (by sample) now at family level
# (this pools together != ASVs that belong to the same Family)
fam2 <- fam %>% 
  group_by(state_colony_frgmt_ms, uncX_FAM) %>% 
  summarize(
            ab_fam = sum(abundance),
            rel_ab_fam = sum(rel_abund),
            .groups = "drop") 

# Sanity check: all rel_ab by family sum up to 1 for each sample
fam2 %>% group_by(state_colony_frgmt_ms) %>% 
  summarise(check_rel_ab_fam = sum(rel_ab_fam)) # YES!

# Put back the sample metadata
comp_fam <- fam %>% 
  select(state, origin_ms, state_colony_frgmt_ms) %>% 
  unique() %>% 
  inner_join(., fam2, by = "state_colony_frgmt_ms")

rm(fam, fam2)


# How many Families in total?
comp_fam$uncX_FAM %>% unique() %>% length() # 137 (from 506 ASVs)



# Pairwise comparisons of FAMILIES rel.ab. - OVERALL symb vs bleached ----------------

# For each ASV, run pairwise comparison (Wilcoxon text) of 
# its rel.ab btw symbiotic and bleached samples
signif_FAMs <- comp_fam %>% 
  nest(data = -uncX_FAM) %>% 
  mutate(test = map(.x = data, 
                    ~wilcox.test(rel_ab_fam~state, data = .x, exact = F) )) %>% 
  mutate(test = map(.x = test,
                    ~broom::tidy(.x))) %>% 
  unnest(test) %>% 
  # filter(p.value < 0.05) # 8 families ...
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p_adj < 0.05)

# Results: nothing significant ... (8 Families, but WITHOUT p adjustment)



# Pairwise comparisons of FAMILIES rel.ab. - RED SEA symb vs bleached ----------------

signif_FAMs_RedSea <- comp_fam %>% 
  filter(origin_ms == "Red Sea") %>% 
  nest(data = -uncX_FAM) %>% 
  mutate(test = map(.x = data, 
                    ~wilcox.test(rel_ab_fam~state, data = .x, exact = F) )) %>% 
  mutate(test = map(.x = test,
                    ~broom::tidy(.x))) %>% 
  unnest(test) %>% 
  # filter(p.value < 0.05) #  2 families ...
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p_adj < 0.05)

# Results: nothing significant ... (only 2 fams WITHOUT p adjustment)

  
  
# Pairwise comparisons of FAMILIES rel.ab. - HONG KONG symb vs bleached ----------------

signif_FAMs_HongKong <- comp_fam %>% 
 filter(origin_ms == "Hong Kong") %>% 
 nest(data = -uncX_FAM) %>% 
 mutate(test = map(.x = data, 
                   ~wilcox.test(rel_ab_fam~state, data = .x, exact = F) )) %>% 
 mutate(test = map(.x = test,
                   ~broom::tidy(.x))) %>% 
 unnest(test) %>% 
 # filter(p.value < 0.05) # only 8 families ...
 mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p_adj < 0.05)

# Results: nothing significant ... (only 6 fams WITHOUT p adjustment)
  
# OVERALL CONCLUSIONS: there are no ASVs nor Families that significantly differ 
# btw symbiotic and bleached samples, neither overall no by origin 
  
  
