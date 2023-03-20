
# 0. Install packages ----------------------------------------
library(here)
library(tidyverse)
library(vegan)
library(patchwork)

# install.packages("BiocManager")
# BiocManager::install("phyloseq")
library(phyloseq)


# 1. Build a phyloseq project -----------------------------------------

# Need:
#   - metadata
#   - feature table
#   - taxonomy
#   - rooted-tree

# Import data (unzipped and converted qiime2 artifacts)
# and re-shaped here to match the format required by {phyloseq}
# Note: for all 3 tables is necessary to set the first column as row names! 
# Note2: {phyloseq} always refers to OTU but we actually use ASVs (but practically for these analysis it doesn't make any difference)

# Import METADATA
metadata <- read_csv("./in/metadata_ms.csv") %>% 
   column_to_rownames(var = "sample_id")


# ... but also keep metadata in tibble format for later use (nMDS ...)
metadata_ms <- read_csv("./in/metadata_ms.csv")

# Import FEATURE TABLE
asv <- read_delim(file = "./in/q2_outputs/for_phyloseq/feature-table.tsv", 
                  delim = "\t", col_names = T,
                  skip = 1, comment = "") # turn off the interpretation of comments altogether.
asv <- asv %>% 
  column_to_rownames(var = "#OTU ID") 

# Import TAXONOMY - taxonomy cleaned in 00_data_structuring.R
# (to avoid repeating the code here to do exactly the same thing ...)
taxonomy <- read_csv("./out/Gfas_16S/useful_tables/taxonomy_all_nonraref.csv")

# Reshape TAXONOMY: duplicate 'Feature ID' so that I can keep it as a var for later
# (just in case)
taxonomy <- taxonomy %>% 
  mutate(Feature_ID = ASV_id) %>% 
  relocate(Feature_ID, .after = ASV_id) # %>% view()

# Reshape TAXONOMY: Make first column as row name 
tax_clean <- taxonomy %>% 
  column_to_rownames(var = "ASV_id")



# Make phyloseq-compatible objects
ASV <- phyloseq::otu_table(as.matrix(asv), taxa_are_rows = TRUE) # technically it's not OTUs but ASV for my data
TAX <- phyloseq::tax_table(as.matrix(tax_clean))
META <- phyloseq::sample_data(metadata)
TREE <- phyloseq::read_tree("./in/q2_outputs/for_phyloseq/tree.nwk")

# Merge the data
ps <- phyloseq::phyloseq(ASV, TAX, META, TREE)
ps


# Clean environment ...
rm(list = setdiff(ls(), c("ps", "metadata", "metadata_ms")))



# 2.1 Rarefy (one time) -------------------------------------------------

# Rarefy to sequencing depth of 2690
# which is the number of features in the second last sample (by nr of reads per sample)
# (see Jupyter notebook for details, and see plot here below)

# Plot rarefaction curve
vegan::rarecurve(t(otu_table(ps)), step = 50, cex = 0.5) # not sure what 'cex = 0.5' mans ...
                     # Can see the last and second last samples with lower nr of species (features)



# Actually rarefy
# Note that {phyloseq} implements rarefaction but highly discourages it!
ps_rarefied = phyloseq::rarefy_even_depth(physeq = ps, 
                                          rngseed = 1, # note: same as previous set.seed()
                                          sample.size = 2690, 
                                          replace = F) 
        # Set:
        #   - rngseed; = random seed, make results reproducible
        #   - replace = F; = sample without replacement
        # Note: sampling with replacement means that there is a chance that the number 
        #       of reads for a given OTU in a given sample could be larger than the 
        #       original count value, as opposed to sampling without replacement 
        #       where the original count value is the maximum possible.

# WARNING MESSAGE:
# `set.seed(1)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(1); .Random.seed` for the full vector
# ...
# 1 samples removedbecause they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 
#   
#   S-60-3-16S
# ...
# 9OTUs were removed because they are no longer 
# present in any sample after random subsampling



# Compare before-after:
# - dropped one sample (as expected)
# - lost 515 - 506 = 9 taxa

ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 515 taxa and 28 samples ]
# sample_data() Sample Data:       [ 28 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 515 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 515 tips and 514 internal nodes ]

ps_rarefied
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 506 taxa and 27 samples ]
# sample_data() Sample Data:       [ 27 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 506 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 506 tips and 505 internal nodes ]




## Subset from ps_rarefied (phyloseq proj) by origin ------------------------

ps_rarefied_RedSea <- phyloseq::subset_samples(ps_rarefied, 
                                               origin_ms == "Red Sea")

ps_rarefied_HongKong <- phyloseq::subset_samples(ps_rarefied, 
                                                 origin_ms == "Hong Kong")

ps_rarefied_RedSea
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 506 taxa and 17 samples ]
# sample_data() Sample Data:       [ 17 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 506 taxa by 10 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 506 tips and 505 internal nodes ]

ps_rarefied_HongKong
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 506 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 506 taxa by 10 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 506 tips and 505 internal nodes ]



# Export RDS ------------------
saveRDS(ps, "./out/RDS_files/ps.rds")
saveRDS(ps_rarefied, "./out/RDS_files/ps_rarefied.rds")
saveRDS(ps_rarefied_RedSea, "./out/RDS_files/ps_rarefied_RedSea.rds")
saveRDS(ps_rarefied_HongKong, "./out/RDS_files/ps_rarefied_HongKong.rds")



# 2.2 Rarefy with iteration ---------------------

# devtools::install_github("vmikk/metagMisc")
library(metagMisc)

## Examples/tryout functions -----------------

seeds <- c(1:10)

# This produces a list of x (nr of iter) ps
ps_rarefied.10 <- metagMisc::phyloseq_mult_raref(ps, 
                                                 SampSize = 2690,
                                                 iter = 10,
                                                 replace = F,
                                                 seeds = seeds)
ps_rarefied.10$`1`
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 506 taxa and 27 samples ]
# sample_data() Sample Data:       [ 27 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 506 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 506 tips and 505 internal nodes ]

ps_rarefied.10[2] # diff way of calling elem of list
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 508 taxa and 27 samples ]
# sample_data() Sample Data:       [ 27 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 508 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 508 tips and 507 internal nodes ]

# etc ...



# Now let's try the average method, which
# calculates the average and outputs one single ps
ps_rarefied.10avg <- metagMisc::phyloseq_mult_raref_avg(ps, 
                                                 SampSize = 2690,
                                                 iter = 10,
                                                 replace = F,
                                                 seeds = seeds)

ps_rarefied.10avg
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 509 taxa and 27 samples ]
# sample_data() Sample Data:       [ 27 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 509 taxa by 10 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 509 tips and 508 internal nodes ]


rm(list = ls(pattern = "ps_rarefied.10")) # Removes all objects whose name matches the pattern



# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ----------------------

## Rarefied ps (1000 iteration) to keep for downstream analysis -------------------

###  !!! WARNING !!! 
###      COMPUTATIONALLY INTESE: DON'T RE-RUN UNLESS NECESSARY !!!
###         -> use .RDS file instead ;) 

# Set seeds for rep results
# Note: vector needs to be same size as nr of iterations
iterations <- 1000
seeds <- c(1:iterations) # Note that for alphas I use exactly the same seeds
                         # to for loop `phyloseq::rarefy_even_depth()``  

ps_rarefied.1000avg <- metagMisc::phyloseq_mult_raref_avg(ps, 
                                                          SampSize = 2690,
                                                          iter = iterations,
                                                          replace = F,
                                                          seeds = seeds)
ps_rarefied.1000avg

# # Export RDS ------------------
# saveRDS(ps_rarefied.1000avg, "./out/RDS_files/ps_rarefied.1000avg.rds")




ps_rarefied.1000avg_RedSea <- phyloseq::subset_samples(ps_rarefied.1000avg, 
                                               origin_ms == "Red Sea")
ps_rarefied.1000avg_RedSea
# otu_table()   OTU Table:         [ 509 taxa and 17 samples ]
# sample_data() Sample Data:       [ 17 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 509 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 509 tips and 508 internal nodes ]

ps_rarefied.1000avg_HongKong <- phyloseq::subset_samples(ps_rarefied.1000avg, 
                                                 origin_ms == "Hong Kong")
ps_rarefied.1000avg_HongKong
# otu_table()   OTU Table:         [ 509 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 509 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 509 tips and 508 internal nodes ]



# Export RDS ------------------
# saveRDS(ps, "./out/RDS_files/ps.rds")
# saveRDS(ps_rarefied, "./out/RDS_files/ps_rarefied.rds")
# saveRDS(ps_rarefied_RedSea, "./out/RDS_files/ps_rarefied_RedSea.rds")
# saveRDS(ps_rarefied_HongKong, "./out/RDS_files/ps_rarefied_HongKong.rds")

saveRDS(ps_rarefied.1000avg, "./out/RDS_files/ps_rarefied.1000avg.rds")
saveRDS(ps_rarefied.1000avg_RedSea, "./out/RDS_files/ps_rarefied.1000avg_RedSea.rds")
saveRDS(ps_rarefied.1000avg_HongKong, "./out/RDS_files/ps_rarefied.1000avg_HongKong.rds")

# xxxxxxxxxxxxxxxxxxxxxxxxxxxx --------------

# Import RDS (to avoid re-running iterated rar.) ------------------------
ps_rarefied.1000avg <- readRDS( "./out/RDS_files/ps_rarefied.1000avg.rds")
ps_rarefied.1000avg_RedSea <- readRDS("./out/RDS_files/ps_rarefied.1000avg_RedSea.rds")
ps_rarefied.1000avg_HongKong <- readRDS("./out/RDS_files/ps_rarefied.1000avg_HongKong.rds")




# Plot rarefaction curve in {ggplot2} ----------------------------
# because rarecurve plots in base r is ugly ...

# For the supplementary, to see if we reach saturation 
# (curves reach asymptotic-like shape or still steep? YES)
# and compare before and after rarefaction


## Rarefaction curve - NON-rarefied data -------------------

# Use ps

rarecurve_data <- vegan::rarecurve(t(otu_table(ps)), step = 50)

rarecurve_data_tibble <- map_dfr(rarecurve_data, bind_rows) %>% 
  bind_cols(sample_id = rownames(t(otu_table(ps)), .)) %>% 
  pivot_longer(-sample_id) %>% 
  drop_na() %>% #view()
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>% 
  rename(n_ASVs = value) %>% #view()
  left_join(., metadata_ms, by = "sample_id") %>% 
  group_by(sample_id) %>% 
  mutate(label = ifelse(n_seqs == max(n_seqs), state_colony_frgmt_ms, NA)) %>% 
  ungroup() #%>%   view()

# Palette
source("colors_GfasMS.R")
palette_GfasMS

ggplot(data = rarecurve_data_tibble, 
       aes(x = n_seqs, y = n_ASVs, group = sample_id)) +
  geom_vline(xintercept = 2690, linetype = 2) +
  geom_line(aes(color = colony_ms)) +
  geom_text(aes(label = label), 
            hjust = 0,
            size = 3) + 
  scale_color_manual(values = palette_GfasMS) +
  labs(
    x = "number of sequences",
    y = "number of ASVs",
    color = "Colony:"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(
    override.aes=list(
      # shape = 18,
      linetype = 1,
      size = 3)))


ggsave("./out/Gfas_16S/supplementary/rarecurve_nonrarefied_color.png",
       bg = "white",
       dpi = 330, 
       units = "cm", width = 15, height = 15)




# Clean up ----------------
# rm(list = setdiff(ls(), c("ps", "ps_rarefied", "metadata_ms")))
rm(rarecurve_data)
rm(rarecurve_data_tibble)
