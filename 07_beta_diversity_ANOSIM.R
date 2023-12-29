######################################################
#
#   ANOSIM - Analysis of similarities 
#   Here performed at ASV level on RAREFIED data
# 
#   Why:
#   Test if Red Sea and Hong Kong are significantly different 
#   (they are with PERMANOVA)
# 
#   This analysis was added during peer-review process as suggested by reviewer
#
######################################################

library(tidyverse)
library(here)
library(phyloseq)


# Import data -----------------------


# rarefied ps (phyloseq) 1000 iterations (same used for PERMANOVA and PERMDISP)

# ps_rarefied <- readRDS("./out/RDS_files/ps_rarefied.rds") #OLD version 1x iteration
ps_rarefied.1000avg <- readRDS( "./out/RDS_files/ps_rarefied.1000avg.rds")

# metadata
metadata_ms <- read_csv("./in/metadata_ms.csv")
# metadata_physeq <- data.frame(sample_data(ps_rarefied.1000avg))



# Data preparation ----------------------------------------------------------------

# OTU table is what I need, but ...
# make it matrix (instead of phyloseq class object) then
# make it tibble for easier wrangling
ASV_table <- as(otu_table(ps_rarefied.1000avg), "matrix") %>% # OLD: ps_rarefied
  t() %>% # transpose
  as_tibble(., rownames = "sample_id")

ASV_table_meta <- inner_join(metadata_ms, ASV_table, by = "sample_id") # OLD: metadata_ms

# Keep only cols with ASV names (Qiime2 generated alphanumeric names)
#  = only ASVs, not metadata cols

# ASV_table_meta[, 11:516] # check that the right cols are selected = only ASVs, not metadata cols
ASV_table_meta %>% names() %>% nchar() # ASV names have 32 characters
for_anosim <- ASV_table_meta %>%
  select(names(.)[nchar(names(.)) == 32])

# sanity check: retained the right number of cols? 
#  (knowing from visual inspection that first 10 cols are metadata)
length(names(ASV_table_meta)) - 10 == length(for_anosim) # TRUE :)



# ANOSIM =======================================================================

# ANOSIM tests whether distances between groups are greater than within groups. 
# (PERMANOVA tests whether distance differ between groups.)

# non-parametric test, but sensitive to unbalanced designs 
# and my data is not balanced:
#  Red Sea:   n = 17
#  Hong Kong: n = 10


# Analysis of similarities (ANOSIM) provides a way to test statistically 
# whether there is a significant difference between two or more groups of sampling units. 
# 
# Function `anosim()` operates directly on a dissimilarity matrix. 
# A suitable dissimilarity matrix is produced by functions dist or vegdist. 
# The method is philosophically allied with NMDS ordination (monoMDS), in that 
# it uses only the rank order of dissimilarity values.
# If two groups of sampling units are really different in their 
# species composition, then compositional dissimilarities between the groups 
# ought to be greater than those within the groups. 


# 'Between' for dissimilarities between classes and 
# 'class_name' for corresponding dissimilarity within class


# ANOSIM - by origin ------------------------------------------------------------

# Use 'origin' variable for grouping
ANOSIM_origin <- vegan::anosim(x = for_anosim, grouping = ASV_table_meta$origin, 
                              distance = "bray", permutations = 9999)
ANOSIM_origin
# ANOSIM statistic R: 0.2015 
# Significance: 0.0186 *


summary(ANOSIM_origin)
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0911 0.1366 0.1864 0.2453 
# 
# Dissimilarity ranks between and within classes:
#         0%    25%   50%    75% 100%   N
# Between  9 111.50 186.5 286.75  349 170
# HKU      4  24.00  47.0  77.00  116  45
# JLU      1 143.75 205.5 252.25  349 136

plot(ANOSIM_origin)

# Export plot
png("out/Gfas_16S/beta_diversity/rar1000_Bray_ANOSIM.png", 
    units = "cm", width = 10, height = 10,
    res = 300)
plot(ANOSIM_origin)
dev.off()

# How to interpret plot:
# - Only compare "Between" (mean distance among all sites/samples) vs.
#    individual groups (here "HKU" and "JLU")
# - NOT meaningful to compare groups with each other (e.g. cannot directly compare "HKU" vs "JLU")
# - notch = 95 % CI around the median. If no overlap = significant difference


