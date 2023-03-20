#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Look at "core" microbiome
#   At ASV level on NON-RAREFIED data
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(here)
library(phyloseq)

# Import data and source color palette ----------------------------------------

# Import {phyloseq} project
ps <- readRDS("./out/RDS_files/ps.rds")

# Color palette
source("colors_GfasMS.R")


# Import taxonomy table
# taxonomy_seqs <- read_csv("./out/Gfas_16S/taxonomy_seqs.csv") # need to use taxonomy from non-rarefied!!
taxonomy_seqs <- read_csv("./out/Gfas_16S/useful_tables/taxonomy_all_nonraref.csv") %>% 
  rename(ASV = ASV_id)


# Import metadata
metadata_ms <- read_csv("./in/metadata_ms.csv")


# Data wrangling ----------------------

# Get ASV table out of "ps" (nonrarefied, phyloseq object)
ASV_table <- as(otu_table(ps), "matrix") %>% 
  t() %>%
  as_tibble(., rownames = "sample_id")

# Rename samples to match ms format/naming
data <- metadata_ms %>% 
  select(sample_id, state, origin_ms, colony_ms, state_colony_frgmt_ms, levels_colony_ms) %>% 
  inner_join(., ASV_table, by = "sample_id") %>% 
  select(-sample_id) %>% 
  mutate(levels_colony_ms = factor(levels_colony_ms))

# Make long 
data <- data %>% pivot_longer(., cols = where(is.numeric), 
                               names_to = "ASV", 
                               values_to = "abundance") 

# Add abundance transformations
data <- data %>% 
  mutate(
    abundance_logx1 = log(abundance + 1),
    pres_abs = if_else(abundance == 0, 0, 1)
  )



# Summary table(s) ----------------------------------

# Create table to summarize the ASV occurrence across the whole data set and
## 1. Order ASV by overall presence/absence: var 'order_bynrofsamples' ------------
ASV_occurrence_overall <- data %>%   
  group_by(ASV) %>% 
  summarise(
    nr_of_samples = sum(pres_abs),
    perc_of_samples = round(sum(pres_abs)/nrow(ASV_table) * 100, 1),
    overall_abund = sum(abundance),
    perc_overallabund = round(sum(abundance)/sum(ASV_table[,-1]) * 100, 1), # sum(data$abundance) = 72630
  ) %>% 
  ungroup() %>% 
  arrange(desc(nr_of_samples)) %>% 
  mutate(order_bynrofsamples = row_number()) %>% 
  relocate(order_bynrofsamples, .after = ASV) %>% 
  inner_join(., taxonomy_seqs, by = "ASV") #%>% view()

## 2. order ASV by overall abundance (relative to the tot nr of reads in the data set) ----------
ASV_occurrence_overall <- ASV_occurrence_overall %>% 
  arrange(desc(overall_abund)) %>%  # same as using 'perc_overallabund'
  mutate(
    order_byabund = row_number()
  ) %>%
  relocate(order_byabund, .after = perc_of_samples)

## Create easier names for ASV based on their abundance (ASV_001 = most abundant) --------
# Make table to reference ASV original id (32 chr from q2) to a numbered name
# as 'ASV_n', where n = the ranking of that ASV by overall abundance (across whole data set)
ASV_numbered <- ASV_occurrence_overall %>% 
  mutate(
    ASV_nr = paste("ASV", str_pad(
      order_byabund, width = 3, side = "left", pad = "0"), 
      sep = "_")
  ) %>% 
  relocate(ASV_nr, .after = ASV)

# Subset for "ASV_numbers" 
ASV_numbers <- ASV_numbered %>% select(ASV, ASV_nr)

write_csv(ASV_numbers, "./out/Gfas_16S/core_mb/nonrarefied/ASV_numbers.csv")


## Calculate % and abundance by symbiotic state - SYMBIOTIC -----------------------
ASV_occurrence_symb <- data %>%   
  filter(state == "symbiotic") 

ASV_occurrence_symb <- ASV_occurrence_symb %>% 
  group_by(ASV) %>% 
  summarise(
    nr_of_samples_symb = sum(pres_abs),
    perc_of_samples_symb = round(sum(pres_abs)/length(unique(ASV_occurrence_symb$state_colony_frgmt_ms)) * 100, 1),
    overall_abund_symb = sum(abundance),
    perc_overallabund_symb = round(sum(abundance)/sum(ASV_occurrence_symb$abundance) * 100, 1), # sum(data$abundance) = 72630
  ) %>% 
  ungroup() 

## Calculate % and abundance by symbiotic state - BLEACHED -----------------------
ASV_occurrence_blea <- data %>%   
  filter(state == "bleached") 

ASV_occurrence_blea <- ASV_occurrence_blea %>% 
  group_by(ASV) %>% 
  summarise(
    nr_of_samples_blea = sum(pres_abs),
    perc_of_samples_blea = round(sum(pres_abs)/length(unique(ASV_occurrence_blea$state_colony_frgmt_ms)) * 100, 1),
    overall_abund_blea = sum(abundance),
    perc_overallabund_blea = round(sum(abundance)/sum(ASV_occurrence_blea$abundance) * 100, 1), # sum(data$abundance) = 72630
  ) %>% 
  ungroup() 



## Put OVERALL, SYMBIOTIC, and BLEACHED together -----------------------
ASV_occurrence_combined <- left_join(ASV_occurrence_overall, ASV_occurrence_symb,
                                     by = "ASV") %>% 
  left_join(., ASV_occurrence_blea, by = "ASV") %>%
  select(
    "Kingdom", "Phylum", "Class", "Order", "Family", "uncX_FAM", "Genus", "Species",
    "Sequence", "ASV", # ASV_nr 
    "nr_of_samples", "nr_of_samples_symb", "nr_of_samples_blea",
    "perc_of_samples", "perc_of_samples_symb", "perc_of_samples_blea",
    "overall_abund", "overall_abund_symb", "overall_abund_blea",
    "perc_overallabund", "perc_overallabund_symb", "perc_overallabund_blea"
  ) %>% 
  inner_join(., ASV_numbers, by = "ASV") %>% 
  relocate(ASV_nr, .after = ASV) %>% 
  arrange(desc(perc_of_samples))


### Export tables -------------------

# Export ASV_occurrence_combined
ASV_occurrence_combined %>% 
  write_csv(., "./out/Gfas_16S/core_mb/nonrarefied/ASV_occurrence_summary_all.csv")

# Export only the ones that are overall present in >= 50 % of samples (at least in 13 samples)
ASV_occurrence_combined %>% 
  filter(perc_of_samples >= 50) %>% 
  write_csv("./out/Gfas_16S/core_mb/nonrarefied/ASV_occurrence_summary_50percsamples.csv")






# Look at cumulative sum of abundance (%) -------------------------------------
# aka how many (and which) ASVs contribute to X % of the tot nr of reads?
ASV_numbered_cumsum <- ASV_numbered %>% 
  mutate(perc_ovallabund_cumsum = cumsum(perc_overallabund) ) %>% 
  relocate(perc_ovallabund_cumsum, .after = perc_overallabund) 

write_csv(ASV_numbered_cumsum, "./out/Gfas_16S/core_mb/nonrarefied/ASV_occurence_overall.csv")


# Plot OVERALL ASV accumulation (cum. ab.) curve -------------------------------
ASV_numbered_cumsum %>% 
  filter(perc_overallabund > 0.1) %>% 
  ggplot(., 
       aes(x = order_byabund,
         # x =fct_reorder(ASV_nr, order_byabund), 
           y = perc_ovallabund_cumsum)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 100), 
                     labels = scales::percent_format(scale = 1),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 150, 10), 
                     expand = c(0, 0)) +
  xlab("nr of ASVs") +
  ylab("cumulative abundance\n(across the whole data set)") +
  labs(caption = "Showing only ASVs with overall abundance > 0.01 %") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.major.x = element_line(colour = "grey75", linetype = "dotted")
  )

ggsave("./out/Gfas_16S/core_mb/nonrarefied/ASV_cum_abund.png",
       units = "cm",
       width = 17, height = 10)



# PLOT OVERALL + GROUPS (state x origin) ASV accumulation curve (cum. ab.) -------------------------------

## Create the 4 tibbles (one per group) -----------
data2 <- data %>% select(-c(pres_abs, abundance_logx1)) # dez not needed

ASV_RS_symb <- data2 %>% 
  filter(origin_ms == "Red Sea" & state == "symbiotic")

ASV_RS_blea <- data2 %>% 
  filter(origin_ms == "Red Sea" & state == "bleached")

ASV_HK_symb <- data2 %>% 
  filter(origin_ms == "Hong Kong" & state == "symbiotic")

ASV_HK_blea <- data2 %>% 
  filter(origin_ms == "Hong Kong" & state == "bleached")

rm(data2)

## Function to apply same transformation ---------------
# calculate rel. ab. (as %) > order by desc rel. ab. > cumsum()

ASVcumsum <- function(data) {
  data_cumsum <- data %>% 
    filter(abundance > 0) %>% 
    group_by(ASV) %>% 
    summarise(abundance_tot = sum(abundance)) %>% 
    ungroup() %>% 
    mutate(
      rel_ab = abundance_tot / sum(abundance_tot) * 100) %>% 
    arrange(desc(abundance_tot)) %>% 
    mutate(
      rank = row_number(),
      rel_ab_cumsum = cumsum(rel_ab)
    )
  
  data %>% 
    select(ASV, origin_ms, state) %>% 
    mutate(
      origin_state = paste(origin_ms, state, sep = " ")
    ) %>% 
    left_join(data_cumsum, ., by = "ASV")
}


##  Apply function --------------
ASV_RS_symb_cumsum <- ASVcumsum(ASV_RS_symb) 
ASV_RS_blea_cumsum <- ASVcumsum(ASV_RS_blea) 
ASV_HK_symb_cumsum <- ASVcumsum(ASV_HK_symb) 
ASV_HK_blea_cumsum <- ASVcumsum(ASV_HK_blea) 



### Make long by binding all tibbles together -----------

#### Adjust overall tibble to match vars and naming with the others -----------
ASV_all_cumsum <- ASV_numbered_cumsum %>% 
  rename(
    rank = order_byabund,
    rel_ab_cumsum = perc_ovallabund_cumsum
  ) %>% 
  mutate(
    origin_state = "Overall"
  ) #%>% view()


#### Select same vars from all tibbles -------------
keepvars <- function(x) {
  select(x, ASV, rank, rel_ab_cumsum, origin_state)
}

ASV_all_cumsum <- keepvars(ASV_all_cumsum)
ASV_RS_symb_cumsum <- keepvars(ASV_RS_symb_cumsum)
ASV_RS_blea_cumsum <- keepvars(ASV_RS_blea_cumsum)
ASV_HK_symb_cumsum <- keepvars(ASV_HK_symb_cumsum)
ASV_HK_blea_cumsum <- keepvars(ASV_HK_blea_cumsum)

names(ASV_all_cumsum) == names(ASV_RS_symb_cumsum) # TRUE TRUE TRUE TRUE


#### rbind() into one long tibble --------------------------
ASV_cumsum_bind <- rbind(ASV_all_cumsum, 
                         ASV_RS_symb_cumsum, ASV_RS_blea_cumsum,
                         ASV_HK_symb_cumsum, ASV_HK_blea_cumsum)


#### Plot ASV accumulation curve for all 5 groups -------------

(
ggplot_acccurve_5gr <- ggplot(data = ASV_cumsum_bind,
                                   
            aes(x = rank, y = rel_ab_cumsum,
                color = origin_state,
                fill = origin_state)) +
       geom_point(shape = 21,
                  size = 2) + 
       scale_color_manual(breaks = c("Red Sea symbiotic", "Red Sea bleached", 
                                     "Hong Kong symbiotic", "Hong Kong bleached",
                                     "Overall"),
                          values = c(Red_Sea, Red_Sea,
                                     Hong_Kong, Hong_Kong,
                                     "black")) +
       scale_fill_manual(breaks = c("Red Sea symbiotic", "Red Sea bleached", 
                                    "Hong Kong symbiotic", "Hong Kong bleached",
                                    "Overall"),
                          values = c(Red_Sea, "transparent", #"white",
                                     Hong_Kong, "transparent", # "white",
                                     "black")) +
       scale_y_continuous(limits = c(0, 100), 
                          labels = scales::percent_format(scale = 1),
                          expand = c(0, 0)) +
       scale_x_continuous(limits = c(0, 60), 
                          breaks = seq(0, 150, 10), 
                          expand = c(0, 0)) +
       xlab("nr of ASVs") +
       ylab("cumulative abundance") +
       theme_classic() +
       theme(
         axis.ticks = element_blank(),
         panel.grid.major.y = element_line(colour = "grey", linetype = "dotted"),
         panel.grid.major.x = element_line(colour = "grey75", linetype = "dotted"),
         # legend.position = "none"
         legend.title = element_blank(),
         legend.position = "bottom"
       )
)


ggsave("./out/Gfas_16S/core_mb/nonrarefied/ASV_cum_abund_5groups_new.png",
       dpi = 330,
       units = "cm",
       width = 16, height = 9)


# Export RDS for composite plot (with alpha diversity) ------------------
saveRDS(ggplot_acccurve_5gr, "./out/RDS_files/ggplot_acccurve_5gr.rds")
saveRDS(ASV_cumsum_bind, "./out/RDS_files/ASV_cumsum_bind.rds")


# Plot: ASVs Occurrence ~ Abundance (supplementary) -------------------------
# How well abundance and presence scale?

ggplot(data = ASV_occurrence_overall,
       aes(x = overall_abund,
           y = nr_of_samples)) +
  geom_point() +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(title = "Relationship between ASV occurrence and abundance\n(non-rarefied data)") +
  xlab("ASV abundance (across the whole data set)") +
  ylab("nr of samples\nin which the ASV is found") +
  theme_classic() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 11)
  )

ggsave("./out/Gfas_16S/supplementary/nonrarefied_presence_abundance.png",
         dpi = 300,
         units = "cm", height = 10, width = 15 )






# Extract list of 28 "core" ASVs ----------------------------

# The ones found across all groups (but not necessarily in all samples ... !)

uniqueASVs <- function (x) {
  x %>% 
    filter(abundance != 0) %>% 
    pull(ASV) %>% 
    unique()
}

ASV_RS_symb_unique <- uniqueASVs(ASV_RS_symb)
ASV_RS_blea_unique <- uniqueASVs(ASV_RS_blea)
ASV_HK_symb_unique <- uniqueASVs(ASV_HK_symb)
ASV_HK_blea_unique <- uniqueASVs(ASV_HK_blea)


# Get 'ASV' (q2 names)
coreASV <- intersect(ASV_RS_symb_unique, ASV_RS_blea_unique) %>% 
  intersect(., ASV_HK_symb_unique) %>% 
  intersect(., ASV_HK_blea_unique)

# Get ASV_nr 
coreASV_occurrence_combined <- ASV_occurrence_combined %>% 
  filter(ASV %in% coreASV) 
  # select(ASV, ASV_nr, Sequence:Species) %>% view()

# Export table
write_csv(coreASV_occurrence_combined, "./out/Gfas_16S/core_mb/nonrarefied/core28ASVs.csv")
  



# PLOT: tiles for overview of ALL ASVs across all samples ---------

## Add order of ASV (by presence) to 'data' for ordering ASV in the plot ------
plot_all <- ASV_occurrence_overall %>%
  select(ASV, order_bynrofsamples) %>%
  left_join(., data, by = "ASV") %>%
  arrange(state_colony_frgmt_ms, order_bynrofsamples) #%>% view()


## Plot ALL samples, ORDERED PRESENCE/ABSENCE ------------------------------------------------
ggplot(data = plot_all, aes(x = forcats::fct_reorder(state_colony_frgmt_ms, as.numeric(levels_colony_ms)),
                        y = forcats::fct_reorder(ASV, desc(order_bynrofsamples)),
                        alpha = factor(pres_abs),
                        fill = colony_ms)) +
  geom_tile() +
  # coord_fixed(ratio = 1) + # square tiles (1:1 x:y)
  # scale_fill_gradient(low = "black", high = "white") +
  scale_alpha_manual(values = c(0, 1), breaks = c("0", "1")) +
  scale_fill_manual(values = palette_GfasMS) +
  # scale_fill_manual(values = c("white", "black"), breaks = c("0", "1")) +
  scale_x_discrete(position = "top") +
  labs(title = "ASVs presence/absence (rarefied data)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "NA", color = "grey69"),
    strip.placement = "outside",
    legend.position = "none",
    panel.background = element_blank()
  ) + 
  facet_wrap(~ forcats::fct_rev(state), 
             scales = "free_x")


ggsave("./out/Gfas_16S/core_mb/nonrarefied/ASV_occurence_overall_pres_abs.png",
         dpi = 300,
         units = "cm", height = 45, width = 10 )




# Plot ALL samples, ORDERED ABUNDANCE ------------------------------------------------
ggplot(data = plot_all, aes(x = forcats::fct_reorder(state_colony_frgmt_ms, as.numeric(levels_colony_ms)),
                        y = forcats::fct_reorder(ASV, desc(order_bynrofsamples)),
                        alpha = abundance_logx1, # abundance
                        fill = colony_ms)) +
  geom_tile() +
  # coord_fixed(ratio = 1) + # square tiles (1:1 x:y)
  # scale_fill_gradient(low = "black", high = "white") +
  # scale_alpha_manual(values = c(0, 1), breaks = c("0", "1")) +
  scale_alpha(range = c(0.00001, 1)) +
  scale_fill_manual(values = palette_GfasMS) +
  # scale_fill_manual(values = c("white", "black"), breaks = c("0", "1")) +
  scale_x_discrete(position = "top") +
  labs(title = "ASVs log(x+1) transf. abundance (rarefied data)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "NA", color = "grey69"),
    strip.placement = "outside",
    # legend.position = "none",
    panel.background = element_blank()
  ) + 
  facet_wrap(~ forcats::fct_rev(state), 
             scales = "free_x")


ggsave("./out/Gfas_16S/core_mb/nonrarefied/ASV_occurence_overall_logx1.png",
       dpi = 300,
       units = "cm", height = 30, width = 13 )



# Summary tables of ASVs occurence by group ----------------------

# Look at ASV occurrence in sub group: state and/or origin

# Produce 6 summary tables (and export to ./*/core_mb/nonrarefied/ ) 
# 1.v Summary of ASV occurrence across samples - overall [done above line 100-115]
# 2. Summary of ASV occurrence across samples - symbiotic (all)
# 3. Summary of ASV occurrence across samples - bleached (all)
# 4. Summary of ASV occurrence across samples - symbiotic Red Sea
# 5. Summary of ASV occurrence across samples - bleached Red Sea
# 6. Summary of ASV occurrence across samples - symbiotic Hong Kong
# 7. Summary of ASV occurrence across samples - bleached Hong Kong



### Write function to make summary tables - 1 factor ---------------------------------------
summarize_ASV_occurrence1 <- function(data, filterby, writeout) {
  
  a <- data %>%
    filter(state == filterby) 
  
  # tot nr of samples in the data (sub)set
  tot_nr_samples <- data %>% 
    filter(state == filterby) %>% 
    pull(state_colony_frgmt_ms) %>% unique() %>% length()
  
  # tot_abundance = tot nr of reads in the data (sub)set = sum of all abund of all ASV
  tot_abund <- a %>% pull(abundance) %>% sum()
  
  summary <- a %>% 
    group_by(ASV) %>% 
    summarise(
      nr_of_samples = sum(pres_abs),
      perc_of_samples = round(sum(pres_abs)/tot_nr_samples * 100, 1),
      abundance = sum(abundance),
      perc_abund = round(sum(abundance)/tot_abund * 100, 1),
    ) %>% 
    ungroup() %>% 
    arrange(desc(nr_of_samples)) %>% 
    mutate(order_bynrofsamples = row_number()) %>% 
    arrange(desc(abundance)) %>%  # same as using 'perc_abund'
    mutate(
      order_byabund = row_number(),
      perc_abund_cumsum = cumsum(perc_abund)
    ) %>%
    inner_join(., taxonomy_seqs, by = "ASV") %>% 
    inner_join(., ASV_numbers, by = "ASV") %>% 
    filter(nr_of_samples != 0)
  
  if (writeout == "yes") {
    write_csv(x = summary, 
              file = paste0("./out/Gfas_16S/core_mb/nonrarefied/ASV_occurrence_", filterby, ".csv") )
  } else {
    view(summary)
  }
}

# Run function
summarize_ASV_occurrence1(data = data, filterby = "symbiotic", writeout = "yes")
summarize_ASV_occurrence1(data = data, filterby = "bleached", writeout = "yes")


### Write function to make summary tables - 2 factors ---------------------------------------

# NEEDS TO BE ADJUSTED (ADD VARS AS ABOVE) *****************************************
summarize_ASV_occurrence2 <- function(data, f_state, f_origin, writeout) {

      a <- data %>%
        filter(state == f_state & origin_ms == f_origin) 
      
      # tot nr of samples in the data (sub)set
      tot_nr_samples <- data %>% 
        filter(state == f_state & origin_ms == f_origin) %>% 
        pull(state_colony_frgmt_ms) %>% unique() %>% length()
      
      # tot_abundance = tot nr of reads in the data (sub)set = sum of all abund of all ASV
      tot_abund <- a %>% pull(abundance) %>% sum()
      
      summary <- a %>% 
        group_by(ASV) %>% 
        summarise(
          nr_of_samples = sum(pres_abs),
          perc_of_samples = round(sum(pres_abs)/tot_nr_samples * 100, 1),
          abundance = sum(abundance),
          perc_abund = round(sum(abundance)/tot_abund * 100, 1),
        ) %>% 
        ungroup() %>% 
        arrange(desc(nr_of_samples)) %>% 
        mutate(order_bynrofsamples = row_number()) %>% 
        relocate(order_bynrofsamples, .after = ASV) %>% 
        inner_join(., taxonomy_seqs, by = "ASV")
      
      if (writeout == "yes") {
        write_csv(x = summary, 
                  file = paste0("./out/Gfas_16S/core_mb/nonrarefied/ASV_occurrence_", 
                                f_state, "_", f_origin, ".csv") )
      } else {
        return(summary)
      }
}

# Use function
summarize_ASV_occurrence2(data = data, 
                          f_state = "symbiotic", f_origin = "Red Sea", writeout = "yes") 
summarize_ASV_occurrence2(data = data, 
                          f_state = "bleached", f_origin = "Red Sea", writeout = "yes")
summarize_ASV_occurrence2(data = data, 
                          f_state = "symbiotic", f_origin = "Hong Kong", writeout = "yes")
summarize_ASV_occurrence2(data = data, 
                          f_state = "bleached", f_origin = "Hong Kong", writeout = "yes")









# Heatmap of CORE microbiome - RELATIVE ABUNDANCE (percent) -------------------------------------------

### Select "Core" ASVs: 70 % and 60 - 70% cutoff ---------------------------
# Keep only ASVs found in >70 % of the samples (overall, whole data set)
core70 <- ASV_occurrence_combined %>% 
  filter(perc_of_samples >= 70) %>% 
  select(ASV, ASV_nr)

core6070 <- ASV_occurrence_combined %>% 
  filter(perc_of_samples < 70 & perc_of_samples >= 60) %>% 
  select(ASV, ASV_nr)



### Calculate RELATIVE ABUNDANCE ---------------------------
data <- data %>% 
  group_by(state_colony_frgmt_ms) %>%
  summarise(
    tot_reads_bysample = sum(abundance)
  ) %>%
  left_join(data, ., by = "state_colony_frgmt_ms") %>% 
  mutate(
    rel_abund_bysample = abundance / tot_reads_bysample * 100
  ) %>% 
  relocate(rel_abund_bysample, .after = abundance) 
  # group_by(state_colony_frgmt_ms) %>%                 # sanity check
  # summarise(tot_relab = sum(rel_abund_bysample)) %>% view()




### Add taxonomy information to ASV_nr ---------------------------

# Create tibble with variable to append taxonomic info
tax4plot <- ASV_numbered %>% 
  mutate(
    append_taxon = ifelse(Genus == "Unclassified" | Genus == "uncultured",
                          uncX_FAM, Genus),
    ASV_nr_taxon = paste(ASV_nr, append_taxon, sep = " ")
  ) %>% 
  select(ASV, ASV_nr, uncX_FAM, Genus, append_taxon, ASV_nr_taxon)

# Join with data
data <- left_join(data, tax4plot, by = "ASV")



### Subset 'data' (ASV table + metadata -> made long) ---------------------

data70 <- data %>% 
  # drop_na() %>% 
  semi_join(., core70, by = c("ASV", "ASV_nr")) %>% 
  left_join(., ASV_numbers, by = c("ASV", "ASV_nr"))# all rows from x with a match in y
# pull(ASV) #%>% # unique()  view()    

data6070 <- data %>% 
  # drop_na() %>%
  semi_join(., core6070, by = c("ASV", "ASV_nr")) %>% 
  left_join(., ASV_numbers, by = c("ASV", "ASV_nr"))



## Plots to patch --------------------

# # Font
# windowsFonts(my_font = windowsFont("Calibri")) 
# theme_set(theme_classic(base_family = "my_font"))




### Set upper limit for scale_fill_gradient ----------------
upper_lim <- bind_rows(data70, data6070) %>% 
  pull(rel_abund_bysample) %>% 
  max()



# Plot: BUBBLEplot of CORE microbiome - RELATIVE ABUNDANCE (percent) -------------------------------------------

# Add labels for grouping by rel. ab. cutoff
data70 <- data70 %>% 
  mutate(cutoff = "**> 70 %**<br>of all samples") 

data6070 <- data6070 %>% 
  mutate(cutoff = "**60-70 %**<br>of all samples") 

data607070 <- rbind(data70, data6070)


data607070 <- data607070 %>% mutate(State = stringr::str_to_title(state))

ggplot(data = data607070, 
       aes(y = forcats::fct_rev(ASV_nr_taxon), #forcats::fct_reorder(family, -fam_order), 
           x = forcats::fct_reorder(state_colony_frgmt_ms, as.numeric(levels_colony_ms))  
       )) +
  geom_point(aes(#size = rel_abund_bysample,
                 size=ifelse(rel_abund_bysample==0, NA, rel_abund_bysample),
                 color = colony_ms), 
             shape = 16, alpha = 0.21) +
  geom_point(aes(#size = rel_abund_bysample, 
                 size=ifelse(rel_abund_bysample==0, NA, rel_abund_bysample),
                 shape = State, 
                 color = colony_ms)) +
  scale_shape_manual(breaks = c("Symbiotic", "Bleached"), 
                     values = c(16, 21)) + 
  scale_color_manual(values = palette_GfasMS) + 
  scale_size_area(name = "Relative abundance",
                  max_size = 10, # 12
                  # n.breaks = 6,
                  breaks = c(0, 1, 10, 20, 30),
                  labels = scales::label_percent(
                    scale = 1,
                    accuracy = 1,
                    suffix = " %")
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 0, size = 9),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9),
    strip.background.x = element_rect(fill = "NA", color = "NA"), # 
    strip.background.y = element_rect(fill = "NA", color = "black"),
    strip.placement = "outside",
    strip.text.y = ggtext::element_markdown(size = 8),
    strip.text.x = element_text(size = 11, vjust = 3, face = "bold"),
    legend.position = "bottom"
  ) +
  guides(shape = "none", color = "none") +
  # facet_wrap(~forcats::fct_rev(State) + forcats::fct_reorder(colony_ms, as.numeric(levels_colony_ms)), 
  #            scales = "free_x",
  #            nrow = 1, strip.position = "top")
  facet_grid(cutoff ~ forcats::fct_rev(State), #+forcats::fct_reorder(colony_ms, as.numeric(levels_colony_ms)),
             scales = "free", space = "free")


ggsave("./out/Gfas_16S/core_mb/nonrarefied/bubble_core_relab_bystate.png",
       bg = "white",
       dpi = 330,
       units = "cm", width = 25, height = 9)

ggsave("./out/Gfas_16S/core_mb/nonrarefied/bubble_core_relab_bystate.svg",
       # bg = "white",
       # dpi = 330,
       units = "cm", width = 25, height = 9)



# Plot: NO transparency BUBBLEplot of CORE microbiome - RELATIVE ABUNDANCE (percent) -------------------------------------------

# # Add labels for grouping by rel. ab. cutoff
# data70 <- data70 %>% 
#   mutate(cutoff = "**> 70 %**<br>of all samples") 
# 
# data6070 <- data6070 %>% 
#   mutate(cutoff = "**60-70 %**<br>of all samples") 
# 
# data607070 <- rbind(data70, data6070)
# 
# 
# data607070 <- data607070 %>% mutate(State = stringr::str_to_title(state))

ggplot(data = data607070, 
       aes(y = forcats::fct_rev(ASV_nr_taxon), #forcats::fct_reorder(family, -fam_order), 
           x = forcats::fct_reorder(state_colony_frgmt_ms, as.numeric(levels_colony_ms))  
       )) +
  geom_point(aes(#size = rel_abund_bysample,
    size = ifelse(rel_abund_bysample==0, NA, rel_abund_bysample),
    color = colony_ms), 
    shape = 16) + #, alpha = 0.21) +
  # geom_point(aes(#size = rel_abund_bysample, 
  #   size = ifelse(rel_abund_bysample==0, NA, rel_abund_bysample),
  #   shape = State, 
  #   color = colony_ms)) +
  # scale_shape_manual(breaks = c("Symbiotic", "Bleached"), 
  #                    values = c(16, 21)) + 
  scale_color_manual(values = palette_GfasMS) + 
  scale_size_area(name = "Relative abundance",
                  max_size = 10, # 12
                  # n.breaks = 6,
                  breaks = c(0, 1, 10, 20, 30),
                  labels = scales::label_percent(
                    scale = 1,
                    accuracy = 1,
                    suffix = " %")
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 0, size = 9),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9),
    strip.background.x = element_rect(fill = "NA", color = "NA"), # 
    strip.background.y = element_rect(fill = "NA", color = "black"),
    strip.placement = "outside",
    strip.text.y = ggtext::element_markdown(size = 8),
    strip.text.x = element_text(size = 11, vjust = 3, face = "bold"),
    legend.position = "bottom"
  ) +
  guides(shape = "none", color = "none") +
  # facet_wrap(~forcats::fct_rev(State) + forcats::fct_reorder(colony_ms, as.numeric(levels_colony_ms)), 
  #            scales = "free_x",
  #            nrow = 1, strip.position = "top")
  facet_grid(cutoff ~ forcats::fct_rev(State), #+forcats::fct_reorder(colony_ms, as.numeric(levels_colony_ms)),
             scales = "free", space = "free")


ggsave("./out/Gfas_16S/core_mb/nonrarefied/bubble_core_relab_bystate_notransparency.png",
       bg = "white",
       dpi = 330,
       units = "cm", width = 25, height = 9)

ggsave("./out/Gfas_16S/core_mb/nonrarefied/bubble_core_relab_bystate_notransparency.svg",
       # bg = "white",
       # dpi = 330,
       units = "cm", width = 25, height = 9)







# Plot: HEATMAP of CORE microbiome - RELATIVE ABUNDANCE (percent) -------------------------------------------

# To have square tiles and panels by state (and also colony maybe?)
#  maybe necessary to make each plot separateldy ...

## Make subsets for each block -------------------

#### Top block ---------------
rss70 <- data70 %>% 
  filter(origin_ms == "Red Sea" & state == "symbiotic") 

rsb70 <- data70 %>% 
  filter(origin_ms == "Red Sea" & state == "bleached") 

hks70 <- data70 %>% 
  filter(origin_ms == "Hong Kong" & state == "symbiotic") 

hkb70 <- data70 %>% 
  filter(origin_ms == "Hong Kong" & state == "bleached") 

#### Bottom block ---------------
rss6070 <- data6070 %>% 
  filter(origin_ms == "Red Sea" & state == "symbiotic") 

rsb6070 <- data6070 %>% 
  filter(origin_ms == "Red Sea" & state == "bleached") 

hks6070 <- data6070 %>% 
  filter(origin_ms == "Hong Kong" & state == "symbiotic") 

hkb6070 <- data6070 %>% 
  filter(origin_ms == "Hong Kong" & state == "bleached") 






## Function for plot block (DRY!) -----------------

heatmaps <- function(data, title, 
                     x_text = element_text(angle = 90, hjust = 0, size = 5),
                     y_text = element_text(vjust = 0.5, size = 5)) {
  
  ggplot(data = data, 
         aes(x = state_colony_frgmt_ms,
             y = forcats::fct_rev(ASV_nr_taxon))) +
    geom_tile(aes(fill = rel_abund_bysample)
              ) +
    coord_equal() +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "black",
                        limits = c(0, upper_lim),
                        labels = scales::label_percent(scale = 1)) +
    labs(title = title) +
    theme(
      axis.text.x = x_text,
      axis.text.y = y_text,
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5), 
      legend.position = "bottom",
      legend.key = element_rect(color="black"),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      panel.border = element_rect(colour = "black", fill = NA)
    ) +
    guides(
      fill = guide_legend(
        title = "Relative abundance",
        title.position = "left", #  "top",
        title.vjust = 0.55,
        title.hjust = 0,
        label.position = "right",
        nrow = 1)
    )
}

(RSS70 <- heatmaps(data = rss70, y_text = element_blank(), title = "Red Sea\nsymbiotic") )
RSB70 <- heatmaps(data = rsb70, y_text = element_blank(), title = "Red Sea\nbleached")
HKS70 <- heatmaps(data = hks70, y_text = element_blank(), title = "Hong Kong\nsymbiotic")
HKB70 <- heatmaps(data = hkb70, y_text = element_text(), title = "Hong Kong\nbleached")

RSS6070 <- heatmaps(data = rss6070, x_text = element_blank(),
                    y_text = element_blank(), title = element_blank())
RSB6070 <- heatmaps(data = rsb6070, x_text = element_blank(),
                    y_text = element_blank(), title = element_blank())
HKS6070 <- heatmaps(data = hks6070, x_text = element_blank(),
                    y_text = element_blank(), title = element_blank())
HKB6070 <- heatmaps(data = hkb6070, x_text = element_blank(),
                    y_text = element_text(), title = element_blank())


### Patch together ------------

library(patchwork)

layout <- "
ABCD
EFGH
"

(
  pwheat <- (RSS70 + RSB70 + HKS70 + HKB70 + RSS6070 + RSB6070 + HKS6070 + HKB6070) + 
    plot_layout(
      guides = "collect",
      design = layout
      ) &
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5),
      legend.key = element_rect(color = "#000000"),
      legend.key.size = unit(0.35, 'cm'))
)


### Export .PNG ------------------
ggsave("./out/Gfas_16S/core_mb/nonrarefied/heat_core_relab_patch.png",
       bg = "white",
       dpi = 330, 
       units = "cm", 
       width = 18, height = 7)


### Export .SVG ------------------

# install.packages("svglite")
library(svglite)

ggsave("./out/Gfas_16S/core_mb/nonrarefied/heat_core_relab_patch.svg",
       bg = "white",
       # dpi = 330, 
       units = "cm", 
       width = 18, height = 7)


# Finish off details in Inkscape ;)
