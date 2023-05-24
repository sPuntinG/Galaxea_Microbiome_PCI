
# Look at the data produced by SymPortal.

# Data downloaded from here: https://symportal.org/data_explorer/
# should become publicly available in a year from date of data submission
# (21.01.2022 -> 21.01.2023)


library(tidyverse)
library(here)
library(ggtext)


# Output folder ------------------------------------------

# Import SymPortal data - ABSOLUTE abundance ------------------------------------------
data_absolute <- read_delim(file = "./in/SymPortal_downloaded_data/20220121_puntin/post_med_seqs/191_20220121_03_DBV_20220121T073608.seqs.absolute.abund_and_meta.txt", 
                            delim = "\t")

# Keep only relevant info
data_absolute <- data_absolute %>% 
  select(sample_name, c(A1:D5m))


# Get nr of counts per sample - pivot_longer()
absolute_abundance <- data_absolute %>% 
  pivot_longer(
    cols = -sample_name,
    names_to = "absolute_abundance",
    values_to = "nr_reads"
  ) %>% 
  group_by(sample_name) %>% 
  summarize(
    abs_reads_nr = sum(nr_reads)
  ) %>% 
  drop_na() # %>% view()


rm(data_absolute)


# Import SymPortal data - RELATIVE abundance ------------------------------------------

# File with post-MED and relative abundance + meta
data <- read_delim(file = "./in/SymPortal_downloaded_data/20220121_puntin/post_med_seqs/191_20220121_03_DBV_20220121T073608.seqs.relative.abund_and_meta.txt", 
           delim = "\t")


# Create table for clade assignment (to join later to data) ---------------------------------------

ITS_seq_type <- data %>% select(c(A1:D5m)) %>% names()

clades <- tibble(ITS_seq_type) %>% 
  mutate(ITS_clade = case_when(
    str_detect(ITS_seq_type, "A") ~ "clade_A",
    str_detect(ITS_seq_type, "B") ~ "clade_B",
    str_detect(ITS_seq_type, "C") ~ "clade_C",
    str_detect(ITS_seq_type, "D") ~ "clade_D"
  )) 

rm(ITS_seq_type)


# Data re-shaping ---------------------------------------

# Keep only relevant info
data <- data %>% 
  select(sample_name, c(A1:D5m))


# Add metadata info to match MS naming and grouping system
meta <- read_csv("./in/metadata_ms.csv")

# data: "sample_name"
# meta: "new_name"
data <- meta %>% 
  select(new_name, state, colony_ms) %>% 
  rename(sample_name = new_name) %>% 
  inner_join(., data, by = "sample_name") # %>% view()


# Make long
data_long <- data %>% 
  pivot_longer(
    cols = A1:D5m,
    names_to = "ITS_seq_type",
    values_to = "rel_abund"
  )

# Remove zeroes (when ITS seq type is not present)
# which makes table unnecessarily long
data_long <- data_long %>% 
  filter(rel_abund != 0)  # drops ~3500 empty rows


# Add clade info
data_long <- data_long %>% 
  inner_join(., clades, by = "ITS_seq_type") %>% 
  relocate(ITS_clade, .after = ITS_seq_type)


# Add ABSOLUTE abundance (nr of reads per sample) -------------------------------------
data_long <- data_long %>% 
  inner_join(., absolute_abundance, by = "sample_name") %>% 
  relocate(abs_reads_nr, .after = sample_name) #%>% view()


# Rename better - match MS sample nomenclature
data_long <- data_long %>% 
  mutate(
    replicate_id = str_extract_all(sample_name, "\\d+$"),
    sample_name = paste(colony_ms, replicate_id, sep = "_")
  ) #%>% view()




# Clean up
rm(list = setdiff(ls(), "data_long"))



# Plot - see if matches with Symportal (sanity check) -------------------------
ggplot(data_long, aes(x = sample_name, y = rel_abund, fill = ITS_seq_type)) +
  geom_col(color = "#000000") + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Cool but too many ITS_seq_types to be able to make sense of this ...
# To do:
# 0. Keep only Symbiotic samples = what we want to discuss in the ms
# 1. Calculate rel_abund by clade 
# 2. identify the most abundant ones
# 3. pool all the rest as "Others"
# 4. plot again


# 0. Keep only Symbiotic samples = what we want to discuss in the ms
data_symb <- data_long %>% 
  filter(state == "symbiotic") 


# 1.0. Calculate rel_abund by clade -----------------------------------
rel_abund_clade <- data_symb %>% 
  group_by(sample_name, abs_reads_nr, ITS_clade) %>% 
  summarise(rel_abund_clade = sum(rel_abund)) 

# Write summary table
write_csv(rel_abund_clade, "./out/Gfas_ITS_SymPortal/rel_abund_clade.csv")


# 1.1. Plot rel_abund by clade, facet by COLONY -----------------------------------

# Create palette for clades
palette3_clades <- c(
  "clade_C" = "royalblue", # limegreen",
  "clade_D" = "orangered3",
  "clade_A" = "gold"
)



# Same plot as above but faceted by colony
rel_abund_clade %>% 
  mutate(colony = str_extract(sample_name, "^[:graph:]{1,3}")) %>% 
  ggplot(., aes(x = sample_name, y = rel_abund_clade, 
                              fill = forcats::fct_shift(factor(ITS_clade), -1))) +
  geom_col(color = "#000000",
           width = 1) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = palette3_clades) +
  labs(y = "Relative abundance") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = NA, color = NA),
    strip.text.x = element_text(face = "bold"),
    strip.placement = "outside"
  ) +
  facet_grid(~colony, 
             scales = "free", space = "free_x"#, switch = "x"
             )

ggsave("./out/Gfas_ITS_SymPortal/rel_abund_clade_bycolony.png", 
       bg = "white",
       dpi = 310, 
       units = "cm", width = 12, height = 12)



# 2.0. Identify most abundant ITS_seq_types -----------------------------------

# Use this to check for each clade (comment out the others)
data_symb %>% 
  # filter(ITS_clade == "clade_A") %>% 
  filter(ITS_clade == "clade_C") %>%
  # filter(ITS_clade == "clade_D") %>%
  arrange(desc(rel_abund)) # %>% view()


# A way to visualize the types with higher rel_abund 
ggplot(data_symb, aes(x = sample_name, y = rel_abund, color = ITS_clade)) +
  geom_text(aes(label = ITS_seq_type), 
            size = 3, 
            position = position_jitterdodge(
              dodge.width = 0.5, 
              jitter.width = 0.3,
              seed = 21),
            show.legend = F) +
  scale_color_manual(values = palette3_clades) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "grey96"),
    panel.grid.minor.y = element_blank()
  )

ggsave("./out/Gfas_ITS_SymPortal/viz_rel_abund_type.png", bg = "white",
       dpi = 310, 
       units = "cm", width = 15, height = 9)



# 2.1. Plot most abundant ITS_seq_types -----------------------------------

# Only want to focus on the most abundant types (by clade)

# Set cutoff for rel_abund: smaller values are lumped as "others"
cutoff <- 0.05 # only >5 % by sample, else lumped as "others"

data_symb_top <- data_symb %>% 
  mutate(ITS_seq_type = if_else(
    rel_abund > cutoff,  
    ITS_seq_type, 
    paste("other", ITS_clade, sep = "_")) ) %>% 
  group_by(sample_name, abs_reads_nr, ITS_seq_type) %>% 
  summarise(
    rel_abund = sum(rel_abund) ) %>%  
  arrange(sample_name, desc(rel_abund)) %>%
  ungroup() %>% 
  mutate(ITS_seq_type = str_replace_all(ITS_seq_type, 
                                        c("_clade_C" = " *Cladocopium*", 
                                          "_clade_D" = " *Durusdinium*", 
                                          "_clade_A" = " *Symbiodinium*"))) #%>% view()




# Write table for SM - note that name is based on set cutoff value! B)
write_csv(data_symb_top, paste0("./out/Gfas_ITS_SymPortal/rel_abund_clade_top", cutoff * 100, "perc.csv") )


# Make palette for plot
palette_5blues <- c("#B7F6FF", 
                    "#5DD8F1", 
                    "#09BFFF", 
                    "#4472C4", 
                    "#3B3A72" ) 



palette12 <- c(
  "C1"  = palette_5blues[1],  
  "other *Cladocopium*" = palette_5blues[5], 
  "C1c" = palette_5blues[2], 
  "D1" = "orangered2",
  "C1b" = palette_5blues[3], 
  "C41f" =  palette_5blues[4], 
  "C41" =  palette_5blues[4], 
  "D4" = "orangered3",
  "A1" = "gold",
  "C39" =  palette_5blues[4], 
  "other *Durusdinium*" = "orangered4",
  "other *Symbiodinium*" = "gold4"
)




# Facet by COLONY -----------------------------------

toplot <- data_symb_top %>% 
  mutate(colony = str_extract(sample_name, "^[:graph:]{1,3}")) %>%
  group_by(sample_name) %>%
  mutate(label = ifelse(rel_abund == max(rel_abund), abs_reads_nr, NA)) %>% 
  ungroup() %>% 
  arrange(desc(rel_abund)) %>% 
  mutate(order = row_number()) %>% 
  mutate(order2 = ifelse(str_detect(ITS_seq_type, "Cladocopium"), order + 100, order)) %>% # 100 is just a big number
  view()
  
# ggplot(toplot, 
#        aes(x = sample_name, y = rel_abund, # rel_abund,
#            fill = forcats::fct_reorder(ITS_seq_type, -order2) # rel_abund
#            )) +
#   geom_col(color = "#000000", width = 1) + 
#   geom_text(aes(label = label, y = 0.13), 
#             angle = 90,
#             hjust = "right",
#             size = 3
#             ) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values = palette12,
#                     breaks = c("C1", "C1c",  
#                                "D1" ,"C1b", "C41f", "C41", "D4", "A1", "C39", 
#                                "other *Durusdinium*", "other *Symbiodinium*",
#                                "other *Cladocopium*")) + #, 
#                     # guide = guide_legend(reverse = F)) +
#   labs(y = "Relative abundance",
#        fill = "ITS2 sequences") +
#   theme_classic() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
#     axis.title.y = element_text(margin = margin(r = 9)),
#     legend.key.size = unit(0.4, 'cm'), #change legend key size
#     legend.title = element_text(size = 9), #change legend title font size
#     legend.text = ggtext::element_markdown(size = 8),
#     axis.ticks.x = element_blank(),
#     strip.background = element_rect(fill = NA, color = NA),
#     strip.text.x = element_text(face = "bold"),
#     strip.placement = "outside"
#   ) +
#   facet_grid(~colony, 
#              scales = "free", space = "free_x"#, switch = "x"
#   )
# 
# 
# ggsave("./out/Gfas_ITS_SymPortal/rel_abund_clade_top5perc_faceted.png", 
#        bg = "111111",
#        dpi = 330,
#        units = "cm", width = 15, height = 9)



### Plot with legend order (by "clade") --------------------------------------
ggplot(toplot, 
       aes(x = sample_name, y = rel_abund, # rel_abund,
           fill = forcats::fct_reorder(ITS_seq_type, -order2) # rel_abund
       )) +
  geom_col(color = "#000000", width = 1) + 
  geom_text(aes(label = label, y = 0.13), 
            angle = 90,
            hjust = "right",
            size = 3
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = palette12,
                    breaks = c("A1", "other *Symbiodinium*",
                               "D1", "D4",   
                               "other *Durusdinium*",
                               "C1", "C1c",  
                               "C1b", "C39", "C41f", "C41",
                               "other *Cladocopium*")) + #, 
  # guide = guide_legend(reverse = F)) +
  labs(y = "Relative abundance",
       fill = "ITS2 sequences") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 9)),
    legend.key.size = unit(0.4, 'cm'), #change legend key size
    legend.title = element_text(size = 9), #change legend title font size
    legend.text = ggtext::element_markdown(size = 8),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text.x = element_text(face = "bold"),
    strip.placement = "outside"
  ) +
  facet_grid(~colony, 
             scales = "free", space = "free_x"#, switch = "x"
  )


ggsave("./out/Gfas_ITS_SymPortal/rel_abund_clade_top5perc_faceted.png", 
       bg = "111111",
       dpi = 330,
       units = "cm", width = 15, height = 9)




