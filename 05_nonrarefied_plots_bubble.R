
# Plots for the manuscript (and extra bar plot and bubble plot)

source("04_organize_Family_level.R")


# Color palette for origin
source("colors_GfasMS.R")




# BUBBLE Plot by individual samples (for MS!) ---------------------------------------------

## Prep data -----------------
# Order Families by abundance 
fam_order_tab <- fam %>% 
  group_by(family) %>% 
  summarise(
    tot_abund = sum(rel_abund_bysample)
  ) %>% 
  arrange(desc(tot_abund)) %>% 
  mutate(
    fam_ordered = row_number()) %>% 
  mutate(
    fam_order = if_else(family == "Others", 19, as.numeric(fam_ordered) )
  ) %>% 
  select(family, fam_order) # %>% view()


# Add family order to main tibble 
fam <- fam %>% 
  left_join(., fam_order_tab, by = "family") #%>% view()



## Plot with samples name --------------------------

fam <- fam %>% mutate(State = stringr::str_to_title(state)) 

ggplot(data = fam, 
       aes(y = forcats::fct_reorder(family, -fam_order), 
           x = forcats::fct_reorder(state_colony_frgmt_ms, levels_colony_ms)  
       )) +
  geom_point(aes(size = rel_abund_bysample, color = colony_ms), shape = 16, alpha = 0.21) +
  geom_point(aes(size = rel_abund_bysample, shape = State, color = colony_ms)) +
  scale_shape_manual(breaks = c("Symbiotic", "Bleached"), 
                     values = c(16, 21)) + 
  scale_color_manual(values = palette_GfasMS) + 
  scale_size_area(name = "Relative abundance",
                  max_size = 12,
                  n.breaks = 4,
                  # breaks = c(20, 30, 40),
                  labels = scales::label_percent(
                    accuracy = 1,
                    suffix = " %")
                  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9),
    strip.background = element_rect(fill = "NA", color = "NA"), # 
    strip.placement = "outside",
    strip.text = element_text(size = 11, vjust = 3, face = "bold"),
    legend.position = "bottom"
  ) +
  guides(shape = "none", color = "none") +
  facet_wrap(~forcats::fct_rev(State), 
             scales = "free_x",
             nrow = 1, strip.position = "top")
             

ggsave("./out/Gfas_16S/Bubble_plot/17fams_indivsamples_wnames.png", 
       bg = "white",
       dpi = 330, 
       units = "cm", width = 25, height = 15)

library(svglite)
ggsave("./out/Gfas_16S/Bubble_plot/17fams_indivsamples_wnames.svg", 
       # bg = "white",
       dpi = 330, 
       units = "cm", width = 25, height = 15)

ggsave("./out/Gfas_16S/Bubble_plot/17fams_indivsamples_wnames.eps", 
       units = "cm", width = 25, height = 15, dpi = 600,
       # bg = "white",
       # plot = state_diamond,
       device = cairo_ps, fallback_resolution = 600)



## Plot with samples name NO TRANSAPARENCY --------------------------

fam <- fam %>% mutate(State = stringr::str_to_title(state))

ggplot(data = fam, 
       aes(y = forcats::fct_reorder(family, -fam_order), 
           x = forcats::fct_reorder(state_colony_frgmt_ms, levels_colony_ms)  
       )) +
  geom_point(aes(size = rel_abund_bysample, color = colony_ms), shape = 16) + #, alpha = 0.21) +
  scale_color_manual(values = palette_GfasMS) + 
  scale_size_area(name = "Relative abundance",
                  max_size = 12,
                  n.breaks = 4,
                  # breaks = c(20, 30, 40),
                  labels = scales::label_percent(
                    accuracy = 1,
                    suffix = " %")
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9),
    strip.background = element_rect(fill = "NA", color = "NA"), # 
    strip.placement = "outside",
    strip.text = element_text(size = 11, vjust = 3, face = "bold"),
    legend.position = "bottom"
  ) +
  guides(shape = "none", color = "none") +
  facet_wrap(~forcats::fct_rev(State), 
             scales = "free_x",
             nrow = 1, strip.position = "top")


ggsave("./out/Gfas_16S/Bubble_plot/17fams_indivsamples_wnames_notransrency.png", 
       bg = "white",
       dpi = 330, 
       units = "cm", width = 25, height = 15)

library(svglite)
ggsave("./out/Gfas_16S/Bubble_plot/17fams_indivsamples_wnames_notransrency.svg", 
       dpi = 330, 
       units = "cm", width = 25, height = 15)

ggsave("./out/Gfas_16S/Bubble_plot/17fams_indivsamples_wnames.eps", 
       units = "cm", width = 25, height = 15, dpi = 600,
       device = cairo_ps, fallback_resolution = 600)



# Plot by GROUP (state x colony) ------------------------------------------------------ 

# 1. N by group
n_bygroup <- fam %>% 
  group_by(group) %>% 
  summarise(n = length(unique(sample_id))) %>% 
  ungroup() #%>% view()

# 2. Nr of reads/count by family by each group
fam_count_bygroup <- fam %>% 
  group_by(group, family) %>% 
  summarise(
    count = sum(count_bysample)) %>% 
  ungroup() # %>% view()

# 3. Tot nr of reads/counts by group (all fams pooled together)
tot_count_bygroup <- fam %>% 
  group_by(group) %>% 
  summarise(
    tot_count = sum(count_bysample)) %>% 
  ungroup()# %>% view()

# 4. metadata
meta <- fam %>% 
  select(group, state, origin_ms, colony_ms, levels_colony_ms)


# Put all together
fam_group_plot <- inner_join(n_bygroup, fam_count_bygroup, by = "group") %>% 
  inner_join(., tot_count_bygroup, by = "group") %>% 
  # inner_join(., meta, by = "group") %>%
  mutate(
    rel_abund = count / tot_count,
    rel_abund_perc = rel_abund * 100) # %>% view()
# group_by(group) %>% # sanity check
# summarise(tot_freq = sum(rel_freq)) %>% 
# view() # sanity check passed :D

# Add metadata
fam_group_plot <- left_join(fam_group_plot, meta, by = "group") %>% unique()

# Add Family order
fam_group_plot <- left_join(fam_group_plot, fam_order_tab, by = "family") 

rm(n_bygroup, fam_count_bygroup, tot_count_bygroup, meta, fam_order_tab)


# BAR plot - not included in manuscript ..............................................

# Pelette with 17 fams + "Others" = 18 colors                 
pal18 <- c(
  "Rhodobacteraceae" = "#812599", #
  "Alteromonadaceae" = "gold1",  #
  "Moraxellaceae" = "#75A5F8", #
  "Bacillaceae" = "#EE8041", #
  "unclassif_Gammaproteobacteria" = "red4", #
  "Hyphomonadaceae" = "#00CCCC", #
  "Endozoicomonadaceae" = "#bfef45", #
  "unclassif_Cellvibrionales" = "blue4", #
  "Amoebophilaceae" = "plum", # "thistle1",  #
  "Colwelliaceae" = "#aaffc3",  #
  "Flavobacteriaceae" =  "green", #
  "Micavibrionaceae" = "#FF5FBB",  #
  "Saprospiraceae" = "firebrick3", #
  "Saccharospirillaceae" = "cyan",  #
  "Cryomorphaceae" =  "#2773F6",  #
  "Cyclobacteriaceae" =  "violet", #
  "Vibrionaceae" = "darkseagreen1",
  "Others" = "grey91" )


ggplot(data = fam_group_plot, 
       aes(x = group, 
           y = rel_abund_perc, 
           fill = forcats::fct_reorder(family, -fam_order))) + # forcats::fct_rev(family)
  geom_col(color = "black") +
  scale_y_continuous(expand = c(.01, 0)) +
  scale_fill_manual(name = "Bacterial families", values = pal18) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 15),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 16)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(~forcats::fct_rev(state), 
             nrow = 1, 
             scales = "free_x")

ggsave("./out/Gfas_16S/bubble_plot/bar_17fams.png", dpi = 210, units = "cm", width = 27, height = 21)


