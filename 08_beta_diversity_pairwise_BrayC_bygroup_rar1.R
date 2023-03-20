
##########################################################################
# Alternative to PERMDISP and PERMPANOVA to look at within group dispersion:
# plot pairwise distance by group
#  with group = colony * state
##########################################################################

library(tidyverse)
library(here)
# library(phyloseq)

dist_bray <- readRDS("./out/RDS_files/dist_bray.rds")
metadata_ms <- read_csv("./in/metadata_ms.csv")

source("colors_GfasMS.R")

library("ggtext")

# Transform distance matrix into pairwise list -----------------------------

# Make 'dist' object into data.frame
dist_bray_df <- as.data.frame(as.matrix(dist_bray))

# Make tibble (tidy object)
dist_bray_tbl <- dist_bray_df %>% as_tibble(rownames = "sample_id_L")

# Pivot longer (so to have 3 cols)
dist_bray_long <- dist_bray_tbl %>% 
  pivot_longer(
    -sample_id_L,
    names_to = "sample_id_R",
    values_to = "pairwise_distance"
  ) 


# Rename samples to match ms format 
# 1. create two tibbles, one for each col to rename
siL <- metadata_ms %>% 
  select(sample_id, state_colony_frgmt_ms) %>% 
  rename(sample_id_L = sample_id)


siR <- metadata_ms %>% 
  select(sample_id, state_colony_frgmt_ms) %>% 
  rename(sample_id_R = sample_id)

# 2. join each to dist_bray_long separately
dist_bray_paired <- inner_join(dist_bray_long, siL, by = "sample_id_L") %>% 
  rename(sample_L = state_colony_frgmt_ms) %>% 
  inner_join(., siR, by = "sample_id_R") %>% 
  rename(sample_R = state_colony_frgmt_ms) %>% 
  select(sample_L, sample_R, pairwise_distance) #%>% view()

# So now the samples have been properly renamed 
# But we still have too many pairwise:
length(dist_bray_paired$pairwise_distance) # 729
# because we still have all the same-samples comparison (distance = 0)
# and the "mirror" comparison (top and bottom triangle of original matrix)

# What we should have   n(n - 1)/2 pairwise comparisons
#                   = 27(27 - 1)/2 pairwise comparisons
# (= all unique pairwise comparison among 27 samples)
n <- dist_bray_tbl$sample_id_L %>% length()
n * (n - 1)/2 # 351 unique pairwise comparisons

# Get rid of duplicated and sample-sample pairs
dist_bray_paired <- dist_bray_paired %>% 
  filter(sample_L < sample_R)  

dist_bray_paired$pairwise_distance %>% length() # 351 correct!

# Now need to select only the sample pairs that belong to the same group
# as state x colony
dist_bray_paired_grouped <- dist_bray_paired %>% 
  mutate(
    gr_L = str_extract(sample_L, "^\\w_\\w\\w\\d"),
    gr_R = str_extract(sample_R, "^\\w_\\w\\w\\d")
  ) %>% 
  filter(gr_L == gr_R) # %>% view()

# Add state and colony_ms (for plotting)
dist_bray_paired_grouped <- dist_bray_paired_grouped %>% 
  separate(gr_L, sep = "_", into = c("state", "colony_ms")) %>% 
  rename(state_colony = gr_R) # same as gr_L


# Identify groups with n = 3
n3 <- dist_bray_paired_grouped %>% 
  group_by(state_colony) %>% # could also be 'gr_R' as they're exactly the same
  count() %>% 
  filter(n == 3) %>% #view()
  pull(state_colony)


# N3 --------------------------

## Exportable N3 tibble ----------------

pwcBD_n3_tibble <- dist_bray_paired_grouped %>%
  filter(state_colony %in% n3)

saveRDS(pwcBD_n3_tibble, "./out/RDS_files/pwcBD_n3_tibble.rds")



## Test difference between groups (only groups with n = 3) -------------------------------------

symb <- dist_bray_paired_grouped %>% 
  filter(state_colony %in% n3) %>% 
  filter(state == "S")

blea <- dist_bray_paired_grouped %>% 
  filter(state_colony %in% n3) %>% 
  filter(state == "B")

### Normality check ---------
# (want p-value > 0.05)
shapiro.test(symb$pairwise_distance) # p-value = 0.3086  OK
shapiro.test(blea$pairwise_distance) # p-value = 0.04454 NOT OK

# since data is NOT normally distributed canNOT use F-test (assumes normality)
# need to use 'Mann-Whitney U test' aka 'Wilcoxon test'.

### Wilcoxon rank sum exact test ----------
w_n3 <- dist_bray_paired_grouped %>% 
  filter(state_colony %in% n3) %>% 
  stats::wilcox.test(data = ., pairwise_distance ~ state,
              alternative = "greater", # test if pairwise_distance is greater among bleached than symbiotic
              paired = F)

w_n3$p.value   # P = 0.0001020651 **
w_n3$statistic # W = 102
w_n3[1:7]

# Very significant! P < 0.001

saveRDS(w_n3, "./out/RDS_files/w_n3.rds")



## Plot by group and only groups with n = 3 - Narrow -------------------------------------
(
  pwcBD_n3 <- dist_bray_paired_grouped %>%
    filter(state_colony %in% n3) %>% 
    ggplot(., aes(
      x = forcats::fct_rev(state), 
      y = pairwise_distance,
      color = colony_ms,
      shape = state)) +
    geom_boxplot(show.legend = F, color = "grey69", outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(aes(size = state),
                                 width = 0.2,
                                 stroke = 1.5) +
    scale_color_manual(values = palette_GfasMS) +
    scale_shape_manual(values = c(16, 1), #16, 21
                       breaks = c("S", "B"),
                       labels = c("Symbiotic", "Bleached")) + 
    scale_size_manual(values = c(1.75, 1), 
                      breaks = c("S", "B")) +
    scale_x_discrete(breaks = c("S", "B"), labels = c("Symbiotic", "Bleached")) +
    ylab("Bray-Curtis pairwise distance") +
    labs(caption = paste0(w_n3$method, ",<br> *P* = ",
                           format(round(w_n3$p.value, 4), scientific = F),
                           ", *W* = ", w_n3$statistic),
         color = "Colony") +
    theme_classic() +
    theme(
      text = element_text(family = "pfont"),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      plot.caption = element_markdown(size = 7, hjust = 1),
      legend.position = "none",
      aspect.ratio = 3
    )
)  
  # guides(state = "none", size = "none", shape = "none")

ggsave("./out/Gfas_16S/supplementary/bray_pairwisedist_bygroup_n3_narrow.png", 
       dpi = 330, 
       bg = "white",
       units = "cm", width = 6, height = 10)

saveRDS(pwcBD_n3, "./out/RDS_files/pwcBD_n3_ggplot.rds")



# Nx --------- 
## Test difference between groups (all groups, including when n = 1) -------------------------------------

symb <- dist_bray_paired_grouped %>% 
  # filter(state_colony %in% n3) %>% # keep all groups now!
  filter(state == "S")

blea <- dist_bray_paired_grouped %>% 
  # filter(state_colony %in% n3) %>% # keep all groups now!
  filter(state == "B")

# Normality check
# (want p-value > 0.05)
shapiro.test(symb$pairwise_distance) # p-value = 0.4873  OK
shapiro.test(blea$pairwise_distance) # p-value = 0.09702 OK

# since data is normally distributed can use F-test (or also Levene ...)

# Homogeneity of variance check
# (F-test, want p-value > 0.05)
var.test(pairwise_distance ~ state, data = dist_bray_paired_grouped) # p-value = 0.3699 OK


# Two-sample t-test 
tt_n1n3 <- t.test(symb$pairwise_distance, blea$pairwise_distance, 
                alternative = "less",  # test if: mean of symb is LESS than blea
                # alternative = "two.sided",
                var.equal = T) 

tt_n1n3$p.value # 0.0003512971

tt_n1n3[1:10] %>% as.data.frame() 


## Plot by group all (also n = 1) ------------------------------------------------
dist_bray_paired_grouped %>% 
  # filter(state_colony %in% n3) %>% 
  ggplot(., aes(
    x = forcats::fct_rev(state), 
    y = pairwise_distance,
    color = colony_ms,
    shape = state )) +
  geom_boxplot(aes(group = state), color = "grey69", outlier.shape = NA,
               show.legend = F) +
  ggbeeswarm::geom_quasirandom(aes(size = state),
                               width = 0.2,
                               stroke = 1.5) +
  # geom_point(aes(size = state),
  #   position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4, seed = 321),
  #            stroke = 1.5) + # size = 2, 
  scale_color_manual(values = palette_GfasMS) +
  scale_shape_manual(values = c(16, 1), 
                     breaks = c("S", "B"),
                     labels = c("Symbiotic", "Bleached")) + 
  scale_size_manual(values = c(1.75, 1), 
                    breaks = c("S", "B")) +
  scale_x_discrete(breaks = c("S", "B"), labels = c("Symbiotic", "Bleached")) +
  ylab("Bray-Curtis pairwise distance") +
  labs(caption = paste0(tt_n1n3$method, ", **P = ",
                         format(round(tt_n1n3$p.value, 4), scientific = F),
                         "**"),
       color = "Colony") +
  theme_classic() +
  theme(
    text = element_text(family = "pfont"),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.caption = element_markdown(size = 7, hjust = 1),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.spacing = unit(0.05, 'cm'),
    legend.key.size = unit(0.3, 'cm'),
    legend.background = element_rect(size = 0.3, 
                                     linetype = "solid", 
                                     colour = "grey69"),
    legend.justification = "top",
    aspect.ratio = 2
  ) +
  guides(state = "none", size = "none", shape = "none")

ggsave("./out/Gfas_16S/supplementary/bray_pairwisedist_bygroup_n1n3_narrow.png", 
       dpi = 330, 
       bg = "white",
       units = "cm", width = 10, height = 10)



# clean up ---------------------------
rm(list = setdiff(ls(), c("ps", "ps_rarefied", "metadata_ms", "palette_GfasMS"))) 
