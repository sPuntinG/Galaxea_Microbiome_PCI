
library(tidyverse)
library(here)
# library(phyloseq)

ps_rarefied <- readRDS("./out/RDS_files/ps_rarefied.rds")

metadata_ms <- read_csv("./in/metadata_ms.csv")

source("colors_GfasMS.R")


# Observed, Chao, ACE, Shannon, Simpson, InvSimpson, Fisher -------------------------

alphas <- phyloseq::estimate_richness(ps_rarefied) %>% # ps_rarefied
  as_tibble(rownames = "sample_id")  %>% 
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-")) %>% 
  inner_join(., metadata_ms, by = "sample_id") %>% 
  select(-c("new_name", "sample_type", "origin")) %>% 
  relocate(c("colony", "state", "origin_ms", "colony_ms"), .after = sample_id)

# So now I have all these alpha diversity indices
# ("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
# in one tibble ('alphas"), by sample + metadata.



# Pielou's evenness ------------------------------

# Calculate Pielou's evenness - requires {microbiome} ..............................
library(microbiome)

pielous <- microbiome::evenness(ps_rarefied, 'pielou') %>% 
  as_tibble(rownames = "sample_id") %>% 
  rename(Pielou = pielou)

# Add Pileou's (evenness) index to the table
alphas <- alphas %>% 
  left_join(., pielous, by = "sample_id") #%>% view()


# Faith phylogenitic distance - requires {btools} -----------------------------------

# install.packages("remotes")
# remotes::install_github("twbattaglia/btools")
library(btools)

faith <- btools::estimate_pd(ps_rarefied) %>% 
  as_tibble(rownames = "sample_id") %>% 
  rename(
    Faith_pd = PD) %>% 
  select(-SR) # sample richness, already have in table

alphas <- alphas %>% 
  left_join(., faith, by = "sample_id") # %>% view()


# All ALPHA's table for supplementary ---------------------
# Make a nice and tidy
alphas <- alphas %>% 
  mutate(origin_ms = factor(origin_ms, levels = c("Red Sea", "Hong Kong"))) %>% 
  mutate(state = factor(state, levels = c("symbiotic", "bleached"))) %>%
  select(sample_id, state_colony_frgmt_ms, origin_ms, colony_ms, state,
         Observed, Pielou, Chao1, Shannon, Simpson, InvSimpson, Fisher, Faith_pd) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  arrange(origin_ms, colony_ms, state) %>% 
 write_csv(., "./out/Gfas_16S/useful_tables/rarefied_alpha_div_metrics.csv")


## Export as RDS ---------------
saveRDS(alphas, "./out/RDS_files/rarefied_alpha_div_metrics.rds")


# All ALPHA's - by group ---------------------

# Make a 3 summary, after grouping by:
# 1. state
# 2. state, origin
# 3. state, origin, colony

# 1. state
alphas %>% 
  group_by(state) %>% 
  summarise(
    n = length(sample_id),
    Observed = mean(Observed), 
    Chao1  = mean(Chao1), 
    Shannon = mean(Shannon), 
    Simpson = mean(Simpson), 
    InvSimpson = mean(InvSimpson), 
    Fisher = mean(Fisher)
  ) %>% 
  ungroup()


# 2. state, origin
alphas %>% 
  group_by(state, origin_ms) %>% 
  summarise(
    n = length(sample_id),
    Observed = mean(Observed), 
    Chao1  = mean(Chao1), 
    Shannon = mean(Shannon), 
    Simpson = mean(Simpson), 
    InvSimpson = mean(InvSimpson), 
    Fisher = mean(Fisher)
  ) %>% 
  ungroup()


# 3. state, origin, colony
alphas %>% 
  group_by(state, origin_ms, colony_ms) %>% 
  summarise(
    n = length(sample_id),
    Observed = mean(Observed), 
    Chao1  = mean(Chao1), 
    Shannon = mean(Shannon), 
    Simpson = mean(Simpson), 
    InvSimpson = mean(InvSimpson), 
    Fisher = mean(Fisher)
  ) %>% 
  ungroup() # %>% write_csv(., "./R_outputs/alphas_cleaned_rarefied_grouped.csv")


# PLOT ALPHA diversities by STATE (x axis) and ORIGIN --------------------------------------

# Function to automate plotting to explore all alpha diversity indices
alpha_plots <- function(data, y) {
  y <- enquo(y)
  
  ggplot(data = data, aes(x = state, y = !!y)) +
    geom_boxplot(aes(group = state),
                 outlier.shape = NA,
                 width = 0.4,
                 color = "gray33") + 
    ggbeeswarm::geom_quasirandom(
      aes(color = colony_ms, shape = state, size = state), # group = origin_ms, 
      width = 0.2,
      stroke = 1.3  ) +
    scale_color_manual(values = palette_GfasMS) +
    scale_shape_manual(values = c(16, 1), breaks = c("symbiotic", "bleached")) +
    scale_size_manual(values = c(2, 1.5), breaks = c("symbiotic", "bleached")) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.ticks.x = element_blank(),
      text = element_text(size = 13),
      panel.grid = element_blank(),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.spacing = unit(0.05, 'cm'),
      legend.key.size = unit(0.3, 'cm'),
      legend.justification = "top",
      aspect.ratio = 3,
      plot.margin = margin(t = 0, r = 0, l = 0.7, b = 0, "cm")
    )
}

alpha_plots(data = alphas, y = Faith_pd)

obs <- alpha_plots(data = alphas, y = Observed)
chao <- alpha_plots(data = alphas, y = Chao1) 
piel <- alpha_plots(data = alphas, y = Pielou) 
shan <- alpha_plots(data = alphas, y = Shannon)
simp <- alpha_plots(data = alphas, y = Simpson)  
# invsimp <- alpha_plots(data = alphas_cleaned, y = InvSimpson)
# fish <- alpha_plots(data = alphas_cleaned, y = Fisher)
fait <- alpha_plots(data = alphas, y = Faith_pd)
  
library(patchwork)
patch <- (obs + chao + piel) / (shan + simp + fait)

out <- patch + plot_layout(guides = "collect")

ggsave("./out/Gfas_16S/alpha_diversity/rarefied_alpha_div_metrics_by_state.png", 
       bg = "white",
       dpi = 330,
       units = "cm", width = 23, height = 20)


# After plotting each alpha diversity it does not look like 
# there's any significant difference between symb and bleached


# Plot for MS: Obs, Shannon, Simpson ----------------------------------------

ms <- (obs + shan + simp)

ms_out <- ms + 
  plot_annotation(
    theme = theme(
      legend.justification = "top")) +
  plot_layout(guides = "collect")


ggsave("./out/Gfas_16S/alpha_diversity/rarefied_alpha_MS.png", 
       bg = "white",
       dpi = 330,
       units = "cm", width = 19, height = 13)



# Stat tests -----------------------
# Test difference is ALPHA diversity indexes between Symb and Bleached

library("ggpubr")

# Visual inspection of normality of distribution
# All good except 3 ... (Pielou, Shannon and Simpson)
obs <- ggpubr::ggqqplot(alphas, x = "Observed", title = "Observed")
chao <- ggpubr::ggqqplot(alphas, x = "Chao1", title = "Chao1")
piel <- ggpubr::ggqqplot(alphas, x = "Pielou", title = "Pielou")   # very bad
shan <- ggpubr::ggqqplot(alphas, x = "Shannon", title = "Shannon")  # very bad
simp <- ggpubr::ggqqplot(alphas, x = "Simpson", title = "Simpson")  # worst
fait <- ggpubr::ggqqplot(alphas, x = "Faith_pd", title = "Faith_pd")


patch <- (obs + chao + piel) / (shan + simp + fait)

out <- patch + plot_layout(guides = "collect") +
  plot_annotation(
    title = 'Inspect normality of distribution',
    theme = theme(
      plot.title = element_text(size = 20, face = "bold")) )

ggsave("./out/Gfas_16S/alpha_diversity/rarefied_alpha_div_metrics_ggqqplot.png", 
       bg = "white",
       dpi = 330,
       units = "cm", width = 27, height = 20) 



# Homogeneity of variance - F-test ......................................

var.test(data = alphas, Observed ~ state) # 0.6128
var.test(data = alphas, Chao1 ~ state)    # 0.6315
var.test(data = alphas, Pielou ~ state)   # 0.1031
var.test(data = alphas, Shannon ~ state)  # 0.1401
var.test(data = alphas, Simpson ~ state)  # 0.03601 *
var.test(data = alphas, Faith_pd ~ state) # 0.3657


# Homogeneity of variance - BARTLETT test ......................................

# H0 = the variances in each of the groups (samples) are the same
# -> want p.value > 0.05 -> can proceed with t-test

# values ~ groups
bartlett.test(data = alphas, Observed ~ state) # 0.6181
bartlett.test(data = alphas, Chao1 ~ state)    # 0.6368
bartlett.test(data = alphas, Pielou ~ state)   # 0.1047
bartlett.test(data = alphas, Shannon ~ state)  # 0.1421
bartlett.test(data = alphas, Simpson ~ state)  # 0.03673 *
bartlett.test(data = alphas, Faith_pd ~ state) # 0.3696

# Same results as with the F-test (above)
# can do regular t-test for all but NOT SIMPSON (will need Welch test)


# t-test ----------------------------------------------------------------------

# t-test is used when the population has n < 30

library(rstatix)
rstatix::t_test(data = alphas, Observed ~ state,
                alternative = "less", # "two.sided"
                paired = F,
                var.equal = T)                        # 0.086 # 0.172  
          t.test(data = alphas, Observed ~ state,
                alternative = "less", # "two.sided"
                paired = F,
                var.equal = T) # t = -1.4059, df = 25, p-value = 0.08603

rstatix::t_test(data = alphas, Chao1 ~ state,
                alternative = "less", # "two.sided"
                paired = F,
                var.equal = T)                        # 0.0851 # 0.170
rstatix::t_test(data = alphas, Pielou ~ state,
                alternative = "less", # "two.sided"
                paired = F,
                var.equal = T)                        # 0.0767 # 0.153
rstatix::t_test(data = alphas, Shannon ~ state,
                alternative = "less", # "two.sided"
                paired = F,
                var.equal = T)                        # 0.0632 # 0.126
rstatix::t_test(data = alphas, Simpson ~ state,
                alternative = "less", # "two.sided"              
                paired = F,
                var.equal = FALSE)                    # 0.0854 # 0.171 ***unequal variance here!***
rstatix::t_test(data = alphas, Faith_pd ~ state,
                alternative = "less", # "two.sided"                 
                paired = F,
                var.equal = T)                        # 0.0826 # 0.165


# None of these is significant
# Not even when I do it one-sided (= group 1 is "less" than group 2), 
# which should be more sensitive, and would make sense to ask
# if symbiotic samples have smaller (alpha diversity) that bleached
# (based on the Anna Karenina Pr. (AKP) of microbial ecology)


# Wilcoxon test (Mann-Whitney) ----------------------------------------------------------------------

# Non-parametric 
# necessary for Pielou, Shannon, and Simpson

# Pielou
stats::wilcox.test(data = alphas, 
                   Pielou ~ state,
                   alternative = "two.sided") # W = 65, p-value = 0.2152

# Shannon
stats::wilcox.test(data = alphas, 
                   Shannon ~ state,
                   alternative = "two.sided") # W = 66.5, p-value = 0.244

# Simpson
stats::wilcox.test(data = alphas, 
                   Simpson ~ state,
                   alternative = "two.sided") # W = 66, p-value = 0.2318 

 
# Again, NOT significant difference between symbiotic states
# for all alpha diversity metrics considered 

# alternative hypothesis: true location shift is not equal to 0
# is not rejected (p > 0.05)

# Note that p.value is halved if 'alternative = "less"' but 
# is still always above 0.05
