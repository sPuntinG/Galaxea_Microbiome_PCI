
library(tidyverse)
library(here)
library(phyloseq)
# library(btools)

source("colors_GfasMS.R")

# ************
# the following part is VERY COMPUTATIONALLY INTENSIVE (for my pc ... !)
#  -> skip to the next xxxxxxxxx where I re-import the data that I create here 
#     and that I need for the following steps

# xxxxxxxxxxxx ----------


# Need to repeat the iterated rarefaction because the phyloseq object with the 
#  averaged values that I obtained from `metagMisc::phyloseq_mult_raref_avg()` 
#  is not compatible with phyloseq::estimate_richness() because this requires count data 
#  and `metagMisc::phyloseq_mult_raref_avg()` produces relative abudnaces

# Import non-rarefied phyloseq object: "ps" -----------------------------------
ps <- readRDS("./out/RDS_files/ps.rds")

# Iterate rarefaction with for loop ------------------------------------

iterations <- 1000 # Note that I used the same seeds for the ps used for beta diversity!

rarefied <- list("phyloseq")  # 1. output
for (i in 1:iterations) {            # 2. sequence
  rarefied[[i]] <- phyloseq::rarefy_even_depth(physeq = ps, 
                                               rngseed = i, # note: same as previous set.seed()
                                               sample.size = 2690, 
                                               replace = F)      # 3. body
}
rarefied

# Observed, Chao, ACE, Shannon, Simpson, InvSimpson, Fisher -------------------------
# Alphas from phyloseq::estimate_richness() 

## Iterate calculations of alpha diversity over list of rarefied ps ------------------------

# Check that phyloseq function works 
phyloseq::estimate_richness(rarefied[[1]]) # good!

# Another for loop to iterate estimate_richness
phyalphas <- list("tbl_df") # "data.frame")  # 1. output
for (i in 1:length(rarefied) ) {            # 2. sequence
  phyalphas[[i]] <- phyloseq::estimate_richness(rarefied[[i]]);  
  phyalphas[[i]] <- as_tibble(phyalphas[[i]], rownames = "sample_id") # 3. body
}
phyalphas

phyalphas[[1]]
phyalphas[[1]] %>% class()

## Bind all data.frames together to make single (long) table -----------------------

# Check that rbind() is what I need (I always confuse with cbind()!)
rbind(phyalphas[[1]], phyalphas[[2]]) #%>% view() # correct

phyalphas %>% length() # 5
phyalphas[[1]] %>% nrow() # 27 -> all phyalphas_rbind should have 5 * 27 = 135 rows

phyalphas_rbind <- bind_rows(phyalphas) #, .id = "column_label") this if I want to keep the [[i]] info as column

# Sanity check: n must be = as nr of iteration
phyalphas_rbind %>% 
  group_by(sample_id) %>% 
  summarise(n = length(Observed))


## Calculate mean (summarize) -----------------
phyalphas_rbind <- phyalphas_rbind %>% 
  group_by(sample_id) %>% 
  summarise_if(is.numeric, mean) #%>% view()


# Pielou's evenness ------------------------------

# Calculate Pielou's evenness - requires {microbiome} 
# library(microbiome)

# Check that {microbiome} function works 
pielous <- microbiome::evenness(rarefied[[1]], 'pielou') %>% 
  as_tibble(rownames = "sample_id") %>% 
  rename(Pielou = pielou)

# for loop to iterate microbiome::evenness
pielous <- list("tbl_df") # "data.frame")  # 1. output
for (i in 1:length(rarefied) ) {            # 2. sequence
  pielous[[i]] <- microbiome::evenness(rarefied[[i]], 'pielou');  
  pielous[[i]] <- as_tibble(pielous[[i]], rownames = "sample_id") # 3. body
}
pielous

pielous[[1]]
pielous[[1]] %>% class()

## Bind all data.frames together to make single (long) table -----------------------
pielous_rbind <- bind_rows(pielous) #, .id = "column_label") this if I want to keep the [[i]] info as column

# Sanity check: n must be = as nr of iteration
pielous_rbind %>% 
  group_by(sample_id) %>% 
  summarise(n = length(pielou))


## Calculate mean (summarize) -----------------
pielous_rbind <- pielous_rbind %>% 
  group_by(sample_id) %>% 
  summarise_if(is.numeric, mean) #%>% view()



# Faith phylogenitic distance - requires {btools} -----------------------------------

# install.packages("remotes")
# remotes::install_github("twbattaglia/btools")

# Make sure {phyloseq} is loaded! Else will through error!!!

# Check that {btools} function works 
btools::estimate_pd(rarefied[[1]]) 

# for loop to iterate microbiome::evenness
faiths <- list("tbl_df") # "data.frame")  # 1. output
for (i in 1:length(rarefied) ) {            # 2. sequence
  faiths[[i]] <- btools::estimate_pd(rarefied[[i]]);  
  faiths[[i]] <- as_tibble(faiths[[i]], rownames = "sample_id") # 3. body
}
faiths

## Bind all data.frames together to make single (long) table -----------------------
faiths_rbind <- bind_rows(faiths) #, .id = "column_label") this if I want to keep the [[i]] info as column

# Sanity check: n must be = as nr of iteration
faiths_rbind %>% 
  group_by(sample_id) %>% 
  summarise(n = length(PD))


## Calculate mean (summarize) -----------------
faiths_rbind <- faiths_rbind %>% 
  group_by(sample_id) %>% 
  summarise_if(is.numeric, mean) #%>% view()





# Combine all alpha div tbl + metadata_ms for correct samples names ---------

metadata_ms <- read_csv("./in/metadata_ms.csv")

# phyalphas_rbind, pielous_rbind, faiths_rbind

# correct sample_id here that for some reason got - replaced by . ...
phyalphas_rbind <- phyalphas_rbind %>% 
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-")) #%>% view()

alphas <- left_join(phyalphas_rbind, pielous_rbind, by = "sample_id") %>% 
  left_join(., faiths_rbind, by = "sample_id") %>% 
  rename(Faith_PD = PD) %>% 
  select(-SR) # same as observed 


## Add metadata ---------------------
alphas <- left_join(alphas, metadata_ms, by = "sample_id") %>% 
  relocate(13:21, .after = sample_id) # %>% view()



# Export as RDS ---------------
saveRDS(alphas, "./out/RDS_files/rarefied1000_alpha_div_metrics.rds")


# xxxxxxxxxxxx ----------

alphas <- readRDS("./out/RDS_files/rarefied1000_alpha_div_metrics.rds")


# All ALPHA's - by group ---------------------

# Make a 3 summary, after grouping by:
# 1. state
# 2. state, origin
# 3. state, origin, colony

# 1. state
alphas %>% 
  group_by(state) %>% mutate(
    n = length(sample_id) ) %>% 
  summarise_if(is.numeric, mean) %>% 
  relocate(n, .after = state) %>%
  write_csv(., "./out/Gfas_16S/alpha_diversity/rar1000_alphas_summarybystate.csv")

# 2. state, origin
alphas %>% 
  group_by(state, origin_ms) %>% 
  summarise_if(is.numeric, mean) %>% 
  ungroup()


# 3. state, origin, colony
alphas %>% 
  group_by(state, origin_ms, colony_ms) %>% 
  mutate(
    n = length(sample_id) ) %>% 
  summarise_if(is.numeric, mean) %>% 
  relocate(n, .after = colony_ms) %>%
  select(-levels_colony_ms) %>% 
  ungroup() %>% 
  write_csv(., "./out/Gfas_16S/alpha_diversity/rar1000_alphas_summarybystateorigincolony.csv")


# PLOT ALPHA diversities by STATE (x axis) and ORIGIN --------------------------------------


# Function to automate plotting to explore all alpha diversity indices
alpha_plots <- function(data, y) {
  y <- enquo(y)
  
  ggplot(data = data, aes(x = fct_rev(state), 
                          y = !!y)) +
    geom_boxplot(aes(group = state),
                 outlier.shape = NA,
                 width = 0.4,
                 color = "gray33") + 
    ggbeeswarm::geom_quasirandom(
      aes(color = colony_ms, shape = state, size = state), # group = origin_ms, 
      width = 0.2,
      stroke = 1.3  ) +
    # geom_point(
    #   position = position_jitterdodge(
    #     jitter.width = 0.3, dodge.width = 0.6, seed = 23),
    #   stroke = 1.3,
    #   aes(color = colony_ms, group = origin_ms, shape = state, size = state)) + 
    scale_color_manual(values = palette_GfasMS) +
    scale_shape_manual(values = c(16, 1), breaks = c("symbiotic", "bleached")) +
    scale_size_manual(values = c(2, 1.5), breaks = c("symbiotic", "bleached")) +
    theme_classic() +
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

alpha_plots(data = alphas, y = Faith_PD)

obs <- alpha_plots(data = alphas, y = Observed)
chao <- alpha_plots(data = alphas, y = Chao1) 
piel <- alpha_plots(data = alphas, y = pielou) 
shan <- alpha_plots(data = alphas, y = Shannon)
simp <- alpha_plots(data = alphas, y = Simpson)  
# invsimp <- alpha_plots(data = alphas_cleaned, y = InvSimpson)
# fish <- alpha_plots(data = alphas_cleaned, y = Fisher)
fait <- alpha_plots(data = alphas, y = Faith_PD)
  
library(patchwork)
patch <- (obs + chao + piel) / (shan + simp + fait)

out <- patch + plot_layout(guides = "collect")

ggsave("./out/Gfas_16S/alpha_diversity/rarefied1000_alpha_div_metrics_by_state.png", 
       bg = "white",
       dpi = 330,
       units = "cm", width = 23, height = 20)


# After plotting each alpha diversity it does not look like 
# there's any significant difference between symb and bleached




# PLOT C.I. 95 % of ALPHA diversities by STATE (x axis) and ORIGIN - for Supplementary --------------------------------------

# Function to automate plotting to explore all alpha diversity indices
alpha_CI95_plots <- function(data, y) {
  y <- enquo(y)
  
  ggplot(data = data, aes(x = fct_rev(state), 
                          y = !!y )) + # log(!!y) #  sqrt(!!y)
        stat_summary(aes(
          color = fct_reorder(colony_ms, levels_colony_ms),
          size = fct_reorder(colony_ms, levels_colony_ms)),
                 # fun.data = mean_cl_normal, # requires {Hmisc}!
                 fun.data = mean_cl_boot, # requires {Hmisc}!
                 fun.args = list(conf.int = 0.95),
                 geom = "errorbar",
                 position = position_dodge(width = 0.5),
                 width = 0.5
                 # size = 1
    ) +
    scale_color_manual(values = palette_GfasMS) +
    scale_size_manual(values = c(1.5, rep(0.7, 4)),
                      guide = "none") +
    labs(color = "Colony") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.ticks.x = element_blank(),
      text = element_text(size = 13),
      panel.grid = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.spacing = unit(0.05, 'cm'),
      legend.key.size = unit(0.5, 'cm'),
      # legend.justification = "top",
      legend.position = "bottom",
      # aspect.ratio = 3,
      plot.margin = margin(t = 0, r = 0, l = 0.7, b = 0, "cm")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = 2)))
}

alpha_CI95_plots(data = alphas, y = Faith_PD)

obs <- alpha_CI95_plots(data = alphas, y = Observed)
chao <- alpha_CI95_plots(data = alphas, y = Chao1) 
piel <- alpha_CI95_plots(data = alphas, y = pielou) 
shan <- alpha_CI95_plots(data = alphas, y = Shannon)
simp <- alpha_CI95_plots(data = alphas, y = Simpson)  
# invsimp <- alpha_CI95_plots(data = alphas, y = InvSimpson)
# fish <- alpha_CI95_plots(data = alphas, y = Fisher)
fait <- alpha_CI95_plots(data = alphas, y = Faith_PD)

library(patchwork)
patch <- (obs + shan + simp + fait + piel) 


  out <- patch + plot_layout(guides = "collect") + 
  plot_layout(guides = "collect", ncol = 5 ) &
  theme(legend.position = "bottom")


ggsave("./out/Gfas_16S/alpha_diversity/rarefied1000_alphas_bystateCI95.png", 
       bg = "white",
       dpi = 330,
       units = "cm", width = 30, height = 13)


ggsave("./out/Gfas_16S/alpha_diversity/rarefied1000_alphas_bystateCI95.svg", 
       bg = "white",
       # dpi = 330,
       units = "cm", width = 30, height = 13)



# Plot for MS: Obs, Shannon, Simpson ----------------------------------------

ms <- (obs + shan + simp)

ms_out <- ms + 
  plot_annotation(
    theme = theme(
      legend.justification = "top")) +
  plot_layout(guides = "collect")


ggsave("./out/Gfas_16S/alpha_diversity/rarefied1000_alpha_MS.png", 
       bg = "white",
       dpi = 330,
       units = "cm", width = 19, height = 13)





# Stat tests -----------------------
# Test difference is ALPHA diversity indexes between Symb and Bleached

# library("ggpubr")

# Visual inspection of normality of distribution
# All good except two ... (Shannon and Simpson)
obs <- ggpubr::ggqqplot(alphas, x = "Observed", title = "Observed")
chao <- ggpubr::ggqqplot(alphas, x = "Chao1", title = "Chao1")
piel <- ggpubr::ggqqplot(alphas, x = "pielou", title = "Pielou")   # very bad
shan <- ggpubr::ggqqplot(alphas, x = "Shannon", title = "Shannon")  # very bad
simp <- ggpubr::ggqqplot(alphas, x = "Simpson", title = "Simpson")  # worst
fait <- ggpubr::ggqqplot(alphas, x = "Faith_PD", title = "Faith_PD")


patch <- (obs + chao + piel) / (shan + simp + fait)

out <- patch + plot_layout(guides = "collect") +
  plot_annotation(
    title = 'Inspect normality of distribution',
    theme = theme(
      plot.title = element_text(size = 20, face = "bold")) )

ggsave("./out/Gfas_16S/alpha_diversity/rarefied1000_alpha_div_metrics_ggqqplot.png", 
       bg = "white",
       dpi = 330,
       units = "cm", width = 27, height = 20) 



# Homogeneity of variance - F-test ......................................

var.test(data = alphas, Observed ~ state) # 0.6217
var.test(data = alphas, Chao1 ~ state)    # 0.6085
var.test(data = alphas, pielou ~ state)   # 0.1134
var.test(data = alphas, Shannon ~ state)  # 0.1437
var.test(data = alphas, Simpson ~ state)  # 0.04222 *
var.test(data = alphas, Faith_PD ~ state) # 0.3703


# Homogeneity of variance - BARTLETT test ......................................

# H0 = the variances in each of the groups (samples) are the same
# -> want p.value > 0.05 -> can proceed with t-test

# values ~ groups
bartlett.test(data = alphas, Observed ~ state) # 0.627
bartlett.test(data = alphas, Chao1 ~ state)    # 0.6138
bartlett.test(data = alphas, pielou ~ state)   # 0.1151
bartlett.test(data = alphas, Shannon ~ state)  # 0.1457
bartlett.test(data = alphas, Simpson ~ state)  # 0.04304 *
bartlett.test(data = alphas, Faith_PD ~ state) # 0.3742

# Same results as with the F-test (above)
# can do regular t-test for all but NOT SIMPSON (will need Welch test)


# t-test ----------------------------------------------------------------------

# t-test is used when the population has n < 30

# "greater" is the alternative that x has a larger mean than y
#  (bleached is greater that symbiotic?)

# library(rstatix)
(tt_obs <- rstatix::t_test(data = alphas, Observed ~ state,
                alternative = "greater", # "two.sided"
                paired = F,
                var.equal = T)  )   # two.sided 0.169 # less 0.916 # greater 0.0843

(tt_chao <- rstatix::t_test(data = alphas, Chao1 ~ state,
                alternative = "greater", # "two.sided", #  #
                paired = F,
                var.equal = T)   )                     # 0.0874 
(tt_piel <- rstatix::t_test(data = alphas, pielou ~ state,
                alternative = "greater", # "two.sided"
                paired = F,
                var.equal = T)   )                     # 0.0804
(tt_shan <- rstatix::t_test(data = alphas, Shannon ~ state,
                alternative = "greater", # "two.sided"
                paired = F,
                var.equal = T)    )                    # 0.0631
(tt_simp <- rstatix::t_test(data = alphas, Simpson ~ state,
                alternative = "greater", # "two.sided"              
                paired = F,
                var.equal = FALSE)   )                 # 0.0816  unequal variance here!***
(tt_fpd <- rstatix::t_test(data = alphas, Faith_PD ~ state,
                alternative = "greater", # "two.sided"                 
                paired = F,
                var.equal = T)   )                     # 0.0802  # 0.16

          
ttests <- rbind(tt_obs, tt_shan, tt_simp, tt_fpd, tt_piel, tt_chao) %>% 
  rename(div_metric = ".y.") %>% 
  mutate(
    test = ifelse(div_metric == "Simpson", "Welch t-test", "unpaired 2-samples t-test"),
    variance = ifelse(div_metric == "Simpson", "unequal", "equal")) %>% 
  relocate(test, .after = div_metric) %>% 
  relocate(variance, .after = test)

write_csv(ttests, "./out/Gfas_16S/alpha_diversity/rar1000_alphas_ttests.csv")

# None of these is significant
# Neither when I do "two.sided" 
# nor when one-sided (group 1 is "greater" than group 2),
# which should be more sensitive, and would make sense to ask
# if symbiotic is smaller (alpha diversity) that bleached
# (based on the Anna Karenina Pr. (AKP) of microbial ecology)


# Wilcoxon test (Mann-Whitney) ----------------------------------------------------------------------

# Non-parametric 
# necessary for Pielou, Shannon, and Simpson ...

# (Note that before repeated rarefaction (rar1000) 
#  I was getting warning: ties occur because I have some identical values
#  bot not anymore :) )
alphas$pielou %>% duplicated() %>% which(. == T) %>% length()  # before repeated rar: 7
alphas$Shannon %>% duplicated() %>% which(. == T) %>% length() 
alphas$Simpson %>% duplicated() %>% which(. == T) %>% length() 


# Pielou
stats::wilcox.test(data = alphas, 
                   pielou ~ state,
                   alternative = "greater" # W = 117, p-value = 0.1099
                   # alternative = "two.sided"
                   ) 

(wcx_piel <- rstatix::wilcox_test(data = alphas, 
                     pielou ~ state,
                     alternative = "greater" #  W = 117, p = 0.11
                     # alternative = "two.sided" #  W = 117, p = 0.22
                     # alternative = "less"
                     )
)

# Shannon
(wcx_shan <- rstatix::wilcox_test(data = alphas, 
                   Shannon ~ state,
                   alternative = "greater" # W = 113, p-value = 0.151
                   # alternative = "two.sided"
                   ) 
)

# Simpson
(wcx_simp <- rstatix::wilcox_test(data = alphas, 
                   Simpson ~ state,
                   alternative = "greater" # W = 123, p-value = 0.0639
                   # alternative = "two.sided"
                   )  
)


wcxtests <- rbind(wcx_shan, wcx_simp, wcx_piel) %>% 
  rename(div_metric = ".y.") %>% 
  mutate(
    test = " Mann-Whitney U-test" #  Wilcox is for dependent samples
    ) %>%
  relocate(test, .after = div_metric) 

write_csv(wcxtests, "./out/Gfas_16S/alpha_diversity/rar1000_alphas_wcxtests.csv")


# Again, NOT significant difference between symbiotic states
# for all alpha diversity metrics considered 

# alternative hypothesis: true location shift is not equal to 0
# is not rejected (p > 0.05)

# Note that p.value is smaller when 'alternative = "greater"' but 
# is still always above 0.05
