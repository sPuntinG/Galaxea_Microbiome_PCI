#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINED PLOTS
# Combine results from different scripts into single
# plots to have consistent layouts (for publication).
#
# This an the other "XX_*.R" are largely redundant with previous  
#  scripts wrt the code for plotting. But I share them anyways in 
#  in case other fellas would be interested in ggplot tricks
#  for multi-panel plots :)
#
# Used rarefied x 1000 data!
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(tidyverse)
library(here)

# Palette
source("colors_GfasMS.R")



# ALPHA diversity + accumulation curve -------------------------

# Import data 

# Alpha diversity metrics
alphas <- readRDS("./out/RDS_files/rarefied1000_alpha_div_metrics.rds ") 
                                # OLD: rarefied_alpha_div_metrics.rds


# Accumulation curve
acc_curve <- readRDS("./out/RDS_files/ASV_cumsum_bind.rds")
# acc_curve_ggplot <- readRDS("./out/RDS_files/ggplot_acccurve_5gr.rds")


## Alpha - statistical tests -----------------------

# Test the difference between symbiotic and bleached, for each alpha div index
# Based on checks of assumptions (see '10_alpha_diversity.R') I can use
#  - t-test for Observed richness and Faith_pd
#  - wilcox.test for Simpson and Shannon

# Here re-do tests so I have the data to insert in the plots 

# t-test Observed
tt_obs <- rstatix::t_test(data = alphas, Observed ~ state,
                alternative = "greater", # "two.sided"
                paired = F,
                var.equal = T)                        

tt_obs$statistic
tt_obs$df
tt_obs$p


# t-test Faith pd
tt_fpd <- rstatix::t_test(data = alphas, Faith_PD ~ state,
                alternative = "greater", # "two.sided"                 
                paired = F,
                var.equal = T)                        
tt_fpd$statistic
tt_fpd$df
tt_fpd$p


# wilcox.test Simpson
wr_sim <- stats::wilcox.test(data = alphas, 
                   Shannon ~ state,
                   alternative = "greater") # W = 66.5, p-value = 0.244
wr_sim$method
wr_sim$statistic
wr_sim$p.value


# wilcox.test Shannon
wr_sha <- stats::wilcox.test(data = alphas, 
                   Simpson ~ state,
                   alternative = "greater") # W = 66, p-value = 0.2318 
wr_sha$method
wr_sha$statistic
wr_sha$p.value



## Plots: alpha div -------------------
alpha_plots <- function(data, y, ylab, test, Pval) {
  y <- enquo(y)
  
  ggplot(data = data, 
         aes(x = fct_rev(state), 
             y = !!y)) +
    geom_boxplot(aes(group = state),
                 outlier.shape = NA,
                 width = 0.4,
                 color = "gray33") + 
    ggbeeswarm::geom_quasirandom(
      aes(color = colony_ms, shape = state, size = state), # group = origin_ms, 
      width = 0.2,
      stroke = 1.3  ) +
    scale_color_manual(values = palette_GfasMS) +
    scale_shape_manual(values = c(16, 1), 
                       breaks = c("symbiotic", "bleached"),
                       labels = c("Symbiotic", "Bleached")) +
    scale_size_manual(values = c(2, 1.5), 
                      breaks = c("symbiotic", "bleached"),
                      labels = c("Symbiotic", "Bleached"),
                      guide = "none") +
    scale_x_discrete(breaks = c("symbiotic", "bleached"), 
                     labels = c("Symbiotic", "Bleached")) +
    labs(
      y = ylab,
      subtitle = paste0(test, "<br>*P* = ", round(Pval, 3)),
      shape = "State",
      color = "Colony"
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      # axis.title.y = element_text(size = 15, face = "bold"),
      # text = element_text(size = 13),
      # axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      # legend.box = "vertical", # works here but not after patchwork
      legend.text = element_text(size = 7, vjust = 0.6),
      legend.title = element_text(size = 7.5, vjust = 0.6, face = "bold"),
      # legend.title = element_blank(),
      legend.spacing = unit(0.05, 'cm'),
      legend.key.size = unit(0.3, 'cm'),
      aspect.ratio = 3,
      plot.margin = margin(t = 0, r = 0, l = 0.7, b = 0, "cm"),
      plot.subtitle = ggtext::element_markdown(size = 6, vjust = 0, hjust = 0)
      ) +
    guides(
      shape = guide_legend(
        override.aes = list(
          size = c(1.9, 1.3))  ) # 
      # color = guide_legend(
      #   override.aes = list(
      #     size = 1.2)  )
    )
}

alpha_plots(data = alphas, y = Faith_PD, ylab = "Faith PD", test = "t-test", Pval = tt_fpd$p)

obs <- alpha_plots(data = alphas, y = Observed, ylab = "Obs. nr of ASVs", test = "t-test", Pval = tt_obs$p)
sha <- alpha_plots(data = alphas, y = Shannon, ylab = "Shannon diversity", test = "Mann-Whitney U-test", Pval = wr_sha$p.value)
sim <- alpha_plots(data = alphas, y = Simpson, ylab = "Simpson evenness", test = "Mann-Whitney U-test", Pval = wr_sim$p.value)  
fpd <- alpha_plots(data = alphas, y = Faith_pd, ylab = "Faith PD", test = "t-test", Pval = tt_fpd$p)



### Patch alphas ------------------------

library(patchwork)

patch <- (obs + sha + sim)

(alpha3 <- patch + plot_layout(guides = "collect") &
    theme(legend.position = "bottom"
          # legend.box = 'vertical'
          )
  )





## Plot: Accumulation curve plot --------------

(
  acccurve <- ggplot(data = acc_curve, #ASV_cumsum_bind,
                                
                                aes(x = rank, y = rel_ab_cumsum,
                                    color = origin_state,
                                    fill = origin_state)) +
    geom_line() +
    geom_point(shape = 21,
               size = 1) + 
    scale_color_manual(breaks = c("Red Sea symbiotic", "Red Sea bleached", 
                                  "Hong Kong symbiotic", "Hong Kong bleached",
                                  "Overall"),
                       values = c(Red_Sea, Red_Sea,
                                  Hong_Kong, Hong_Kong,
                                  "black")) +
    scale_fill_manual(breaks = c("Red Sea symbiotic", "Red Sea bleached", 
                                 "Hong Kong symbiotic", "Hong Kong bleached",
                                 "Overall"),
                      values = c(Red_Sea, "white", # "transparent", #"white",
                                 Hong_Kong, "white", #"transparent", # 
                                 "black")) +
    scale_y_continuous(limits = c(0, 100), 
                       labels = scales::percent_format(scale = 1),
                       # suffix = " %", # Added recently - try 
                       expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 75), 
                       breaks = seq(0, 150, 10), 
                       expand = c(0, 0)) +
    xlab("nr of ASVs") +
    ylab("Cumulative abundance") +
    theme_classic() +
    theme(
      axis.ticks = element_blank(),
      panel.grid.major.y = element_line(colour = "grey96"), # "grey" linetype = "dotted"),
      panel.grid.major.x = element_line(colour = "grey96"), # "grey75" linetype = "dotted"),
      # legend.position = "none"
      legend.text = element_text(size = 7, vjust = 0.6),
      legend.spacing = unit(0.05, 'cm'),
      legend.key.size = unit(0.3, 'cm'),
      legend.title = element_blank(),
      legend.position = "bottom",
      # legend.justification = "centre", # doesn't do anything
      aspect.ratio = 1
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = F)
    )
)





## Patchwork --------------------

layout <- "
ABCD
"

(
  patch <- (alpha3 + acccurve) + 
  plot_layout(
    # guides = "collect",
    design = layout,
    nrow = 1
  ) +
  plot_annotation(tag_levels = 'A') &
    theme(
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.tag = element_text(size = 10),
      # plot.background = element_rect(fill = "yellow"),
      plot.margin = margin(t = 0, r = 0.1, l = 0, b = 0, "cm")
    )
)


## Export .PNG -----------------------------------
ggsave("./out/Gfas_16S/alpha_diversity/rar1000_alphas3_acccurve.png",
       bg = "white",
       dpi = 330,
       units = "cm", width = 21, height = 11)


## Export .SVG -----------------------------------
library(svglite)

ggsave("./out/Gfas_16S/alpha_diversity/rar1000_alphas3_acccurve.svg",
       bg = "white",
       # dpi = 330, 
       units = "cm", 
       width = 21, height = 11)

