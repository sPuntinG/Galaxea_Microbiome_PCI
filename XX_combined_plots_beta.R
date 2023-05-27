#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multi-panel figure - BETA DIVERSITY (nMDS and pw Bray-Curtis boxplot)
# 
# Combine results from different scripts into single
# plots to have consistent layouts (for publication).
#
# Used rarefied x 1000 data!
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(tidyverse)
library(here)

# Palette
source("colors_GfasMS.R")


# BETA diversity (nMDS's + boxplot of B-C pwd) ------------------------

## Import .rds ----------------

# nMDS
all <- readRDS("./out/RDS_files/nMDS_bray1000.rds") 
ordn_bray_all <- readRDS("./out/RDS_files/ordn1000_bray_all.rds") 

RS <- readRDS("./out/RDS_files/nMDS_bray1000_RS.rds") 
ordn_bray_RS <- readRDS("./out/RDS_files/ordn1000_bray_RS.rds") 
HK <- readRDS("./out/RDS_files/nMDS_bray1000_HK.rds") 
ordn_bray_HK <- readRDS("./out/RDS_files/ordn1000_bray_HK.rds")

# pw Bray-Curtis boxplot
pwcBC_n3 <- readRDS("./out/RDS_files/pwcBC1000_n3_tibble.rds") 
w_n3 <- readRDS("./out/RDS_files/rar1000_w_n3.rds") 



## Ordination plots: general settings for plots --------------

theme_set(theme_bw())







### Figure 2 A - NEW (colony by shape) ----------------------------------------------

# Text: test results

# Plot 
(nMDS_all_shapes <- ggplot() +
   stat_ellipse(data = all,
                aes(x = MDS1, y = MDS2,
                    color = origin_ms),
                level = 0.95,
                size = 0.5,
                show.legend = F) + 
   scale_color_manual(breaks = c("Red Sea", "Hong Kong"),
                      values = c(Red_Sea, Hong_Kong),
                      guide = "none") +
   ggnewscale::new_scale_color() +
   geom_point(data = all,
              inherit.aes = F,
              aes(x = MDS1, y = MDS2,
                  color = colony_ms, 
                  shape = colony_ms,
                  fill = state
                  ), 
              stroke = 1)  +
   scale_shape_manual(
     breaks = c("RS1", "RS2", "RS3", "HK1", "HK2"),
     values = c(21:25)
   ) +
   scale_fill_manual(
     breaks = c("symbiotic", "bleached"),
     values = c("black", "white"), # then manually replace in Inkscape ...
     labels =  c("Symbiotic", "Bleached")) +
   scale_color_manual(values = palette_GfasMS,
                      guide = guide_legend(order = 2)) + 
   labs(
     subtitle = paste0("Stress = ", round(ordn_bray_all$stress, 3))
   ) +
   scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.5)) +
   scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.5)) +
   coord_fixed(ratio = 1, # necessary for correct ratio of axes!
               xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
   theme(
     panel.grid = element_blank(),
     plot.subtitle = element_text(hjust = 0.02, vjust = -10, 
                                  size = 6,
                                  color = "black"), # top left
     text = element_text(size = 13),
     axis.text = element_text(size = 10, colour = "black"),
     axis.title = element_text(size = 10),
     legend.title = element_blank() ,
     legend.position = "none"
   )
)


ggsave("./out/Gfas_16S/beta_diversity/rar1000_nMDS_originellipses_all_shapes.png",
       # bg = "white",
       dpi = 330,
       units = "cm", width = 15.5, height = 12.5)





### Figure 2 B - NEW (pw B-C dist, boxplot) ----------------------------------------------

# Text: test results

wlcx_md <- paste0("Mann-Whitney U test<br>",
                  "*W* = ", w_n3$statistic, ", *P* = **", # w_n3$method
                  format(round(w_n3$p.value, 4), scientific = F),
                  "**")

# Plot
(
  pwcBDn3_shapes <- ggplot(pwcBC_n3, aes(
    x = forcats::fct_rev(state), 
    y = pairwise_distance,
    color = colony_ms,
    shape = colony_ms,
    fill = state)) +
    geom_boxplot(inherit.aes = F,
                 aes(x = forcats::fct_rev(state),
                     y = pairwise_distance),
                 show.legend = F, 
                 color = "grey69", 
                 outlier.shape = NA,
                 ) +
    ggbeeswarm::geom_quasirandom(aes(fill = state),
                                 width = 0.2,
                                 stroke = 1) + # 1.5
    scale_color_manual(values = palette_GfasMS) +
    scale_fill_manual(breaks = c("S", "B"),
                      values = c("black", "white"),
                      labels = c("Symbiotic", "Bleached")) +
    scale_shape_manual(
      breaks = c("RS1", "RS2", "RS3", "HK1", "HK2"),
      values = c(21:25)) +
    scale_x_discrete(breaks = c("S", "B"), 
                     labels = c("Symbiotic", "Bleached")) +
    labs(
      y = "Bray-Curtis pairwise distance",
      color = "Colony",
      shape = "Colony",
      fill = "State"
    ) +
    ggtext::geom_richtext(aes(x = 0.5, y = 0.9,
                              label = wlcx_md),
                          stat = "unique", 
                          color = "black", fill = NA,
                          label.color = NA,
                          size = 2,
                          hjust = 0 ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(color = "black", size = 10),
      axis.title = element_text(size = 10),
      panel.grid = element_blank(),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.position = "right",
      legend.justification = "top",
      aspect.ratio = 2
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = c(21, 21))  ))

    )


ggsave("./out/Gfas_16S/beta_diversity/rar1000_pwcBDn3_shapes.png",
       # bg = "white",
       dpi = 330,
       units = "cm", width = 15.5, height = 12.5)



## Panel figures -----------------

### TOP: nMDS_all + pwc_BC_n3 ------------------------

#### .{Patchwork} - good -------------------------------
library(patchwork)

(
  AB <- nMDS_all_shapes + pwcBDn3_shapes 
)


ggsave("./out/Gfas_16S/beta_diversity/rar1000_nMDSall_pwBCn3_legendR_shapes.png",
       # bg = "white",
       dpi = 330,
       units = "cm", width = 18, height = 10) # Keep dimensions to match with CD panels (patchwork below)


ggsave("./out/Gfas_16S/beta_diversity/rar1000_nMDSall_pwBCn3_legendR_shapes.svg",
       # bg = "white",
       dpi = 330,
       units = "cm", width = 18, height = 10) # Keep dimensions to match with CD panels (patchwork below)





##### Shape by colony & Legend outside ----------------------

# facet combine nMDS data
  
RSHK <- bind_rows(RS, HK)


(nMDS_RSHK_out <- ggplot(data = RSHK, 
                         aes(x = MDS1, 
                             y = MDS2,
                             # shape = state,
                             fill = state,
                             shape = colony_ms,
                             group = state, 
                             linetype = state,
                             color = origin_ms)) +
  geom_point() +
  scale_shape_manual(
    breaks = c("RS1", "RS2", "RS3", "HK1", "HK2"),
    values = c(21:25)
    # labels =  c("Symbiotic", "Bleached")
    ) +
  scale_fill_manual(
      breaks = c("symbiotic", "bleached"),
      values = c(Red_Sea, "white"),
      labels =  c("Symbiotic", "Bleached")) +
  stat_ellipse(level = 0.95,
               size = 0.5) +
  scale_linetype_manual(
    breaks = c("symbiotic", "bleached"), 
    values = c("solid", "dotted"),
    labels =  c("Symbiotic", "Bleached")
  ) +
  scale_color_manual(values = c(Red_Sea, Hong_Kong),
                     breaks = c( "Red Sea", "Hong Kong")) +
  coord_fixed(ratio = 1, # necessary for correct ratio of axes!
              xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.background = element_rect(fill = NA),
    legend.spacing = unit(0.1, 'cm'),
    legend.key.size = unit(0.5, 'cm'),
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~fct_rev(origin_ms)) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        linetype = c("solid", "dotted")) )
    )
)

ggsave("./out/Gfas_16S/beta_diversity/rar1000_Bray_nMDS_ellipses_colonyshape_RSHK.png",
       bg = "white",
       dpi = 330,
       units = "cm", width = 17, height = 15)

ggsave("./out/Gfas_16S/beta_diversity/rar1000_Bray_nMDS_ellipses_colonyshape_RSHK.svg",
       bg = "white",
       dpi = 330,
       units = "cm", width = 17, height = 15)





### For SUPPLEMENTARY ------------------------------

#### nMDS with ellipses by STATE -----------------

# show overlapping ellipses = won't find a significant difference with PERMDISP
# nor PERMANOVA (see "PERMANOVA with vegan::ADONIS2" lines 386-435)

ggplot(data = all,  # nMDS_bray
       aes(x = MDS1, 
           y = MDS2,
           shape = state,
           group = state, linetype = state)) +
  geom_point() +
  scale_shape_manual(
    breaks = c("symbiotic", "bleached"), 
    values = c(16, 1),
    labels =  c("Symbiotic", "Bleached")) +
  stat_ellipse(level = 0.95,
               size = 0.5) +
  scale_linetype_manual(
    breaks = c("symbiotic", "bleached"),
    values = c("solid", "dotted"),
    labels =  c("Symbiotic", "Bleached")) +
  scale_x_continuous(
    breaks = seq(-0.5, 0.5, by = 0.5)) +
  scale_y_continuous(
    breaks = seq(-0.5, 0.5, by = 0.5)) +
  coord_fixed(ratio = 1, # necessary for correct ratio of axes!
              xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.justification = "top"
  )


ggsave("./out/Gfas_16S/beta_diversity/rar1000_Bray_nMDS_ellipses_state.png", 
       bg = "white",
       dpi = 330, 
       units = "cm", width = 15.5, height = 12.5)


