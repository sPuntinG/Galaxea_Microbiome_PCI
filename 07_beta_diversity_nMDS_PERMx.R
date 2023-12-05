
library(tidyverse)
library(here)
library(phyloseq)

#*** source("06_phyloseq_proj_and_rarefaction.R")
metadata_ms <- read_csv("./in/metadata_ms.csv")

# 1. Rarefied - single subsampling ========================================
ps_rarefied <- readRDS("./out/RDS_files/ps_rarefied.rds")


# Transformed vs non-transformed, weightedUniFrac vs Bray-C data


## Transform data -------------------------------------------------

# Note: phyloseq::transform_sample_counts() works sample-by-sample

# To relative abundance (https://joey711.github.io/phyloseq/preprocess.html)
ps_rar_rel = phyloseq::transform_sample_counts(ps_rarefied, function(x) x / sum(x) )

# To square root
ps_rar_sqrt = phyloseq::transform_sample_counts(ps_rarefied, function(x) sqrt(x) )

# To log(x + 1)
# log() = natural log
ps_rar_logx1 = phyloseq::transform_sample_counts(ps_rarefied, function(x) log(x + 1) )


# # To visually/manually check ...
# phyloseq::otu_table(ps_rarefied) %>% view()
# phyloseq::otu_table(ps_rar_rel) %>% view()
# phyloseq::otu_table(ps_rar_sqrt) %>% view()
# phyloseq::otu_table(ps_rar_logx1) %>% view()



## Compute BETA diversity distance matrix with {phyloseq} ----------------------------------------

# Make sample-wise distance object (class: "dist")
# summary(dist) # "Distance matrix by lower triangle"

### Weighted UniFrac ----------------------------------------

# Non-transformed data
dist_wuf = phyloseq::distance(ps_rarefied, method = "wunifrac")

# Rel abundance
dist_wuf_relab = phyloseq::distance(ps_rar_rel, method = "wunifrac")

# Sqrt-transformed data
dist_wuf_sqrt = phyloseq::distance(ps_rar_sqrt, method = "wunifrac")

# log(x + 1)-transformed data
dist_wuf_logx1 = phyloseq::distance(ps_rar_logx1, method = "wunifrac")


### Bray-Curtis ----------------------------------------------------

# Non-transformed data
dist_bray = phyloseq::distance(ps_rarefied, method = "bray")

# Rel abundance
dist_bray_relab = phyloseq::distance(ps_rar_rel, method = "bray")

# Sqrt-transformed data
dist_bray_sqrt = phyloseq::distance(ps_rar_sqrt, method = "bray")

# log(x + 1)-transformed data
dist_bray_logx1 = phyloseq::distance(ps_rar_logx1, method = "bray")


## Export RDS (for later scripts) ------------------
saveRDS(dist_bray, "./out/RDS_files/dist_bray.rds")



## Ordination plots (nMDS) - make and compare ----------------------------------------------------------------------

# Create "metaMDS" "monoMDS" object(s)

# Rule of thumb to evaluate nMDS plots based on "stress":
# < 0.05 excellent representation in reduced dimensions
# < 0.1 great,
# < 0.2 good/ok
# < 0.3 poor representation


### Create MDS objects -------------------------------------------------
#### Weighted UniFrac -------------------------------------------------

# Non-transformed data
ordn_wuf = phyloseq::ordinate(ps_rarefied, method = "NMDS", distance = dist_wuf)

# Rel abundance
ordn_wuf_relab = phyloseq::ordinate(ps_rar_rel, method = "NMDS", distance = dist_wuf_relab)

# Sqrt-transformed data
ordn_wuf_sqrt = phyloseq::ordinate(ps_rar_sqrt, method = "NMDS", distance = dist_wuf_sqrt)

# log(x + 1)-transformed data
ordn_wuf_logx1 = phyloseq::ordinate(ps_rar_logx1, method = "NMDS", distance = dist_wuf_logx1)


#### Bray-Curtis ------------------------------------------------

# Non-transformed data
ordn_bray = phyloseq::ordinate(ps_rarefied, method = "NMDS", distance = dist_bray)

# Rel abundance
ordn_bray_relab = phyloseq::ordinate(ps_rar_rel, method = "NMDS", distance = dist_bray_relab)

# Sqrt-transformed data
ordn_bray_sqrt = phyloseq::ordinate(ps_rar_sqrt, method = "NMDS", distance = dist_bray_sqrt)

# log(x + 1)-transformed data
ordn_bray_logx1 = phyloseq::ordinate(ps_rar_logx1, method = "NMDS", distance = dist_bray_logx1)



## ORDINATION PLOTS -------------------------------

#### General settings for plots ------------------------------

# Fonts
windowsFonts(pfont = windowsFont("Calibri")) 

# Set theme for the whole script
theme_set(theme_bw(base_family = "pfont"))

# Palette
source("colors_GfasMS.R")
palette_GfasMS


# ORDINATION PLOTS: function .................................

my_nMDS <- function(data, ordn, title, ...) {
  phyloseq::plot_ordination(data, ordn, color = "colony_ms", shape = "state") +
    scale_shape_manual(breaks = c("bleached", "symbiotic"), values = c(1, 16)) +
    scale_color_manual(values = palette_GfasMS) + # breaks = brks, 
    labs(
      title = title,
      subtitle = paste0("Stress = ", round(ordn$stress, 3)) ) +
    stat_ellipse(aes(group = origin, linetype = origin)) + # , linetype = 2
    theme(
      panel.grid = element_blank(),
      plot.subtitle = element_text(hjust = 1)
    )
}


#### Weighted Uni Frac ----------------------

# Non-transformed data
a <- my_nMDS(
  data = ps_rarefied,
  ordn = ordn_wuf,
  title = "Weighted UniFrac" )

# Relative abundance
aa <- my_nMDS(
  data = ps_rar_rel,
  ordn = ordn_wuf_relab,
  title = "Weighted UniFrac, relative abundance" )


# Sqrt-transformed data
b <- my_nMDS(
  data = ps_rar_sqrt,
  ordn = ordn_wuf_sqrt,
  title = "Weighted UniFrac, sqrt()"
)


# log(x + 1)-transformed data
cc <- my_nMDS(
  data = ps_rar_logx1,
  ordn = ordn_wuf_logx1,
  title = "Weighted UniFrac, log(x + 1)"
)



#### Bray-Curtis -----------------------------------------

# Non-transformed data
d <- my_nMDS(
  data = ps_rarefied,
  ordn = ordn_bray,
  title = "Bray-Curtis"
)

# Sqrt-transformed data
dd <- my_nMDS(
  data = ps_rar_rel,
  ordn = ordn_bray_relab,
  title = "Bray-Curtis, relative abundance"
)

# Sqrt-transformed data
ee <- my_nMDS(
  data = ps_rar_sqrt,
  ordn = ordn_bray_sqrt,
  title = "Bray-Curtis, sqrt()"
)


# log(x + 1)-transformed data
f <- my_nMDS(
  data = ps_rar_logx1,
  ordn = ordn_bray_logx1,
  title = "Bray-Curtis, log(x + 1)"
)


# Patchwork .............................................
library(patchwork)

pw <- ( a + aa + b + cc ) / ( d + dd + ee + f )


(
  pw +
  plot_layout(guides = "collect")
)

ggsave("./out/Gfas_16S/beta_diversity/rarefied_wUniF_BrayC_nMDS_8.png", 
       bg = "white",
       dpi = 330, 
       units = "cm", width = 19, height = 27) # width = 25, height = 17

# Now we can visually compare diversity and transformations:
# Untransformed Bray-Curtis seems the best
# -> use this one for the Galaxea coral model MANUSCRIPT

rm(list = setdiff(ls(), c("metadata_ms", 
                          "palette_GfasMS", "Red_Sea", "Hong_Kong")))





# 2. Rarefied - 1000x subsampling ========================================
ps_rarefied.1000avg <- readRDS("./out/RDS_files/ps_rarefied.1000avg.rds")


## Transform data -------------------------------------------------

# Note: phyloseq::transform_sample_counts() works sample-by-sample

# To relative abundance (https://joey711.github.io/phyloseq/preprocess.html)
ps_rar1000_rel = phyloseq::transform_sample_counts(ps_rarefied.1000avg, function(x) x / sum(x) )

# To square root
ps_rar1000_sqrt = phyloseq::transform_sample_counts(ps_rarefied.1000avg, function(x) sqrt(x) )

# To log(x + 1)
# log() = natural log
ps_rar1000_logx1 = phyloseq::transform_sample_counts(ps_rarefied.1000avg, function(x) log(x + 1) )


## Compute BETA diversity distance matrix with {phyloseq} ----------------------------------------

# Make sample-wise distance object (class: "dist")
# summary(dist) # "Distance matrix by lower triangle"

### Weighted UniFrac ----------------------------------------

# Non-transformed data
dist1000_wuf = phyloseq::distance(ps_rarefied.1000avg, method = "wunifrac")

# Rel abundance
dist1000_wuf_relab = phyloseq::distance(ps_rar1000_rel, method = "wunifrac")

# Sqrt-transformed data
dist1000_wuf_sqrt = phyloseq::distance(ps_rar1000_sqrt, method = "wunifrac")

# log(x + 1)-transformed data
dist1000_wuf_logx1 = phyloseq::distance(ps_rar1000_logx1, method = "wunifrac")


### Bray-Curtis ----------------------------------------------------

# Non-transformed data
dist1000_bray = phyloseq::distance(ps_rarefied.1000avg, method = "bray")

# Rel abundance
dist1000_bray_relab = phyloseq::distance(ps_rar1000_rel, method = "bray")

# Sqrt-transformed data
dist1000_bray_sqrt = phyloseq::distance(ps_rar1000_sqrt, method = "bray")

# log(x + 1)-transformed data
dist1000_bray_logx1 = phyloseq::distance(ps_rar1000_logx1, method = "bray")


## Export RDS (for later scripts) ------------------
saveRDS(dist1000_bray, "./out/RDS_files/dist1000_bray.rds")



## Ordination plots (nMDS) - make and compare ----------------------------------------------------------------------

# Create "metaMDS" "monoMDS" object(s)

# Rule of thumb to evaluate nMDS plots based on "stress":
# < 0.05 excellent representation in reduced dimensions
# < 0.1 great,
# < 0.2 good/ok
# < 0.3 poor representation


### Create MDS objects -------------------------------------------------
#### Weighted UniFrac -------------------------------------------------

# Non-transformed data
ordn1000_wuf = phyloseq::ordinate(ps_rarefied.1000avg, method = "NMDS", distance = dist1000_wuf)

# Rel abundance
ordn1000_wuf_relab = phyloseq::ordinate(ps_rar1000_rel, method = "NMDS", distance = dist1000_wuf_relab)

# Sqrt-transformed data
ordn1000_wuf_sqrt = phyloseq::ordinate(ps_rar1000_sqrt, method = "NMDS", distance = dist1000_wuf_sqrt)

# log(x + 1)-transformed data
ordn1000_wuf_logx1 = phyloseq::ordinate(ps_rar1000_logx1, method = "NMDS", distance = dist1000_wuf_logx1)


#### Bray-Curtis ------------------------------------------------

# Non-transformed data
ordn1000_bray = phyloseq::ordinate(ps_rarefied.1000avg, method = "NMDS", distance = dist1000_bray)

# Rel abundance
ordn1000_bray_relab = phyloseq::ordinate(ps_rar1000_rel, method = "NMDS", distance = dist1000_bray_relab)

# Sqrt-transformed data
ordn1000_bray_sqrt = phyloseq::ordinate(ps_rar1000_sqrt, method = "NMDS", distance = dist1000_bray_sqrt)

# log(x + 1)-transformed data
ordn1000_bray_logx1 = phyloseq::ordinate(ps_rar1000_logx1, method = "NMDS", distance = dist1000_bray_logx1)



## ORDINATION PLOTS -------------------------------

#### General settings for plots ------------------------------

# Fonts
windowsFonts(pfont = windowsFont("Calibri")) 

# Set theme for the whole script
theme_set(theme_bw(base_family = "pfont"))

# Palette
source("colors_GfasMS.R")
palette_GfasMS


# ORDINATION PLOTS: function .................................

my_nMDS <- function(data, ordn, title, ...) {
  phyloseq::plot_ordination(data, ordn, color = "colony_ms", shape = "state") +
    scale_shape_manual(breaks = c("bleached", "symbiotic"), values = c(1, 16)) +
    scale_color_manual(values = palette_GfasMS) + # breaks = brks, 
    labs(
      title = title,
      subtitle = paste0("Stress = ", round(ordn$stress, 3)) ) +
    stat_ellipse(aes(group = origin, linetype = origin)) + # , linetype = 2
    theme(
      panel.grid = element_blank(),
      plot.subtitle = element_text(hjust = 1)
    )
}


#### Weighted Uni Frac ----------------------

# Non-transformed data
a <- my_nMDS(
  data = ps_rarefied.1000avg,
  ordn = ordn1000_wuf,
  title = "Weighted UniFrac" )

# Relative abundance
aa <- my_nMDS(
  data = ps_rar1000_rel,
  ordn = ordn1000_wuf_relab,
  title = "Weighted UniFrac, relative abundance" )


# Sqrt-transformed data
b <- my_nMDS(
  data = ps_rar1000_sqrt,
  ordn = ordn1000_wuf_sqrt,
  title = "Weighted UniFrac, sqrt()"
)


# log(x + 1)-transformed data
cc <- my_nMDS(
  data = ps_rar1000_logx1,
  ordn = ordn1000_wuf_logx1,
  title = "Weighted UniFrac, log(x + 1)"
)



#### Bray-Curtis -----------------------------------------

# Non-transformed data
d <- my_nMDS(
  data = ps_rarefied.1000avg,
  ordn = ordn1000_bray,
  title = "Bray-Curtis"
)

# Sqrt-transformed data
dd <- my_nMDS(
  data = ps_rar1000_rel,
  ordn = ordn1000_bray_relab,
  title = "Bray-Curtis, relative abundance"
)

# Sqrt-transformed data
ee <- my_nMDS(
  data = ps_rar1000_sqrt,
  ordn = ordn1000_bray_sqrt,
  title = "Bray-Curtis, sqrt()"
)


# log(x + 1)-transformed data
f <- my_nMDS(
  data = ps_rar1000_logx1,
  ordn = ordn1000_bray_logx1,
  title = "Bray-Curtis, log(x + 1)"
)


# Patchwork .............................................
library(patchwork)

pw <- ( a + aa + b + cc ) / ( d + dd + ee + f )


(
  pw +
    plot_layout(guides = "collect")
)

ggsave("./out/Gfas_16S/beta_diversity/rarefied1000_wUniF_BrayC_nMDS_8.png", 
       bg = "white",
       dpi = 330, 
       units = "cm", width = 19, height = 27) # width = 25, height = 17


# Pretty much the same as with ps-rarefied (single subsampling) :D





# PLOT nMDS for MANUSCRIPT =============================================

# Plot nMDS in {ggplot2} (without {phyloseq})

# The trick here is to extract the values calculated for the nMDS matrix 
# and convert them to a tidy object (tibble) and then just proceed
# with regular tidyverse manipulations (ggplot!)

ordn1000_bray # is already metaMDS() type so should do ...
str(ordn1000_bray) # see that it's a list of 35 objects
ordn1000_bray$points # position of each point in the MDS space
vegan::scores(ordn1000_bray) # same as above but from {vegan}

# Base R plot - just FYI
# plot(ordn_bray)
plot(ordn1000_bray)

# Create tibble for ggplot
nMDS_bray1000 <- ordn1000_bray$points %>% 
  as_tibble(rownames = "sample_id") %>% # make into tidy object
  inner_join(., metadata_ms, by = "sample_id") # add metadata (use the one for the MS!)





## nMDS with ellipses by origin -------------------------------------------

# install.packages("ggnewscale")
library('ggnewscale')

ggplot() +
  stat_ellipse(data = nMDS_bray1000,
               aes(x = MDS1, y = MDS2,
                   color = origin_ms),
               level = 0.95,
               size = 0.5,
               show.legend = F) +
  scale_color_manual(breaks = c("Red Sea", "Hong Kong"),
                     values = c(Red_Sea, Hong_Kong),
                     guide = "none") +
  ggnewscale::new_scale_color() +
  geom_point(data = nMDS_bray1000,
             inherit.aes = F,
             aes(x = MDS1, y = MDS2,
                 color = colony_ms, 
                 shape = state,
                 size = state), 
             stroke = 1.3)  +
  scale_color_manual(values = palette_GfasMS,
                     guide = guide_legend(order = 2)) + 
  labs(
    subtitle = paste0("Stress = ", round(ordn1000_bray$stress, 3))) +
  scale_shape_manual(
    breaks = c("symbiotic", "bleached"), 
    labels = c("Symbiotic", "Bleached"), 
    values = c(16, 1),
    guide = guide_legend(order = 1)) +  
  scale_size_manual(breaks = c("symbiotic", "bleached"), 
                    values = c(2, 1.5), 
                    # guide = guide_legend(order = 1)
                    guide = "none"
                    ) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.5)) +
  coord_fixed(ratio = 1, # necessary for correct ratio of axes!
              xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme(
    panel.grid = element_blank(),
    plot.subtitle = element_text(hjust = 0.98, vjust = -10, size = 9),
    text = element_text(size = 13),
    legend.title = element_blank(),
    legend.justification = "top",
    legend.spacing = unit(0.05, 'cm'),
    legend.key.size = unit(0.3, 'cm')
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(
        size = c(2.1, 1.5))  )  )


ggsave("./out/Gfas_16S/beta_diversity/Bray_nMDS_rarefied1000_ellipses.png", 
       bg = "white",
       dpi = 330, 
       units = "cm", width = 15.5, height = 12.5)





# Calculate dist for ORIGIN (RS and HK) separately ------------

## Subset from ps_rarefied by 'colony_ms' ----------------------

# Red Sea: ps_rarefied_RS
ps_rarefied.1000avg_RedSea <- readRDS("./out/RDS_files/ps_rarefied.1000avg_RedSea.rds")

# Hong Kong: ps_rarefied_HK
ps_rarefied.1000avg_HongKong <- readRDS("./out/RDS_files/ps_rarefied.1000avg_HongKong.rds")



## Calculate Bray-Curtis distance ------------------

# Red Sea
dist1000_bray_RS = phyloseq::distance(ps_rarefied.1000avg_RedSea, method = "bray")

# Hong Kong
dist1000_bray_HK = phyloseq::distance(ps_rarefied.1000avg_HongKong, method = "bray")


## Calculate ordination points ------------------

# Red Sea
ordn1000_bray_RS = phyloseq::ordinate(ps_rarefied.1000avg_RedSea, 
                                  method = "NMDS", 
                                  distance = dist1000_bray_RS)
plot(ordn1000_bray_RS)

# Hong Kong
ordn1000_bray_HK = phyloseq::ordinate(ps_rarefied.1000avg_HongKong, 
                                  method = "NMDS", 
                                  distance = dist1000_bray_HK)
plot(ordn1000_bray_HK)

### Export ordn RDS -----------------------
saveRDS(ordn1000_bray_RS, "./out/RDS_files/ordn1000_bray_RS.rds") 
saveRDS(ordn1000_bray_HK, "./out/RDS_files/ordn1000_bray_HK.rds")


## Convert to tibbles (for ggplot) ----------------------

# Red Sea
nMDS_bray1000_RS <- ordn1000_bray_RS$points %>% 
  as_tibble(rownames = "sample_id") %>% 
  inner_join(., metadata_ms, by = "sample_id") 

# Hong Kong
nMDS_bray1000_HK <- ordn1000_bray_HK$points %>% 
  as_tibble(rownames = "sample_id") %>% 
  inner_join(., metadata_ms, by = "sample_id") 


## Export as .RDS (for plotting in 'XX_combined_plots.R') -------------

# Red Sea
saveRDS(nMDS_bray1000_RS, file = "./out/RDS_files/nMDS_bray1000_RS.rds")

# Hong Kong
saveRDS(nMDS_bray1000_HK, file = "./out/RDS_files/nMDS_bray1000_HK.rds")

# all (not subsampled)
saveRDS(nMDS_bray1000, file = "./out/RDS_files/nMDS_bray1000.rds")

# all ordn1000_bray
saveRDS(ordn1000_bray, file = "./out/RDS_files/ordn1000_bray_all.rds")





## Plot ellipses by ORIGIN - for SUPPLEMENTARY ------------------------------

# Dispersion is clearly different -> found significant difference with PERMDISP

ggplot(data = nMDS_bray1000, aes(x = MDS1, 
                             y = MDS2,
                             color = origin_ms,
                             shape = state,
                             group = origin_ms)) +
  geom_point() +
  scale_color_manual(breaks = c("Red Sea", "Hong Kong"), values = c(Red_Sea, Hong_Kong)) +
  scale_shape_manual(breaks = c("symbiotic", "bleached"), values = c(16, 1)) +
  stat_ellipse() +
  scale_x_continuous(limits = c(-0.9, 0.8), breaks = seq(-0.5, 0.5, by = 0.5)) +
  scale_y_continuous(limits = c(-0.9, 0.8), breaks = seq(-0.5, 0.5, by = 0.5)) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.justification = "top"
  )


ggsave("./out/Gfas_16S/beta_diversity/rarefied1000_Bray_nMDS_ellipses_origin.png", 
       bg = "white",
       dpi = 330, 
       units = "cm", width = 15.5, height = 12.5)



# +++++++++++++++++++++++++++++++ -------------------
# BELOW THIS: MOVE TO DEDICATED SCRIPT! ========================================




# PERMDISP with vegan::betadisper() ----------------------------------------------------------------------

dist1000_bray
metadata_physeq <- data.frame(phyloseq::sample_data(ps_rarefied.1000avg))

## By STATE -----------------------------

by_state <- vegan::betadisper(d = dist1000_bray,
                              group = metadata_physeq$state)
# Check
plot(by_state)
boxplot(by_state)
vegan::scores(by_state, display = "centroids")
vegan::scores(by_state, display = "sites")

# Test
anova(by_state) 
# Response: Distances
#           Df  Sum Sq  Mean Sq F value Pr(>F)
# Groups     1 0.00169 0.001688  0.0437  0.836
# Residuals 25 0.96459 0.038584   

# -> does not reject the Ho ... 
# so no significant difference in dispersion between 
# symbiotic and bleached 



## By ORIGIN  -----------------------------

by_origin <- vegan::betadisper(d = dist1000_bray, #  dist_bray,
                               group = metadata_physeq$origin_ms)
# Check
plot(by_origin)
boxplot(by_origin)
# scores(by_origin, display = "centroids")
# scores(by_origin, display = "sites")

# Test
anova(by_origin) 
# Response: Distances
#           Df  Sum Sq Mean Sq F value    Pr(>F)    
# Groups     1 0.33712 0.33712  118.73 5.516e-11 ***
# Residuals 25 0.07098 0.00284 

# As expected, dispersion is significanlty diff by ORIGIN!



# By COLONY ................................
boxplot(vegan::betadisper(d = dist1000_bray,
                          group = metadata_physeq$colony_ms))
anova(vegan::betadisper(d = dist1000_bray, #dist_bray,
                        group = metadata_physeq$colony_ms))
#           Df  Sum Sq  Mean Sq F value Pr(>F)
# Groups     4 0.24776 0.061941  2.0261 0.1258
# Residuals 22 0.67257 0.030571


# not significant because we can't compare state **within** colony ...
# (we are simply looking at whether colonies (including both symbiotic and bleached) 
#  have significantly different dispersion)



## By STATE for RED SEA only --------------------------------------

# Extract metadata form phyloseq
metadata_RS <- data.frame(phyloseq::sample_data(ps_rarefied.1000avg_RedSea))

# betadisper
by_state_RS <- vegan::betadisper(d = dist1000_bray_RS,
                                 group = metadata_RS$state)

# Check
plot(by_state_RS)
boxplot(by_state_RS)

# Test
anova(by_state_RS) 
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     1 0.004185 0.0041848  0.7379 0.4038
# Residuals 15 0.085063 0.0056709                

# NOT SIGNIFICANT



## By STATE for HONG KONG only --------------------------------------

# Extract metadata form phyloseq
metadata_HK <- data.frame(phyloseq::sample_data(ps_rarefied.1000avg_HongKong))

# betadisper
by_state_HK <- vegan::betadisper(d = dist1000_bray_HK, #  dist_bray_HK,
                                 group = metadata_HK$state)


# Check
plot(by_state_HK)
boxplot(by_state_HK)

# Test
anova(by_state_HK) # Pr(>F) = OLD 0.05488, NOW = 0.05652 
#           Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Groups     1 0.016926 0.0169263  4.9617 0.05652 .
# Residuals  8 0.027291 0.0034114                  

# NOT SIGNIFICANT but ALMOST ...
round(0.05652, 2) # 0.06
round(0.05652, 3) # 0.057




# PERMANOVA with vegan::ADONIS2 - whole data set ----------------------------------------------------------------------

# Cannot use original metadata because it does not match in rows with 
# ps.rarefied 'coz we dropped C1 (from beginning) and one S-60 after rarefying 
metadata_physeq <- data.frame(sample_data(ps_rarefied.1000avg)) # OLD ps_rarefied

# Specify which "dist"
# (for later to save typing)
distance <- dist1000_bray  # dist_bray

### `state + origin + state * origin` -------------------------------------------- 
## Test effect of state, origine and the combination of the two:
# distance ~ state + origin + state * origin 
vegan::adonis2(distance ~ state + origin + state * origin, 
               data = metadata_physeq, permutations = 9999) 

#              Df SumOfSqs      R2      F Pr(>F)    
# state         1   0.3199 0.04072 1.2676 0.2405    
# origin        1   1.3982 0.17797 5.5406 0.0007 ***
# state:origin  1   0.3339 0.04250 1.3231 0.2145    
# Residual     23   5.8042 0.73881                  
# Total        26   7.8562 1.00000                  


# Export results - note that they change every time coz random permutations!
# (but always < 0.0001)
vegan::adonis2(distance ~ state + origin + state * origin, 
               data = metadata_physeq, permutations = 99999) %>% 
  .[1:5] %>% as_tibble(., rownames = " ") %>% 
  mutate_if(is.numeric, round, 5) %>% # nr of decimals for P value = nr of digits of permutations
  mutate_if(is.numeric, format, scientific = FALSE) %>% 
  write_csv("./out/Gfas_16S/useful_tables/rar1000_PERMANOVA_adonis2_state-origin.csv", na = "")

# CONCLUSION:
# Bacterial communities composition differ significantly based on sample ORIGIN


### `origin` (Red Sea vs Hong Kong) -------------------------------------------- 
## Test effect of state, origine and the combination of the two:
# distance ~ state + origin + state * origin 
vegan::adonis2(distance ~ origin, 
               data = metadata_physeq, permutations = 9999) 

#          Df SumOfSqs      R2      F Pr(>F)    
# origin    1   1.4106 0.17956 5.4713  6e-04 ***
# Residual 25   6.4454 0.82044                  
# Total    26   7.8559 1.00000   

# CONCLUSION:
# Bacterial communities composition differ significantly based on sample ORIGIN,
#   (in line with previous test)


### 'state' (Symbiotic vs Bleached) -------------------------------------------- 
## Just out of curiosity: try PERMANOVA with only 'state' (Symbiotic vs Bleached)
vegan::adonis2(distance ~ state, 
               data = metadata_physeq, permutations = 99999) %>% 
  .[1:5] %>% as_tibble(., rownames = " ") %>% 
  mutate_if(is.numeric, round, 5) %>% 
  mutate_if(is.numeric, format, scientific = FALSE) %>% # view()
  write_csv("./out/Gfas_16S/useful_tables/rar1000_PERMANOVA_adonis2_state.csv", na = "")
# Pr(>F) = OLD 0.3515, NOW = 0.34837
# ` `          Df SumOfSqs R2       F         `Pr(>F)` 
# 1 state    " 1"  0.31988  0.04072 "1.06113" "0.34837"
# 2 Residual "25"  7.53631  0.95928 "     NA" "     NA"
# 3 Total    "26"  7.85619  1.00000 "     NA" "     NA"

# Again, not significant by state!



## PERMANOVA by STATE, for RED SEA only -------------

distance <- dist1000_bray_RS # dist_bray_RS

metadata_RS <- data.frame(phyloseq::sample_data(ps_rarefied.1000avg_RedSea)) 

vegan::adonis2(distance ~ state, 
               data = metadata_RS, permutations = 99999) %>% 
  .[1:5] %>% as_tibble(., rownames = " ") %>%
  mutate_if(is.numeric, round, 5) %>% # P nr of decimals = permutations digits
  mutate_if(is.numeric, format, scientific = FALSE) %>% # view()
  write_csv("./out/Gfas_16S/useful_tables/rar1000_PERMANOVA_adonis2_RedSea_state.csv", na = "")

# ` `        Df    SumOfSqs R2      F         `Pr(>F)` 
# 1 state    " 1"  0.26541  0.04968 "0.78413" "0.55784"
# 2 Residual "15"  5.07711  0.95032 "     NA" "     NA"
# 3 Total    "16"  5.34251  1.00000 "     NA" "     NA"


# Not significant (as expected!)
# so, when looking at the Red Sea samples alone, the 
# effect of state is still not significant



## PERMANOVA by STATE, for HONG KONG only -------------

distance <- dist1000_bray_HK # dist_bray_HK

metadata_HK <- data.frame(phyloseq::sample_data(ps_rarefied.1000avg_HongKong)) 

vegan::adonis2(distance ~ state, 
               data = metadata_HK, permutations = 99999) %>% 
  .[1:5] %>% as_tibble(., rownames = " ") %>%
  mutate_if(is.numeric, round, 5) %>% 
  mutate_if(is.numeric, format, scientific = FALSE) %>%
  write_csv("./out/Gfas_16S/useful_tables/rar1000_PERMANOVA_adonis2_HongKong_state.csv", na = "")

# 1 state    1     0.37605  0.34089 "4.1375" "0.00842"
# 2 Residual 8     0.72711  0.65911 "    NA" "     NA"
# 3 Total    9     1.10316  1.00000 "    NA" "     NA"

# Now it is SIGNIFICANT! (as expected from plot ... !)

# Conclusions: state has a significant effect but only for Hong Kong samples
# (not for Red Sea and not when considering them together)