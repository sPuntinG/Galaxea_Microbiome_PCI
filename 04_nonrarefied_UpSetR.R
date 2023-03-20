
library(tidyverse)
library(here)

# Notes:
# Problem exporting the plot with ggsave() because UpSetR is not a ggplot object.
# I tried to use the GitHub version instead of the CRAN, but still not better ...
# (for this I used remotes::install_github("hms-dbmi/UpSetR")

# install.packages("UpSetR") # This is from CRAN but I need the GitHub version
# remotes::install_github("hms-dbmi/UpSetR") # this one worked
# devtools::install_github("hms-dbmi/UpSetR") # this one wasn't working
library(UpSetR)


# Import data 
decont_long <- readRDS("./out/RDS_files/decont_long.rds")
metadata <- read_csv("./in/metadata_ms.csv")

# Color palette
source("colors_GfasMS.R")
 

# Prep data for UpSetR --------------------------

# Filter to AVS level
data <- decont_long %>% 
  filter(taxon_level == "ASV_id") #%>% view()

rm(decont_long)


# Create upsetr group for origin x state
data <- data %>% 
  mutate(group_upsetr = paste(origin_ms, state, sep = "_")) %>% 
  # mutate(group_upsetr = str_replace(group_upsetr, "Red Sea", "RS")) %>% 
  # mutate(group_upsetr = str_replace(group_upsetr, "Hong Kong", "HK")) %>% 
  mutate(group_upsetr = str_replace(group_upsetr, " ", "")) #%>% view()

data$group_upsetr %>% unique()


# Input from list: create lists

RS_symb <- data %>% 
  filter(group_upsetr == "RedSea_symbiotic") %>% 
  pull(taxon_name) %>% unique()

RS_blea <- data %>% 
  filter(group_upsetr == "RedSea_bleached") %>% 
  pull(taxon_name) %>% unique()


HK_symb <- data %>% 
  filter(group_upsetr == "HongKong_symbiotic") %>% 
  pull(taxon_name) %>% unique()

HK_blea <- data %>% 
  filter(group_upsetr == "HongKong_bleached") %>% 
  pull(taxon_name) %>% unique()

# Make list 
listInput <- list(RS_symb, RS_blea, HK_symb, HK_blea)

# Name list items (necessary for UpSetR)
names(listInput) <- c('Red Sea symbiotic', 'Red Sea bleached',
                      'Hong Kong symbiotic', 'Hong Kong bleached')


# Plot! ---------------------------------------

( upsetr <- UpSetR::upset(fromList(listInput), 
                        sets = rev(c('Red Sea symbiotic', 'Red Sea bleached',
                                     'Hong Kong symbiotic', 'Hong Kong bleached')),
                        keep.order = T, # keep order of sets()
                        order.by = "freq",
                        mb.ratio = c(0.65, 0.35),
                        point.size = 4, line.size = 1,
                        sets.x.label = "Group size\n(nr of ASVs)",
                        mainbar.y.label = "Intersections size\n(nr of ASVs)",
                        # set_size.show = T,
                        # text.scale = 2,
                        text.scale = c(1.1, #intersection size title, 
                                       1.3, #intersection size tick labels, 
                                       1.1, #set size title, 
                                       1.3, #set size tick labels, 
                                       1.6, #set names, 
                                       1.3), # numbers above bars
                  
                        queries = list(
                          list(
                            query = intersects, 
                            params = list("Red Sea symbiotic", "Red Sea bleached",
                                          "Hong Kong symbiotic", "Hong Kong bleached"), 
                            color = "purple",
                            active = T
                          )
                        )
) 
)

### Export .PNG --------------
# 1. Specify device and path
png("./out/Gfas_16S/UpSetR/upsetr.png", 
    width = 21, height = 13, units = "cm",
    pointsize = 13,
    bg = "white",  res = 300
)
# 2. Plot 
upsetr 
# 3. 
dev.off()


### Export .SVG ----------------
# 1. Specify device and path
svg("./out/Gfas_16S/UpSetR/upsetr.svg", 
    width = 8.26772, # inches = 21 cm 
    height = 5.11811, # 13 cm 
    # units = "cm", works only in inches ...
    pointsize = 13,
    bg = "white" 
    # res = 300
    # fallback_resolution = 330
)
# 2. Plot 
upsetr 
# 3. 
dev.off()



# TO DO in Inkscape (much faster, else would need {Complex Upset} ... )
# - move y axis label closer to the axis 
# - change colors in sets and points
