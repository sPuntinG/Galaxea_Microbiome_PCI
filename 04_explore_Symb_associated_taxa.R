
# See if I can find any of the bacterial taxa reported to be 
#  associated with Symbiodiniaceae in other studies, such as:
#  - Lawson 2018 (Symbiodiniaceae core microbiome) https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12599
#  - Nitschke 2020 (Symbiolites) https://www.nature.com/articles/s41396-020-0629-z
#  - Maire 2022 https://www.nature.com/articles/s41396-020-0629-z


library(tidyverse)
library(here)

# Note that this is non-rarefied data, and that
# "decont" means that the ASVs considered lab contaminants have been removed

decont_long <- readRDS("./out/RDS_files/decont_long.rds")

source("colors_GfasMS.R")


# Import data ------------------------------

lawson2018 <- read_csv(
  "./in/Symb_associated_bacteria/Symb_associated_core_OTUs_Lawson2018_cleaned.csv")

nitschke2020 <- read_csv(
  "./in/Symb_associated_bacteria/Symb_associated_OTUs_Nitschke2020.csv")


# Lawson et al. 2018 -------------------------

# Extract vector of names and replace "UC " with "unclassif_" so that it 
# matches with my nomenclature system
lawson <- lawson2018 %>% 
  pull(Taxonomic_ID) %>% 
  unique() %>% 
  str_replace(., "UC ", "unclassif_")

# Retain only what matches with lawson
decont_long_lawson <- decont_long %>% 
  filter(taxon_name %in% lawson) #%>% view()

# Inspect:
decont_long_lawson %>% pull(taxon_level) %>% unique() # "Genus"
         # all matches at the "Genus" level

decont_long_lawson %>% pull(taxon_name) %>% unique() 
         # "Alteromonas" "Pseudohongiella" "Labrenzia" "Marinobacter"
         # so 4 out of the 21 taxa reported by Lawson are found in my data set

decont_long_lawson %>% pull(state_colony_frgmt_ms) %>% unique() %>% length() #27
         # matches in 27 out of 28 samples


## Plot (bubble) by individual sample -----------------------------

# Under the same genus there are several spp and therefore I see repeated rows
# of the same genus, for example Alteromonas
decont_long_lawson %>% 
  filter(taxon_name == "Alteromonas") #%>% view()

# Therefore I need to first 'collapse' these by lumping them together
decont_long_lawson <- decont_long_lawson %>% 
  group_by(sample_id, taxon_name) %>% 
  mutate(
    count_sample = sum(count_sample), 
    rel_abund = sum(rel_abund)
  ) %>% 
  unique() #%>% view()



# Nitschke et al. 2020 -------------------------

# Let's first try with full spp names 
nitschke <- nitschke2020 %>% pull(Species) %>% unique() # There's 8 spp

decont_long %>% 
  filter(taxon_name %in% nitschke) # no matches

# Let's try by genus ...
# (as spp assignment with 16S is not reliable, and it might vary based on the method)
nitschke2020 <- nitschke2020 %>% 
  separate(Species, into = c("Genus", "species"), 
           sep = " ", remove = F) 

nitschke <- nitschke2020 %>% 
  pull(Genus) %>% unique()

decont_long %>% 
  filter(taxon_name %in% nitschke) # so now there are 108 matches


decont_long %>% 
  filter(taxon_name %in% nitschke) %>% 
  pull(taxon_name) %>% unique()
# "Alteromonas" "Labrenzia" "Thalassospira" "Pseudoalteromonas" "Marinobacter" 
# 5 genera 
# Note that these were also found before from Lawson:
#       "Alteromonas" "Labrenzia" "Marinobacter"
#   while these are "new": "Thalassospira" "Pseudoalteromonas" 



# Maire et al. 2021 --------------------------------------

# Import from excel ...in batch (faster)
path_maire <- "./in/Symb_associated_bacteria/41396_2021_902_MOESM4_ESM.xlsx"
tab_names <- readxl::excel_sheets(path = path_maire)
list_all <- lapply(tab_names, 
                   function(x) readxl::read_excel(path = path_maire, sheet = x))

# Clean up table to keep only taxonomy 
# (for now not interested in their relative abundance)
maire2021 <- map_dfr(list_all, bind_rows) %>% 
  drop_na() %>% 
  select(Phylum:Genus) %>% 
  mutate_all(~str_replace(., "^[:alpha:]_[:digit:]__", ""))

maire <- maire2021 %>% pull(Genus) %>% unique() # 15 genera

decont_long %>% 
  filter(taxon_name %in% maire) %>% 
  pull(taxon_name) %>% unique()  # 7 matches
# "uncultured" "Staphylococcus" "Pseudohongiella" "Pseudomonas" "Labrenzia"
# "SM1A02" "Hyphomicrobium" 

decont_long_maire <- decont_long %>% 
  filter(taxon_name %in% maire) #%>% view()

# Note that there are many matches with "uncultured" that is not useful info ...
# remove
decont_long_maire <- decont_long_maire %>% 
  # filter(taxon_name == "uncultured") # 185 rows
  filter(taxon_name != "uncultured")




# All together ----------------------------------------

# Make summary tables
this_study <- decont_long %>% 
  filter(taxon_level == "Genus") %>% 
  pull(taxon_name) %>% unique() #%>% length() # 135

this_study_symb <- decont_long %>% 
  filter(state == "symbiotic") %>% 
  filter(taxon_level == "Genus") %>% 
  pull(taxon_name) %>% unique() #%>% length() # 94

this_study_blea <- decont_long %>% 
  filter(state == "bleached") %>% 
  filter(taxon_level == "Genus") %>% 
  pull(taxon_name) %>% unique() #%>% length() # 102

# Create a tibble for each study
tibble_this <- tibble("Genus" = this_study) %>%
          mutate(study = "this_study")

tibble_this_symb <- tibble("Genus" = this_study_symb) %>%
  mutate(study = "this_study_symb")

tibble_this_blea <- tibble("Genus" = this_study_blea) %>%
  mutate(study = "this_study_blea")


tibble_lawson <- tibble("Genus" = lawson) %>%
  mutate(study = "lawson")

tibble_nitschke <- tibble("Genus" = nitschke) %>%
  mutate(study = "nitschke")

tibble_maire <- tibble("Genus" = maire) %>%
  mutate(study = "maire")


# Put together all tibbles
list_tibbles_6 <- list(tibble_this, tibble_this_symb, tibble_this_blea,
                       tibble_lawson, tibble_nitschke, tibble_maire)

list_tibbles_5 <- list(tibble_this_symb, tibble_this_blea, # tibble_this,
                       tibble_lawson, tibble_nitschke, tibble_maire)

list_tibbles_4 <- list(tibble_this, #tibble_this_symb, tibble_this_blea,
                       tibble_lawson, tibble_nitschke, tibble_maire)

all6 <- map_dfr(list_tibbles_6, bind_rows)
all5 <- map_dfr(list_tibbles_5, bind_rows)
all4 <- map_dfr(list_tibbles_4, bind_rows)


# Clean up some space
rm(list = ls(pattern = "tibble"))

# Filter this study to keep only genera (Genus) that are present in 
# at least one of the other studies
other_studies <- c(lawson, nitschke, maire) %>% unique() %>% # 37 elements
  .[!. == 'uncultured'] #%>% length() # 36

all6 <- all6 %>% 
  filter(Genus %in% other_studies)

all5 <- all5 %>% 
  filter(Genus %in% other_studies)

all4 <- all4 %>% 
  filter(Genus %in% other_studies)


# score (with "n") based on occurrence (5 = found in all studies)
scored5 <- all5 %>% 
  group_by(Genus) %>% 
  summarise(n = n()) %>% #view()
  right_join(all5, ., by = "Genus") %>% 
  mutate(study = factor(study,
                        levels = c("this_study", "this_study_symb",
                                   "this_study_blea",
                                   "lawson", "nitschke", "maire")))

# score (with "n") based on occurrence (4 = found in all studies)
scored4 <- all4 %>% 
  group_by(Genus) %>% 
  summarise(n = n()) %>% #view()
  right_join(all4, ., by = "Genus") %>% 
  mutate(study = factor(study,
                        levels = c("this_study", "this_study_symb",
                                   "this_study_blea",
                                   "lawson", "nitschke", "maire")))


# Make nicer labels for groups (y axis) 
# 1. Make tibble
study <- all6$study %>% unique()
study_label <- c("This study", 
                 "This study: Symbiotic polyps", "This study: Bleached polyps",
                 "Lawson et al. 2018", 
                 "Nitschke et al. 2020", 
                 "Maire et al. 2021")

# in case I wanted to add some info on the other studies ...
# BUT then problem with text size and overlapping ... !
study_label_2 <- c("This study", 
                 "This study: Symbiotic polyps", "This study: Bleached polyps",
                 "Lawson et al. 2018\n(cultured Symbiodiniaceae)", 
                 "Nitschke et al. 2020\n(cultured Symbiodiniaceae, Symbiolite)", 
                 "Maire et al. 2021\n(freshly isolated Symb. from G. fascicularis; intra- extracell.)")

study_order <- c(1, 2, 3, 4, 5, 6)

studylabels <- tibble(study, study_label, study_order)

# 2. Join to data for plotting
left_join(scored5, studylabels, by = "study")


## Plot 5 groups (this symb, this bleach + 3 others) ----------------------
left_join(scored5, studylabels, by = "study") %>% 
  ggplot(data = ., 
         aes(x = forcats::fct_reorder(Genus, -n), 
             y = forcats::fct_reorder(study_label, -study_order))) +
  geom_point(aes(color = n)) +
  scale_color_viridis_c(guide = "none") +
  scale_x_discrete(position = "top") +
  # labs(caption = "Comparison between taxa reported in studies on Symbiodiniaceae associated bacteria and this study.\nNote: only matches at the Genus level reported (no matches found at species level). Color-coded by nr of intersections.") +
  theme_light() +
  theme(
    panel.border = element_rect(colour = NA),
    axis.text = element_text(size = 5, color = "black"), 
    axis.text.x = element_text(size = 5, angle = 45, hjust = 0), 
    axis.ticks = element_blank(),
    plot.caption.position = "plot",
    # plot.caption = element_text(size = 4, hjust = 0),
    axis.title = element_blank(),
    plot.margin = margin(t = 0, r = 2.5, l = 0.5, b = 0, "cm"),
    aspect.ratio = 5/37 # make "square" space: 5 gr. on y and 37 genera on x
  ) 
  
# Plot without comments to other studies
ggsave("./out/Gfas_16S/Symb_associated_taxa/overview_presabs_5.png",
       bg = "white",
       dpi = 330,
       type = 'cairo', # doesn't really help with the weird letter aligning ...
       units = "cm", width = 14, height = 5)
