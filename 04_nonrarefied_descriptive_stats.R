
# Descriptive stats for the Results section
# Non-rarefied data

source("03_check_data_structuring2.R")

source("colors_GfasMS.R")


# Families stats ------------------------------------------------

# How many Families ...
## Abundance ------------
### Overall (across whole data set) ------------

decont_long %>% 
  filter(taxon_level == "Family") %>% 
  pull(taxon_name) %>% unique() %>% sort() %>% 
  length() # 119
           # see that it includes "uncultured" and "Unclassified" which (potentially)
           # lumps different taxa together

decont_long %>% 
  filter(taxon_level == "uncX_FAM") %>% 
  pull(taxon_name) %>% unique() %>% sort() %>% 
  length() # 138
           # there's more now that we use uncX_FAM which is more accurate

fams <- decont_long %>% filter(taxon_level == "uncX_FAM")


### in SYMBIOTIC (all, RS & HK together) ----------

fams %>% 
  filter(state == "symbiotic") %>% 
  pull(taxon_name) %>% unique() %>%
  length() # 104

decont_long %>% 
  filter(state == "symbiotic") %>% 
  filter(taxon_level == "uncX_FAM") %>% 
  pull(taxon_name) %>% unique() %>%
  length() # 104

#### in SYMBIOTIC Red Sea ----------
decont_long %>% 
  filter(state == "symbiotic" & origin_ms == "Red Sea") %>% 
  filter(taxon_level == "uncX_FAM") %>% 
  pull(taxon_name) %>% unique() %>%
  length() # 74

#### in SYMBIOTIC Hong Kong ----------
decont_long %>% 
  filter(state == "symbiotic" & origin_ms == "Hong Kong") %>% 
  filter(taxon_level == "uncX_FAM") %>% 
  pull(taxon_name) %>% unique() %>%
  length() # 66



### in BLEACHED (all, RS & HK together) ----------
decont_long %>% 
  filter(state == "bleached") %>% 
  filter(taxon_level == "uncX_FAM") %>% 
  pull(taxon_name) %>% unique() %>%
  length() # 109

#### in BLEACHED Red Sea ----------
decont_long %>% 
  filter(state == "bleached" & origin_ms == "Red Sea") %>% 
  filter(taxon_level == "uncX_FAM") %>% 
  pull(taxon_name) %>% unique() %>%
  length() # 89

#### in BLEACHED Hong Kong ----------
decont_long %>% 
  filter(state == "bleached" & origin_ms == "Hong Kong") %>% 
  filter(taxon_level == "uncX_FAM") %>% 
  pull(taxon_name) %>% unique() %>%
  length() # 69


## Ranking -------------------------

### Overall ------------------------
overall <- fams %>% 
  group_by(taxon_name) %>% 
  summarise(nr_reads = sum(count_sample)) %>% 
  arrange(desc(nr_reads)) %>% 
  mutate(
    perc = nr_reads / sum(nr_reads) * 100,
    perc_cumsum = cumsum(perc),
    rank = row_number()) %>% 
  # pull(perc) %>% sum() # 100 # sanity check
  rename(Family = taxon_name) %>% 
  mutate_at(vars(starts_with("perc")), round, 2) #%>% view()

# Write out
write_csv(overall,
          "./out/Gfas_16S/nonrarefied_descriptive_stats/Family_overall.csv")

# Plot curve
ggplot(data = overall, 
       aes(x = rank, y = perc_cumsum)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 100), 
                     labels = scales::percent_format(scale = 1),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 150, 10), expand = c(0, 0)) +
  xlab("nr of Families") +
  ylab("Cumulative abundance\n(overall)") +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),
    panel.grid.major.y = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.major.x = element_line(colour = "grey75", linetype = "dotted")
  )

ggsave("./out/Gfas_16S/nonrarefied_descriptive_stats/Family_cumabund_overall.png",
       units = "cm",
       width = 17, height = 10)


### Symbiotic ------------------------
symb <- fams %>% 
  filter(state == "symbiotic") %>% 
  group_by(taxon_name) %>% 
  summarise(nr_reads = sum(count_sample)) %>% 
  arrange(desc(nr_reads)) %>% 
  mutate(
    perc = nr_reads / sum(nr_reads) * 100,
    perc_cumsum = cumsum(perc),
    rank = row_number()) %>% 
  # pull(perc) %>% sum() # 100 # sanity check
  rename(Family = taxon_name) %>% 
  mutate_at(vars(starts_with("perc")), round, 2)  %>% 
  write_csv(.,
          "./out/Gfas_16S/nonrarefied_descriptive_stats/Family_symbiotic.csv")


### Bleached ------------------------
blea <- fams %>% 
  filter(state == "bleached") %>% 
  group_by(taxon_name) %>% 
  summarise(nr_reads = sum(count_sample)) %>% 
  arrange(desc(nr_reads)) %>% 
  mutate(
    perc = nr_reads / sum(nr_reads) * 100,
    perc_cumsum = cumsum(perc),
    rank = row_number()) %>% 
  # pull(perc) %>% sum() # 100 # sanity check
  rename(Family = taxon_name) %>% 
  mutate_at(vars(starts_with("perc")), round, 2) %>%
  write_csv(.,
          "./out/Gfas_16S/nonrarefied_descriptive_stats/Family_bleached.csv")

