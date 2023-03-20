
source("03_check_data_structuring2.R")
source("colors_GfasMS.R")

library("eulerr")


# Euler/Venn diagram ========================================================

# Euler diagram is a Venn with circles size proportional to their "content"

# Count by ASVs
filtered <- decont_long %>% 
  filter(taxon_level == "ASV_id")

symbiotic <- filtered %>% 
  filter(state == "symbiotic") %>% 
  pull(taxon_name) %>% unique()

symbiotic_redsea <- filtered %>% 
  filter(state == "symbiotic" & origin_ms == "Red Sea") %>% 
  pull(taxon_name) %>% unique()

symbiotic_hongkong <- filtered %>% 
  filter(state == "symbiotic" & origin_ms == "Hong Kong") %>% 
  pull(taxon_name) %>% unique()


bleached <- filtered %>% 
  filter(state == "bleached") %>% 
  pull(taxon_name) %>% unique()

bleached_redsea <- filtered %>% 
  filter(state == "bleached" & origin_ms == "Red Sea") %>% 
  pull(taxon_name) %>% unique()

bleached_hongkong <- filtered %>% 
  filter(state == "bleached" & origin_ms == "Hong Kong") %>% 
  pull(taxon_name) %>% unique()



redsea <- filtered %>% 
  filter(origin_ms == "Red Sea") %>% 
  pull(taxon_name) %>% unique()

hongkong <- filtered %>% 
  filter(origin_ms == "Hong Kong") %>% 
  pull(taxon_name) %>% unique()



# Plot - by ASVs - by STATE and by ORIGIN (separatedly) .....................................................

windowsFonts(A = windowsFont("Times New Roman"))
windowsFonts(C = windowsFont("Calibri"))

# by STATE
by_state <- eulerr::euler(list("    Symbiotic\n  (n = 14)" = symbiotic, "Bleached\n(n = 14)" = bleached))

ASV_state <- base::plot(
  by_state,
  col = c("NA", "black"),
  fills = c("black", "white"), 
  labels = list(col = c("white", "black")),
  quantities = list(
    col = c("white", "black", "black"),
    font = 3),
  family = "C" )


(ggplotized_state <- patchwork::wrap_elements(ASV_state) )

ggsave("./out/Gfas_16S/supplementary/Eulerr_ASV_state.png",
       bg = "white", 
       dpi = 330,
       units = "cm", height = 11, width = 11)


# by ORIGIN
by_origin <- eulerr::euler(list("Red Sea\n(n = 18)" = redsea, 
                                "Hong Kong\n(n = 10)" = hongkong))

(
  ASV_origin <- base::plot(
  by_origin,
  fills = c(Red_Sea, Hong_Kong),
  col = "NA",
  labels = list(col = c("black", "black")),
  quantities = list(
    col = c("black", "black", "black"),
    font = 3),
  family = "C" )
)

(ggplotized_origin <- patchwork::wrap_elements(ASV_origin) )

ggsave("./out/Gfas_16S/supplementary/Eulerr_ASV_origin.png",
       bg = "white", 
       dpi = 330,
       units = "cm", height = 11, width = 11)



# Euler 4-ways -----------------------------------------------------------------
fourways_ASV <- eulerr::venn(list(  # eulerr::euler()
  "Red Sea\nBleached" = bleached_redsea,
  "Hong Kong\nBleached" = bleached_hongkong,
  "Red Sea\nSymbiotic" = symbiotic_redsea,
  "Hong Kong\nSymbiotic" = symbiotic_hongkong
)) 

fourways_ASV_plot <- base::plot(fourways_ASV,
                                shape = "ellipse",
                                lwd = 3,
                                quantities = T,
                                labels = F,
                                # fills = c("#D84E1A30", "#2670E430", Red_Sea, Hong_Kong),
                                fills = c("white", "white", Red_Sea, Hong_Kong),
                                col = c(Red_Sea, Hong_Kong, Red_Sea, Hong_Kong) )

fourways_ASV_plot


# Export with ggsave() (I really don't like base R for that!)
ggplotized <- patchwork::wrap_elements(fourways_ASV_plot)

ggsave("./out/Gfas_16S/supplementary/Venn_4ways_ASV_nolabs.png",
       plot = ggplotized,
       bg = "white",
       dpi = 330, 
       units = "cm", width = 11, height = 11)

# Will have to adjust the labels manually later 
# (easier - and maybe not even necessary)
