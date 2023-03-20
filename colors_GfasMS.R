
# Color palette for "origin_ms"
palette_GfasMS <- c( "RS1" = "#F5A413", "RS2" = "#EF4719", "RS3" = "#A50021",
                     "HK1" = "#15AEF3", "HK2" = "#3632D6" )


# # This is silenced so I can source() without interfering with the running script
# RGBs <- col2rgb(palette_GfasMS)
# 
# library("tidyverse")
# 
# 
# middle <- as_tibble(RGBs) %>%
#   mutate(
#     rs = round((RS1 + RS2 + RS3)/3, 0),
#     hk = round((HK1 + HK2)/2, 0),
#     all = round((RS1 + RS2 + RS3 + HK1 + HK2)/5, 0)
#   )
# 
# rs_rgb <- middle %>% pull(rs)
# hk_rgb <- middle %>% pull(hk)
# all_rgb <- middle %>% pull(all)


Red_Sea <- rgb(216, 78, 26, maxColorValue = 255) # , names = "Red_Sea"
Hong_Kong <- rgb(38, 112, 228, maxColorValue = 255) # , names = "Hong_Kong"
all <- rgb(145, 92, 107, maxColorValue = 255)

# scales::show_col(c(Red_Sea, Hong_Kong, all))
