
# Script dedicated to exploring presence and abundance of 
# Endozoicomonadaceae


source("04_organize_Family_level.R")


# How many symb and bleached samples? -----------------------

# Symbiotic
fam %>% 
  filter(state == "symbiotic") %>% 
  pull(state_colony_frgmt_ms) %>% 
  unique() # %>% length() # 14

fam %>% 
  filter(state == "symbiotic" & family == "Endozoicomonadaceae") %>% 
  pull(state_colony_frgmt_ms) %>% 
  unique() %>% length() # 7

  # 7/14 symb samples have Endozoicomonadaceae


# Bleached
fam %>% 
  filter(state == "bleached") %>% 
  pull(state_colony_frgmt_ms) %>% 
  unique()  %>% length() # 14

fam %>% 
  filter(state == "bleached" & family == "Endozoicomonadaceae") %>% 
  pull(state_colony_frgmt_ms) %>% 
  unique() %>% length() # 8

  # 8/14 symb samples have Endozoicomonadaceae


# How many samples with Endozoicomonadaceae ---------------------

# Keep only Endozoicomonadaceae ----------------
endozoicomonad <- fam %>% 
  filter(family == "Endozoicomonadaceae") %>% 
  arrange(desc(state), state_colony_frgmt_ms) # %>% view()

# 15 rows = present in only 15 (out of 28) samples


# Mean relative abundance of Endozoicomonadaceae in the data set --------------------

## Mean rel. ab. overall -------------
endozoicomonad %>% 
  pull(rel_abund_bysample) %>%
  sort() %>%  
  "*"(100) %>% 
  mean(., na.rm = T) %>% # 5.78 %
  round(., 2) %>% 
  # range() %>% # "0.06 %"  "39.12 %"
  paste(., "%") 


## Mean rel. ab. in Symbiotic only ------------------
endozoicomonad %>% 
  filter(state == "symbiotic") %>% 
  pull(rel_abund_bysample) %>%
  sort() %>%  
  "*"(100) %>% 
  mean(., na.rm = T) %>% # "2.25 %"
  # range() %>% # "0.3 %"  "7.58 %"
  round(., 2) %>% 
  paste(., "%")  



## Mean rel. ab. in Bleached only ---------------
endozoicomonad %>% 
  filter(state == "bleached") %>% 
  pull(rel_abund_bysample) %>%
  sort() %>%  
  "*"(100) %>% 
  # mean(., na.rm = T) %>% # "8.87 %"
  range() %>%              # "0.06 %"  "39.12 %"
  round(., 2) %>% 
  paste(., "%")  


