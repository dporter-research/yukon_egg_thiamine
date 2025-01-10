# 2023 Yukon Chinook egg thiamine data tidying
# D. Porter - 2024-11-19

# The purpose of this script is to add site data to the 2023 Yukon Chinook
# egg thiamine dataset. Site data is embedded in the SIN variable, I am just
# extracting it into its own column. I am also adding a "group" factor to
# distinguish the river entry, transit, and spawning sites. I am also adding
# a year column in case this dataset gets incorporated into a larger set.
# I am exporting the cleaned data as csv for general purpose and as rds
# to read in to future R scripts while maintaining factor levels established
# here

library(tidyverse)
library(readxl)
library(here)

data <- read_xlsx(path = here("data", "yukon23eggT.xlsx"), na = c("NA", "X", "x"))

# Correct nmol T/g for dry mass
data <- data |> 
  mutate(nmol_T_g_dw_corrected = (nmol_T_g / ((100 - moisture) / 100))) |> 
  relocate(nmol_T_g_dw_corrected, .after = nmol_T_g)


data <- data |> 
  mutate(year = factor("2023")) |> 
  mutate(site = factor(case_when(
    str_detect(SIN, "Chena") ~ "Chena",
    str_detect(SIN, "WH") ~ "WH",
    str_detect(SIN, "FOYU") ~ "FOYU",
    str_detect(SIN, "PIST") ~ "PIST",
    str_detect(SIN, "RARA") ~ "RARA",
    str_detect(SIN, "Salcha") ~ "Salcha",
    TRUE ~ "Other"  # Default value if none match
  ),
  levels = c("PIST", "RARA", "FOYU", "Chena", "Salcha", "WH")
  )) |> 
  mutate(group = factor(case_when(
    site %in% c("Chena", "Salcha", "WH") ~ "Spawning",
    site %in% "PIST" ~ "Entry",
    site %in% c("FOYU", "RARA") ~ "Transit",
    .default = "Other"
  ),
  levels = c("Entry", "Transit", "Spawning")
  ))

write_csv(data, file = here("data", "yukon_egg_thiamine_2023.csv"))
write_rds(data, file = here("data", "yukon_egg_thiamine_2023.rds"))
