library(tidyverse)
library(ggridges)
library(ggpubr)
library(ggpmisc)
library(psych)
library(drc)

## Global ggplot options -------------------------------------------------------

color_palette <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
                   "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF")

options(ggplot2.discrete.colour = color_palette)
options(ggplot2.discrete.fill = color_palette)

theme_set(theme_bw())

jitter_width = 0.1

## Load data -------------------------------------------------------------------
wh_2023 <- read_rds("data/processed/clean_2023_data.rds") |> 
  filter(site == "WH")

# Load WH hatchery data
hatchery_2023 <- read_csv("data/raw/2023 Whitehorse Egg Survival.csv") |> 
  rename(
    WH_hatchery_id = `Hatchery Fish ID`,
    WH_egg_take_date  = `Egg Take Date`,
    WH_surv_green_eyed = `% Survival Green-Eyed`,
    WH_surv_green_hatch = `% Survival: Green - Hatch`,
    WH_surv_green_pond = `% Survival: Green - Ponding`,
    WH_female_notes = `Female Notes`,
    pre_eyed_morts = `Pre-eyed Morts`,
    tot_eyed_eggs = `Total Eyed Eggs`,
    egg_donations = `Green Or Eyed Egg Donations`,
    total_hatched = `Total at Hatched Stage`,
    hatch_ponding_morts = `Hatch to Ponding Morts`,
    tot_ponding = `Total at Ponding`
  ) |> 
  dplyr::select(
    WH_hatchery_id,
    WH_egg_take_date, 
    WH_surv_green_eyed, 
    WH_surv_green_hatch,
    WH_surv_green_pond, 
    WH_female_notes,
    pre_eyed_morts,
    tot_eyed_eggs,
    egg_donations,
    total_hatched,
    hatch_ponding_morts,
    tot_ponding
  ) |> 
  drop_na(
    WH_hatchery_id
  ) |> 
  mutate(
    tot_egg = tot_eyed_eggs + pre_eyed_morts,
    tot_survived = tot_ponding,
    tot_morts = tot_egg - tot_survived,
    surv_pct = tot_survived / tot_egg
  )

# Join
wh_2023 <- wh_2023 |> 
  left_join(hatchery_2023)

# Make long
wh_2023_long <- wh_2023 |> 
  select()

# Examine distribution of variables --------------------------------------------
## Egg nmol T


ggplot(data = wh_2023, aes(x = nmol_T_egg)) +
  geom_density() +
  geom_rug()




    
##
ggplot(data = wh_2023) +
  geom_point(aes(x = nmol_T_egg, y = surv_pct),
             size = 3)

# 2 param model

model_fit <- drm(
  surv_pct ~ nmol_T_egg,
  weights = tot_egg,
  data = wh_2023,
  fct = LL.2(),
  type = "binomial"
)

summary(model_fit)


