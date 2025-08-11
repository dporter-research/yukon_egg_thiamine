library(tidyverse)
library(ggridges)
library(ggpubr)
library(ggpmisc)
library(psych)

## Global ggplot options -------------------------------------------------------

color_palette <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
                   "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF")

options(ggplot2.discrete.colour = color_palette)
options(ggplot2.discrete.fill = color_palette)

theme_set(theme_bw())

jitter_width = 0.1

## Load data -------------------------------------------------------------------
all_2023 <- read_rds("data/processed/clean_2023_data.rds") |> 
  drop_na(g_egg)

pist <- all_2023 |> 
  filter(site == "PIST")

quantile(pist$g_egg)

mid_50_pist <- pist |> 
  filter(g_egg > 0.06 & g_egg < 0.09) 

describe(mid_50$g_egg)

quantile(all_2023$g_egg)

ggplot(data = mid_50_pist, aes(x = nmol_T_egg, y = nmol_T_g)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_poly_eq(method = "lm", use_label("eq", "R2"),
               formula = y ~ x)

pist_fit <- lm(nmol_T_g ~ nmol_T_egg, data = mid_50_pist)
summary(pist_fit)


# Read in clean 2023 data and filter for eggs collected at PIST
ich_2023_df <- read_rds("data/processed/clean_2023_data.rds") |> 
  drop_na(ich_visual_positive) 

ggplot(data = ich_2023_df, aes(x = length_nose_fork_mm, y = length_mideye_fork_mm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_poly_eq(method = "lm", use_label("eq", "R2"),
               formula = y ~ x)

ich_2023_no0 <- ich_2023_df |> 
  filter(ich_visual_count > 0)

ggplot() +
  geom_point(data = ich_2023_df,
             aes(x = ich_visual_count, y = nmol_T_egg))

ggplot() +
  geom_violin(data = ich_2023_df,
             aes(x = ich_visual_positive, y = g_egg, fill = ich_visual_positive)) +
  geom_jitter(data = ich_2023_df,
              aes(x = ich_visual_positive, y = g_egg),
              width = jitter_width) +
  guides(fill = "none") +
  facet_wrap(~ site)

ggplot() +
  geom_point(data = ich_2023_df,
             aes(x = length_mideye_fork_mm, y = nmol_T_egg))

summary_table <- ich_2023_df |> 
  group_by(site) |> 
  summarize(
    ich_pos = sum(ich_visual_positive, na.rm = TRUE),
    ich_neg = sum(!ich_visual_positive, na.rm = TRUE),
    total_fish = ich_pos + ich_neg,
    percent_pos = ich_pos / total_fish
  )

summary_table


all_ich_data <- read_csv("data/raw/2023_Rapids.csv") |> 
  mutate(ich_visual_positive = case_when(
    `Ich Visual Positive` == "Yes" ~ TRUE,
    `Ich Visual Positive` == "No" ~ FALSE,
    .default = NA
  ))

summary_2 <- all_ich_data |> 
  group_by(Site) |> 
  summarize(
    ich_pos = sum(ich_visual_positive, na.rm = TRUE),
    ich_neg = sum(!ich_visual_positive, na.rm = TRUE),
    total_fish = ich_pos + ich_neg,
    percent_pos = ich_pos / total_fish
  )

summary_2
