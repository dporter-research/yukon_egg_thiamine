library(tidyverse)
library(ggridges)

## Global ggplot options -------------------------------------------------------

color_palette <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
                   "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF")

options(ggplot2.discrete.colour = color_palette)
options(ggplot2.discrete.fill = color_palette)

theme_set(theme_bw())

jitter_width = 0.1

## Load data -------------------------------------------------------------------

# Read in clean 2023 data and filter for eggs collected at PIST
pist_2023_df <- read_rds("data/processed/clean_2023_data.rds") |> 
  drop_na(nmol_T_g) |> 
  filter(site == "PIST")

## Summary table for different genetic identifiers -----------------------------

# Make a long df
combined_counts <- pist_2023_df |> 
  select(SIN, country, broad_group, fine_group) |> 
  pivot_longer(
    cols = -SIN,
    names_to = "grouping_level",
    values_to = "group_name"
  ) |> 
  count(grouping_level, group_name, name = "n_individuals") |> 
  filter(!is.na(group_name))

combined_counts

saveRDS(combined_counts, file = "output/tables/PIST_genetic_group_counts.rds")

# Looks like maybe I can start at broad group

## Density Plot ----------------------------------------------------------------

# Filter out NAs and No Assignments
broad_groups_df <- pist_2023_df |> 
  drop_na(broad_group, g_egg) |> 
  filter(broad_group != "No Assignment") |> 
  select(SIN, broad_group, g_egg, nmol_T_egg, lipid_g, water_g, protein_g_est)

# Make a long df for plotting
long_broad_groups_df <- broad_groups_df |> 
  pivot_longer(
    cols = -c(SIN, broad_group),
    names_to = "Measurement",
    values_to = "Value"
  ) |> 
  drop_na(Value)

# Ridgeline plot
pist_density_plot <- ggplot(data = long_broad_groups_df,
                            aes(x = Value, y = broad_group,
                                fill = broad_group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  guides(fill = "none") +
  facet_wrap(~ Measurement, scales = "free")

pist_density_plot 

ggsave(
  filename = "output/figures/PIST_facetted_densities.png", 
  plot = pist_density_plot,                         
  width = 7,                                        
  height = 5,                                      
  dpi = 300                                        
)
