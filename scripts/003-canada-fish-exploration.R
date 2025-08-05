library(tidyverse)
library(ggridges)
library(multcompView)
library(rstatix)
library(tibble)

# Load clean data, remove points with no egg weight or moisture data
canada_2023_df <- read_rds("data/processed/clean_2023_canada_data.rds") |> 
  drop_na(g_egg, pct_moisture)

## Histogram visualizations ----------------------------------------------------

# Visualize the distribution of egg weights at each site
p_egg_mass_distribution <- ggplot(data = canada_2023_df, 
                                  aes(x = g_egg, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  labs(
    x = "Average egg mass (g)",
    y = "Site",
    title = "Canda-bound Egg Mass Distribution"
  ) +
  theme_bw()

p_egg_mass_distribution 

#ggsave(
#  filename = "output/figures/egg_mass_distribution_canada.png", 
#  plot = p_egg_mass_distribution,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)

## Boxplot visualizations and ANOVAs -------------------------------------------

### Eggs get larger as fish move up the river ----------------------------------
#### Boxplot -------------------------------------------------------------------
p_egg_mass_site <- ggplot() +
  geom_violin(data = canada_2023_df, 
              aes(x = site, y = g_egg, fill = group)) +
  geom_jitter(data = canada_2023_df, 
              aes(x = site, y = g_egg),
              width = 0.2) +
  labs(
    x = "Site",
    y = "Average Egg Mass (g)",
    title = "Average Egg Mass (g)"
  ) +
  theme_bw()

p_egg_mass_site

#### ANOVA --------------------------------------------------------------------- 
egg_mass_aov <- aov(g_egg ~ site, data = canada_2023_df)
summary.aov(egg_mass_aov) # significant, differences between sites exist

#### Tukey ---------------------------------------------------------------------
egg_mass_Tukey <- TukeyHSD(egg_mass_aov) # WH > FOYU = RARA > PIST
egg_mass_Tukey

#### Generate comparison letters -----------------------------------------------
cld_egg_mass <- multcompLetters4(egg_mass_aov, egg_mass_Tukey)

# Create df of letters and positions
cld_egg_mass_df <- as.data.frame.list(cld_egg_mass$site) |> 
  rownames_to_column(var = "site")

# Get y-positions for labels
egg_mass_text_df <- canada_2023_df |> 
  group_by(site) |> 
  summarise(y_pos = max(g_egg) + 0.025) |> 
  left_join(cld_egg_mass_df)

# Add letters to plot
p_egg_mass_site <- p_egg_mass_site +
  geom_text(data = egg_mass_text_df, aes(x = site, y = y_pos, label = Letters))

ggsave(
  filename = "output/figures/egg_mass_site_canada.png", 
  plot = p_egg_mass_site,                         
  width = 7,                                        
  height = 5,                                      
  dpi = 300                                        
)

### Eggs thiamine concentration decreases as fish move up the river ------------
#### Boxplot -------------------------------------------------------------------
p_egg_thiamine_conc_site <- ggplot() +
  geom_violin(data = canada_2023_df, 
              aes(x = site, y = nmol_T_g, fill = group)) +
  geom_jitter(data = canada_2023_df, 
              aes(x = site, y = nmol_T_g),
              width = 0.2) +
  labs(
    x = "Site",
    y = "Egg Thiamine Concentration (nmol/g)",
    title = "Egg Thiamine Concentration"
  ) +
  theme_bw()

p_egg_thiamine_conc_site

#### ANOVA --------------------------------------------------------------------- 
egg_conc_aov <- aov(nmol_T_g ~ site, data = canada_2023_df)
summary.aov(egg_conc_aov) # significant, differences between sites exist

#### Tukey ---------------------------------------------------------------------
egg_conc_Tukey <- TukeyHSD(egg_conc_aov) # WH > FOYU = RARA > PIST
egg_conc_Tukey

#### Generate comparison letters -----------------------------------------------
cld_egg_conc <- multcompLetters4(egg_conc_aov, egg_conc_Tukey)

# Create df of letters and positions
cld_egg_conc_df <- as.data.frame.list(cld_egg_conc$site) |> 
  rownames_to_column(var = "site")

# Get y-positions for labels
egg_conc_text_df <- canada_2023_df |> 
  group_by(site) |> 
  summarise(y_pos = max(nmol_T_g) + 2) |> 
  left_join(cld_egg_conc_df)

# Add letters to plot
p_egg_thiamine_conc_site <- p_egg_thiamine_conc_site +
  geom_text(data = egg_conc_text_df, aes(x = site, y = y_pos, label = Letters))

p_egg_thiamine_conc_site

#ggsave(
#  filename = "output/figures/egg_thiamine_conc_site_canada.png", 
#  plot = p_egg_thiamine_conc_site,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)
