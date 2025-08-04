library(tidyverse)
library(ggridges)
library(ggpubr)
library(multcompView)
library(rstatix)
library(tibble)

# Load clean data, remove points with no egg weight or moisture data
clean_2023_df <- read_rds("data/processed/clean_2023_data.rds") |> 
  drop_na(g_egg, pct_moisture)

# Remove Chena and Salcha just because they are a bit weird?
removals <- c("Chena", "Salcha")

clean_2023_df_nochsal <- clean_2023_df |> 
  filter(!(site %in% removals))


## Histogram visualizations ----------------------------------------------------

# Visualize the distribution of egg weights at each site
p_egg_mass_distribution <- ggplot(data = clean_2023_df, 
                                  aes(x = g_egg, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
    ) +
  labs(
    x = "Average egg mass (g)",
    y = "Site",
    title = "Egg Mass Distribution"
    ) +
  theme_bw()

p_egg_mass_distribution 

#ggsave(
#  filename = "output/figures/egg_mass_distribution_site.png", 
#  plot = p_egg_mass_distribution,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)

# Visualize the distribution of egg weights at each site no Chena, no Salcha
p_egg_mass_distribution_nochalcha <- ggplot(data = clean_2023_df_nochsal, 
                                           aes(x = g_egg, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  labs(
    x = "Average egg mass (g)",
    y = "Site",
    title = "Egg Mass Distribution, Chena and Salcha Removed"
  ) +
  theme_bw()

p_egg_mass_distribution_nochalcha

ggsave(
  filename = "output/figures/egg_mass_distribution_site_nochalcha.png", 
  plot = p_egg_mass_distribution_nochalcha,                         
  width = 7,                                        
  height = 5,                                      
  dpi = 300                                        
)

# Visualize the distribution of egg thiamine grams at each site
ggplot(data = clean_2023_df, aes(x = nmol_T_egg, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  theme_bw()


# Visualize the distribution of egg thiamine grams at each site
ggplot(data = clean_2023_df_nochsal, aes(x = nmol_T_egg, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  theme_bw()

## Boxplot visualizations and ANOVAs -------------------------------------------

# site comparisons list for anova comparisons
site_comparisons <- list(
  c("WH", "FOYU"),
  c("FOYU", "RARA"),
  c("RARA", "PIST"),
  c("WH", "RARA"),
  c("WH", "PIST"),
  c("FOYU", "PIST")
  )

### Eggs get larger as fish move up the river ----------------------------------
#### Boxplot -------------------------------------------------------------------
p_egg_mass_site <- ggplot() +
  geom_violin(data = clean_2023_df_nochsal, 
              aes(x = site, y = g_egg, fill = group)) +
  geom_jitter(data = clean_2023_df_nochsal, 
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
egg_mass_aov <- aov(g_egg ~ site, data = clean_2023_df_nochsal)
summary.aov(egg_mass_aov) # significant, differences between sites exist

#### Tukey ---------------------------------------------------------------------
egg_mass_Tukey <- TukeyHSD(egg_mass_aov) # WH > FOYU = RARA > PIST

#### Generate comparison letters -----------------------------------------------
cld_egg_mass <- multcompLetters4(egg_mass_aov, egg_mass_Tukey)

# Create df of letters and positions
cld_egg_mass_df <- as.data.frame.list(cld_egg_mass$site) |> 
  rownames_to_column(var = "site")

# Get y-positions for labels
egg_mass_text_df <- clean_2023_df_nochsal |> 
  group_by(site) |> 
  summarise(y_pos = max(g_egg) + 0.075) |> 
  left_join(cld_egg_mass_df)

# Add letters to plot
p_egg_mass_site +
  geom_text(data = egg_mass_text_df, aes(x = site, y = y_pos, label = Letters))

### Eggs thiamine concentration decreasesas fish move up the river -------------
#### Boxplot -------------------------------------------------------------------
p_egg_thiamine_conc_site <- ggplot() +
  geom_violin(data = clean_2023_df_nochsal, 
              aes(x = site, y = nmol_T_g, fill = group)) +
  geom_jitter(data = clean_2023_df_nochsal, 
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
egg_thiamine_conc_aov <- aov(nmol_T_g ~ site, data = clean_2023_df_nochsal)
summary.aov(egg_thiamine_conc_aov) # significant, differences between sites exist

#### Tukey ---------------------------------------------------------------------
egg_thiamine_conc_Tukey <- TukeyHSD(egg_thiamine_conc_aov) # WH > FOYU = RARA > PIST
egg_thiamine_conc_Tukey

#### Generate comparison letters -----------------------------------------------
cld_egg_thiamine_conc <- multcompLetters4(egg_thiamine_conc_aov, egg_thiamine_conc_Tukey)

# Create df of letters and positions
cld_egg_thiamine_conc_df <- as.data.frame.list(cld_egg_thiamine_conc$site) |> 
  rownames_to_column(var = "site")

# Get y-positions for labels
egg_thiamine_conc_text_df <- clean_2023_df_nochsal |> 
  group_by(site) |> 
  summarise(y_pos = max(nmol_T_g) + 2) |> 
  left_join(cld_egg_thiamine_conc_df)

# Add letters to plot
p_egg_thiamine_conc_site +
  geom_text(data = egg_thiamine_conc_text_df, aes(x = site, y = y_pos, label = Letters))

### Total egg thiamine per egg does not change as fish move up the river -------
#### Boxplot -------------------------------------------------------------------
p_egg_thiamine_total_site <- ggplot() +
  geom_violin(data = clean_2023_df_nochsal, 
              aes(x = site, y = nmol_T_egg, fill = group)) +
  geom_jitter(data = clean_2023_df_nochsal, 
              aes(x = site, y = nmol_T_egg), 
              width = 0.2) +
  labs(
    x = "Site",
    y = "Total Egg Thiamine (nmol/egg)",
    title = "Total Egg Thiamine"
  ) +
  theme_bw()

p_egg_thiamine_total_site

#### ANOVA --------------------------------------------------------------------- 
egg_thiamine_total_aov <- aov(nmol_T_egg ~ site, data = clean_2023_df_nochsal)
summary.aov(egg_thiamine_total_aov) # no significant, differences between sites

#### Tukey ---------------------------------------------------------------------
egg_thiamine_total_Tukey <- TukeyHSD(egg_thiamine_total_aov) 
egg_thiamine_total_Tukey

#### Generate comparison letters -----------------------------------------------
cld_egg_thiamine_total <- multcompLetters4(egg_thiamine_total_aov, egg_thiamine_total_Tukey)

# Create df of letters and positions
cld_egg_thiamine_total_df <- as.data.frame.list(cld_egg_thiamine_total$site) |> 
  rownames_to_column(var = "site")

# Get y-positions for labels
egg_thiamine_total_text_df <- clean_2023_df_nochsal |> 
  group_by(site) |> 
  summarise(y_pos = max(nmol_T_egg) + 0.5) |> 
  left_join(cld_egg_thiamine_total_df)

# Add letters to plot
p_egg_thiamine_total_site <- p_egg_thiamine_total_site +
  geom_text(data = egg_thiamine_total_text_df, aes(x = site, y = y_pos, label = Letters))

p_egg_thiamine_total_site

#### Save plot -----------------------------------------------------------------
#ggsave(
#  filename = "output/figures/egg_thiamine_total_site.png", 
#  plot = p_egg_thiamine_total_site,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)
