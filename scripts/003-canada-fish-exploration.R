library(tidyverse)
library(ggridges)
library(multcompView)
library(rstatix)
library(tibble)
library(car)
library(rcompanion)

# Set ggplot theme to bw
theme_set(theme_bw())

# Load clean data, remove points with no egg weight or moisture data
canada_2023_df <- read_rds("data/processed/clean_2023_canada_data.rds") |> 
  drop_na(nmol_T_g)

canada_2023_df_egg_mass_trim <- canada_2023_df |> 
  drop_na(g_egg)

canada_2023_df_moisture_trim <- canada_2023_df |> 
  drop_na(pct_moisture)

## Visualizations and Difference Tests -----------------------------------------

### Average Egg Mass -----------------------------------------------------------

#### Distribution --------------------------------------------------------------

p_egg_mass_distribution <- ggplot(data = canada_2023_df_egg_mass_trim, 
                                  aes(x = g_egg, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  labs(
    x = "Average egg mass (g)",
    y = "Site"
  ) 

p_egg_mass_distribution 

#ggsave(
#  filename = "output/figures/egg_mass_distribution_canada.png", 
#  plot = p_egg_mass_distribution,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)

#### Boxplot -------------------------------------------------------------------
p_egg_mass_site <- ggplot() +
  geom_violin(data = canada_2023_df_egg_mass_trim, 
              aes(x = site, y = g_egg, fill = group)) +
  geom_jitter(data = canada_2023_df_egg_mass_trim, 
              aes(x = site, y = g_egg),
              width = 0.2) +
  labs(
    x = "Site",
    y = "Average Egg Mass (g)"
  )

p_egg_mass_site

#### Homoscedasticity of Variance Test -----------------------------------------
leveneTest(g_egg ~ site, data = canada_2023_df_egg_mass_trim)

# Variance is homoscedastic, ANOVA is valid

#### ANOVA --------------------------------------------------------------------- 
egg_mass_aov <- aov(g_egg ~ site, data = canada_2023_df_egg_mass_trim)
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
egg_mass_text_df <- canada_2023_df_egg_mass_trim |> 
  group_by(site) |> 
  summarise(y_pos = max(g_egg) + 0.025) |> 
  left_join(cld_egg_mass_df)

# Add letters to plot
p_egg_mass_site <- p_egg_mass_site +
  geom_text(data = egg_mass_text_df, aes(x = site, y = y_pos, label = Letters))

p_egg_mass_site

#ggsave(
#  filename = "output/figures/egg_mass_site_canada.png", 
#  plot = p_egg_mass_site,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)

### Egg Thiamine Concentration -------------------------------------------------

#### Distribution --------------------------------------------------------------

p_egg_thiamine_conc_distribution <- ggplot(data = canada_2023_df_egg_mass_trim, 
                                           aes(x = nmol_T_g, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  labs(
    x = "Egg Thiamine Concentration (nmol/g)",
    y = "Site"
  ) 

p_egg_thiamine_conc_distribution 

#ggsave(
#  filename = "output/figures/egg_thiamine_concentration_distribution_canada.png", 
#  plot = p_egg_thiamine_conc_distribution,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)

#### Boxplot -------------------------------------------------------------------

p_egg_thiamine_conc_site <- ggplot() +
  geom_violin(data = canada_2023_df, 
              aes(x = site, y = nmol_T_g, fill = group)) +
  geom_jitter(data = canada_2023_df, 
              aes(x = site, y = nmol_T_g),
              width = 0.2) +
  labs(
    x = "Site",
    y = "Egg Thiamine Concentration (nmol/g)"
  )

p_egg_thiamine_conc_site

#### Homoscedasticity of Variance Test -----------------------------------------
leveneTest(nmol_T_g ~ site, data = canada_2023_df)

# The data fail the homoscedasticity of variance test and they are not 
# normally distributed. I need to use a non-parametric test (Kruskal-Wallis)

#### Kruskal-Wallis ------------------------------------------------------------ 
egg_conc_kw <- kruskal.test(nmol_T_g ~ site, data = canada_2023_df)
egg_conc_kw

# significant differences between sites exist
# Since not an ANOVA, can't run a Tukey test, Dunn Test is the alternative for
# Kruskal-Wallis

#### Dunn Test -----------------------------------------------------------------
egg_conc_dunn <- dunn_test(data = canada_2023_df,
                            formula = nmol_T_g ~ site,
                            p.adjust.method = "bonferroni") # WH > FOYU = RARA > PIST

egg_conc_dunn_formatted <- egg_conc_dunn %>%
  mutate(comparison = paste(group1, group2, sep = "-"))

egg_conc_dunn_formatted

#### Generate comparison letters -----------------------------------------------
cld_egg_conc <- cldList(p.adj ~ comparison,
                        data = egg_conc_dunn_formatted,
                        threshold = 0.05) 

# Create df of letters and positions
cld_egg_conc_df <- cld_egg_conc |> 
  rename(site = Group)

# Get y-positions for labels
egg_conc_text_df <- canada_2023_df |> 
  group_by(site) |> 
  summarise(y_pos = max(nmol_T_g) + 2) |> 
  left_join(cld_egg_conc_df)

# Add letters to plot
p_egg_thiamine_conc_site <- p_egg_thiamine_conc_site +
  geom_text(data = egg_conc_text_df, aes(x = site, y = y_pos, label = Letter))

p_egg_thiamine_conc_site

#ggsave(
#  filename = "output/figures/egg_thiamine_conc_site_canada.png", 
#  plot = p_egg_thiamine_conc_site,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)

### Egg percent moisture -------------------------------------------------------

#### Distribution --------------------------------------------------------------

p_egg_pct_moisture_distribution <- ggplot(data = canada_2023_df_moisture_trim, 
                                           aes(x = pct_moisture, y = site, fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = "raincloud",
    alpha = 0.75
  ) +
  labs(
    x = "Egg % Moisture",
    y = "Site"
  ) 

p_egg_pct_moisture_distribution 

# Not a normal distribution, can't use ANOVA or Welch's ANOVA

#ggsave(
#  filename = "output/figures/egg_pct_moisture_distribution_canada.png", 
#  plot = p_egg_pct_moisture_distribution,                         
#  width = 7,                                        
#  height = 5,                                      
#  dpi = 300                                        
#)

#### Boxplot -------------------------------------------------------------------

p_egg_pct_moisture_site <- ggplot() +
  geom_violin(data = canada_2023_df_moisture_trim, 
              aes(x = site, y = pct_moisture, fill = group)) +
  geom_jitter(data = canada_2023_df_moisture_trim, 
              aes(x = site, y = pct_moisture),
              width = 0.2) +
  labs(
    x = "Site",
    y = "Egg % Moisture"
  )

p_egg_pct_moisture_site

#### Homoscedasticity of Variance Test -----------------------------------------
leveneTest(pct_moisture ~ site, data = canada_2023_df_moisture_trim)

# Passes the Levene Test (p > 0.05), so the variance is homoscedastic between groups

#### Kruskal-Wallis ------------------------------------------------------------ 
egg_pct_moisture_kw <- kruskal.test(pct_moisture ~ site,
                                    data = canada_2023_df_moisture_trim)
egg_pct_moisture_kw 

# significant differences between sites exist

#### Dunn's Test ---------------------------------------------------------------
egg_pct_moisture_dunn <- dunn_test(data = canada_2023_df_moisture_trim,
                           formula = pct_moisture ~ site,
                           p.adjust.method = "bonferroni") # WH > FOYU = RARA > PIST

egg_pct_moisture_dunn_formatted <- egg_pct_moisture_dunn %>%
  mutate(comparison = paste(group1, group2, sep = "-"))

egg_pct_moisture_dunn_formatted

#### Generate comparison letters -----------------------------------------------
cld_egg_pct_moisture <- cldList(p.adj ~ comparison,
                        data = egg_pct_moisture_dunn_formatted,
                        threshold = 0.05) 

# Create df of letters and positions
cld_egg_pct_moisture_df <- cld_egg_pct_moisture |> 
  rename(site = Group)

# Get y-positions for labels
egg_pct_moisture_text_df <- canada_2023_df_moisture_trim |> 
  group_by(site) |> 
  summarise(y_pos = max(pct_moisture) + 0.025) |> 
  left_join(cld_egg_pct_moisture_df)

# Add letters to plot
p_egg_pct_moisture_site <- p_egg_pct_moisture_site +
  geom_text(data = egg_pct_moisture_text_df, aes(x = site, y = y_pos, label = Letter))

p_egg_pct_moisture_site

ggsave(
  filename = "output/figures/egg_pct_moisture_site_canada.png", 
  plot = p_egg_pct_moisture_site,                         
  width = 7,                                        
  height = 5,                                      
  dpi = 300                                        
)
