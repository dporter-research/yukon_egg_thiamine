library(tidyverse)
library(here)
library(ggpmisc)
library(mgcv)
library(psych)

yukon_data <- read_rds(here("data", "yukon_egg_thiamine_2023.rds"))

yukon_clean <- yukon_data |> 
  drop_na(nmol_T_g, g_egg, pct_moisture)

# Egg weight vs pct moisture
ggplot(yukon_clean, aes(x = g_egg, y = pct_moisture)) +
  geom_point(aes(color = site))


# Egg weight vs thiamine conc
ggplot(yukon_clean, aes(x = g_egg, y = nmol_T_g)) +
  geom_point(aes(color = site)) +
  geom_smooth(method = "gam",
              formula = y ~ s(x)) 
  #stat_poly_eq(
  #  aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
  #  formula = y ~ x,
  #  parse = TRUE
  #)


# Perform Spearman's rank correlation test
cor.test(yukon_clean$g_egg, yukon_clean$nmol_T_g, method = "spearman")

gam_model <- gam(nmol_T_g ~ s(g_egg), data = yukon_clean)

summary(gam_model)




library(ggplot2)

# Plot egg weight vs. percent moisture
ggplot(yukon_clean, aes(x = g_egg, y = pct_moisture)) +
  geom_point(aes(color = site), size = 4, alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  labs(
    title = "Testing the Dilution Mechanism",
    subtitle = "Larger eggs appear to have higher moisture content",
    x = "Egg Weight (g)",
    y = "Percent Moisture"
  ) +
  theme_bw()



# Plot egg weight vs. DRY WEIGHT thiamine concentration
ggplot(yukon_clean, aes(x = g_egg, y = nmol_T_g_dw_corrected)) +
  geom_point(aes(color = site), alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  labs(
    title = "Confirming the Dilution Hypothesis",
    subtitle = "The trend disappears after correcting for moisture",
    x = "Egg Weight (g)",
    y = "Thiamine Concentration (nmol/g, Dry Weight)"
  ) +
  theme_bw()

library(ggplot2)

# Plot egg weight vs. percent lipid
ggplot(yukon_clean, aes(x = g_egg, y = pct_lipid_wet)) +
  geom_point(aes(color = site), alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x)) + # Using a GAM to fit a flexible curve
  labs(
    title = "Investigating the Role of Lipids",
    subtitle = "Are larger eggs 'fattier'?",
    x = "Egg Weight (g)",
    y = "Percent Lipid (Wet Weight)"
  ) +
  theme_bw()

yukon_clean$dry_matter <- 100 - yukon_clean$pct_moisture
yukon_clean$pct_lipid_dry <- (yukon_clean$pct_lipid_wet / yukon_clean$dry_matter) * 100
yukon_clean$pct_protein_est <- 100 - yukon_clean$pct_lipid_dry


library(ggplot2)

# Plot egg weight vs. your new estimated protein column
ggplot(yukon_clean, aes(x = g_egg, y = pct_protein_est)) +
  geom_point(aes(color = site), alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  labs(
    title = "Testing the Protein Loading Hypothesis",
    x = "Egg Weight (g)",
    y = "Estimated Protein % (Dry Weight)"
  ) +
  theme_bw()


ggplot(yukon_clean, aes(x = pct_protein_est, y = nmol_T_g)) +
  geom_point(aes(color = site), alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  labs(
    x = "Estimated Protein % (Dry Weight)",
    y = "nmol_T_g"
  ) +
  theme_bw()


ggplot(data = yukon_clean, aes(x = factor(group), y = pct_protein_est)) +
  geom_boxplot()


yukon_clean$lipid_g <- yukon_clean$g_egg * (yukon_clean$pct_lipid_dry / 100)
yukon_clean$protein_g <- yukon_clean$g_egg * (yukon_clean$pct_protein_est / 100)
yukon_clean$water_g <- yukon_clean$g_egg * (yukon_clean$pct_moisture / 100)

ggplot(data = yukon_clean, aes(x = factor(group), y = lipid_g)) +
  geom_boxplot(aes(fill = factor(group)))

ggplot(data = yukon_clean, aes(x = factor(group), y = protein_g)) +
  geom_boxplot(aes(fill = factor(group)))

ggplot(data = yukon_clean, aes(x = factor(group), y = water_g)) +
  geom_boxplot(aes(fill = factor(group)))

ggplot(data = yukon_clean, aes(x = factor(group), y = g_egg)) +
  geom_boxplot(aes(fill = factor(group)), outlier.shape = NA) +
  geom_jitter(aes(color = survivor), width = 0.2)

ggplot(data = yukon_clean, aes(x = factor(group), y = nmol_T_egg)) +
  geom_boxplot(aes(fill = factor(group))) +
  geom_jitter(width = 0.2)


describeBy(yukon_clean$nmol_T_egg, group = yukon_clean$group)



yukon_clean_fallout <- yukon_clean |> 
  filter(nmol_T_egg < 0.63)

describeBy(yukon_clean_fallout$nmol_T_egg, group = yukon_clean_fallout$group)

8/46

5/36


dead_sins <- yukon_clean_fallout$SIN


yukon_clean <- yukon_clean |> 
  mutate(survivor = if_else(SIN %in% dead_sins, FALSE, TRUE))


ggplot(data = yukon_clean, aes(x = factor(group), y = nmol_T_egg)) +
  geom_boxplot(aes(fill = factor(group))) +
  geom_jitter(aes(color = survivor), width = 0.2)


