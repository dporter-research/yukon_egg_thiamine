# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Analysis of Yukon River Chinook Salmon Egg Composition
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#==============================================================================
# 1. SETUP: Load libraries and data
#================================C==============================================

# Load all necessary packages at the start
library(tidyverse) # For data manipulation and plotting (includes ggplot2, dplyr)
library(here)      # For easy file path management
library(mgcv)      # For fitting General Additive Models (GAMs)


# Load the raw data from an .rds file
yukon_data <- read_rds(here("data", "yukon_egg_thiamine_2023.rds"))


#==============================================================================
# 2. DATA CLEANING AND PREPARATION
#==============================================================================

# Remove rows with missing values in key measurement columns
yukon_clean <- yukon_data |>
  drop_na(nmol_T_g, g_egg, pct_moisture)

# --- Calculate estimated protein percentage ---
# This is done "by difference," assuming dry matter is composed of lipids and protein.
# This is a strong hypothesis but should be confirmed with direct nitrogen analysis.
yukon_clean <- yukon_clean |>
  mutate(
    # Step 1: Calculate the percentage of the egg that is not water
    pct_dry_matter = 100 - pct_moisture,
    
    # Step 2: Convert lipid percentage from a wet-weight to a dry-weight basis
    pct_lipid_dry = (pct_lipid_wet / pct_dry_matter) * 100,
    
    # Step 3: Estimate protein as the remaining fraction of the dry matter
    pct_protein_est_dry = 100 - pct_lipid_dry,
    
    # --- Calculate total mass of each component ---
    lipid_g = g_egg * (pct_lipid_wet / 100),
    water_g = g_egg * (pct_moisture / 100),
    protein_g_est = g_egg * (pct_dry_matter / 100) * (pct_protein_est_dry / 100)
  )


#==============================================================================
# 3. EXPLORATORY ANALYSIS: How does egg composition change with size?
#==============================================================================

# --- Plot 1: Egg Weight vs. Thiamine Concentration ---
# Investigates the "dilution effect" of thiamine in larger eggs.
ggplot(yukon_clean, aes(x = g_egg, y = nmol_T_g)) +
  geom_point(aes(color = site), alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black") +
  labs(
    title = "Egg Size vs. Thiamine Concentration",
    x = "Egg Weight (g)",
    y = "Thiamine Concentration (nmol/g)"
  ) +
  theme_bw()

# --- Statistical Test 1: GAM for Egg Size vs. Thiamine ---
# Formally models the non-linear relationship seen in the plot.
gam_model <- gam(nmol_T_g ~ s(g_egg), data = yukon_clean)
summary(gam_model) # A significant p-value and edf > 1 confirms a non-linear relationship.

# INTERPRETATION OF RESULTS:
#
# p-value < 2e-16: The relationship is highly significant. Egg weight is a strong
#                  predictor of thiamine concentration.
#
# edf = 2.338: The 'Effective Degrees of Freedom' is greater than 1, confirming
#              the relationship is NON-LINEAR. A simple straight line would be
#              inappropriate.
#
# R-sq.(adj) = 0.466: The model explains 46.6% of the variation in thiamine
#                     concentration, indicating a moderately strong relationship.

# --- Plot 2: Egg Weight vs. Percent Lipid ---
# Tests the hypothesis that larger eggs have a different lipid concentration.
ggplot(yukon_clean, aes(x = g_egg, y = pct_lipid_wet)) +
  geom_point(aes(color = site), alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black") +
  labs(
    title = "Egg Size vs. Lipid Content",
    x = "Egg Weight (g)",
    y = "Percent Lipid (Wet Weight)"
  ) +
  theme_bw()

# Fit the GAM to test the relationship
gam_lipid <- gam(pct_lipid_wet ~ s(g_egg), data = yukon_clean)

# Get the summary to interpret the results
summary(gam_lipid)
# INTERPRETATION OF RESULTS:
#
# p-value < 2e-16: This is a highly significant result. It confirms that
#                  egg weight is a strong predictor of the wet lipid percentage.
#
# edf = 3.466: The 'Effective Degrees of Freedom' is well above 1, providing
#              strong statistical evidence that the relationship is NON-LINEAR.
#              A simple straight line would not have accurately captured the trend.
#
# R-sq.(adj) = 0.45: The model explains 45% of the variation in wet lipid
#                    percentage. This is a moderately strong relationship,
#                    indicating egg size is a meaningful biological factor.

# --- Plot 3: Egg Weight vs. Estimated Protein ---
# Tests the "protein loading" hypothesis.
ggplot(yukon_clean, aes(x = g_egg, y = pct_protein_est_dry)) +
  geom_point(aes(color = site), alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black") +
  labs(
    title = "Egg Size vs. Estimated Protein Content",
    x = "Egg Weight (g)",
    y = "Estimated Protein % (Dry Weight)"
  ) +
  theme_bw()


#==============================================================================
# 4. HYPOTHESIS TESTING: Are there differences across river groups?
#==============================================================================

# --- Boxplot 1: Total Thiamine per Egg ---
# Visualizes the core hypothesis: is thiamine content stable across groups?
ggplot(data = yukon_clean, aes(x = factor(group), y = nmol_T_egg)) +
  geom_boxplot(aes(fill = factor(group)), show.legend = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    title = "Total Thiamine per Egg Remains Stable Across Groups",
    x = "Migration Group",
    y = "Total Thiamine (nmol/egg)"
  ) +
  theme_minimal()


# --- Boxplot 2: Egg Size (g) ---
# Visualizes if egg size changes across the groups.
ggplot(data = yukon_clean, aes(x = factor(group), y = g_egg)) +
  geom_boxplot(aes(fill = factor(group)), show.legend = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    title = "Egg Size Increases Across Groups",
    x = "Migration Group",
    y = "Egg Weight (g)"
  ) +
  theme_minimal()


#==============================================================================
# 5. STATISTICAL ANALYSIS: ANOVA for River Group Differences
#==============================================================================
# This section formally tests the patterns observed in the boxplots above.

# Ensure the 'group' variable is treated as a factor for modeling
yukon_clean$group <- as.factor(yukon_clean$group)

# --- ANOVA Test 1: Total Thiamine per Egg (nmol_T_egg) ---
# Tests the hypothesis that thiamine per egg does NOT change.
# A non-significant result (p > 0.05) is the expected outcome.
anova_thiamine_egg <- aov(nmol_T_egg ~ group, data = yukon_clean)
summary(anova_thiamine_egg)
# INTERPRETATION OF RESULTS:
#
# p-value = 0.271: The p-value is much greater than 0.05. This is NOT a
#                  significant result. We cannot reject the null hypothesis.
#
# CONCLUSION: This supports the core hypothesis. There is no statistical
#             difference in the total amount of thiamine per egg between the
#             "Entry", "Transit", and "Spawning" groups.


# --- ANOVA Test 2: Egg Weight (g_egg) ---
# Tests if egg weight changes significantly across groups.
anova_egg_weight <- aov(g_egg ~ group, data = yukon_clean)
summary(anova_egg_weight)
# INTERPRETATION OF RESULTS:
#
# p-value < 2e-16: The p-value is extremely small, indicating a HIGHLY
#                  significant result. Egg weight is significantly
#                  different among the three groups.


# If the ANOVA is significant, run Tukey's HSD to see which specific groups differ.
TukeyHSD(anova_egg_weight)
# INTERPRETATION OF TUKEY HSD RESULTS:
#
# All "p adj" values are essentially zero. This means that the mean egg weight
# is significantly different between ALL pairs of groups:
#   - "Transit" eggs are heavier than "Entry" eggs.
#   - "Spawning" eggs are heavier than "Entry" eggs.
#   - "Spawning" eggs are heavier than "Transit" eggs.


# --- ANOVA Test 3: Estimated Protein Mass (protein_g_est) ---
# Tests if the total mass of protein changes across groups.
anova_protein_mass <- aov(protein_g_est ~ group, data = yukon_clean)
summary(anova_protein_mass)
# INTERPRETATION OF RESULTS:
#
# p-value < 2e-16: This is a highly significant result. The estimated
#                  total mass of protein is significantly different among the
#                  "Entry", "Transit", and "Spawning" groups.

# Check which specific groups differ.
TukeyHSD(anova_protein_mass)
# INTERPRETATION OF TUKEY HSD RESULTS:
#
# All "p adj" values are zero. This confirms that the estimated protein mass
# is significantly different between ALL pairs of groups, increasing at each
# stage of the migration.


# --- ANOVA Test 4: Lipid Mass (lipid_g) ---
# Tests if the total mass of lipids changes across groups.
anova_lipid_mass <- aov(lipid_g ~ group, data = yukon_clean)
summary(anova_lipid_mass)
# INTERPRETATION OF RESULTS:
#
# p-value = 3.27e-11: This is a highly significant result. The total mass
#                     of lipids is also significantly different among the
#                     three groups.

# Check which specific groups differ.
TukeyHSD(anova_lipid_mass)
# INTERPRETATION OF TUKEY HSD RESULTS:
#
# All "p adj" values are less than 0.05. This confirms that the total lipid
# mass is also significantly different between ALL pairs of groups,
# increasing at each stage of the migration.
