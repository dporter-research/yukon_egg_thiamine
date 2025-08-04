library(tidyverse)


# Load raw data
raw_2023_df <- read_rds("data/raw/yukon_egg_thiamine_2023.rds")

# Calculate estimated protein and total mass of each component
# Note: Protein is estimated by difference, a hypothesis to be confirmed.
clean_2023_df <- raw_2023_df |>
  mutate(
    pct_moisture = pct_moisture / 100,
    pct_lipid_wet = pct_lipid_wet / 100
  ) |> 
  mutate(
    pct_dry_matter = 1 - pct_moisture, # % of egg that is not water
    pct_lipid_dry = (pct_lipid_wet / pct_dry_matter), 
    pct_protein_est_dry = 1 - pct_lipid_dry,
    lipid_g = g_egg * pct_lipid_wet,
    water_g = g_egg * pct_moisture,
    protein_g_est = g_egg * pct_dry_matter * pct_protein_est_dry
  )

# Write clean data
write_rds(clean_2023_df, "data/processed/clean_2023_data.rds")
