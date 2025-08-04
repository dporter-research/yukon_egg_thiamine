library(tidyverse)


# Load raw data
raw_2023_df <- read_rds("data/raw/yukon_egg_thiamine_2023.rds")

# Load PIST genetic assignment
pist_2023_genetics <- read_csv("data/raw/2023_yukon_PIST_genetic_assignment.csv")

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

# Clean PIST genetic assignment data
pist_2023_genetics <- pist_2023_genetics |> 
  mutate(new_ID = paste0("PIST23_", sprintf("%03d", ID))) |> 
  select(new_ID, Country, `Broad Group`, `Fine Group`) |> 
  rename(
    id = new_ID,
    country = Country,
    broad_group = `Broad Group`,
    fine_group = `Fine Group`
  ) |> 
  mutate(
    country = factor(country),
    broad_group = factor(broad_group),
    fine_group = factor(fine_group)
  )

# Join PIST genetic assignment data to whole dataset
clean_2023_df <- clean_2023_df |> 
  left_join(pist_2023_genetics, join_by(SIN == id))

# Add fine_group designations for other reaches of the river
# site == "WH", fine_group == "Canada"
# site == "FOYU" | site == "RARA", fine_group == "Canada"
# site == "Chena" | site == "Salcha", fine_group == "Tanana"
clean_2023_df <- clean_2023_df |>
  mutate(
    fine_group = case_when(
      # Combine the "Canada" rules
      site %in% c("WH", "FOYU", "RARA") ~ "Canada",
      
      # The "Tanana" rule
      site %in% c("Chena", "Salcha") ~ "Tanana",
      
      # For all other rows, keep the existing value
      TRUE ~ fine_group
    )
  )

# Write clean data
write_rds(clean_2023_df, "data/processed/clean_2023_data.rds")

# Write a clean copy of only Canada fish
clean_2023_canada_df <- clean_2023_df |> 
  filter(fine_group == "Canada")

write_rds(clean_2023_canada_df, "data/processed/clean_2023_canada_data.rds")
