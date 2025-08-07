library(tidyverse)
library(readxl)

# Load raw data
raw_2023_df <- read_rds("data/raw/yukon_egg_thiamine_2023.rds")

# Load PIST genetic assignment
pist_2023_genetics <- read_csv("data/raw/2023_yukon_PIST_genetic_assignment.csv") |> 
  mutate(collection_date = mdy(`Sample Date`))

# Load PIST field data for lengths/weights
pist_2023_lw <- read_csv("data/raw/2023_PILOT_FIELD_DATA.csv")

# Load 2023 Rapids data (Ich data for RARA and FOYU)
rapids_2023_ich <- read_csv("data/raw/2023_Rapids.csv") |> 
  mutate(collection_date = mdy(Date),
         time_of_death = ymd_hms(paste(collection_date, `Fish Death Time`)),
         blood_draw_time = ymd_hms(paste(collection_date, `Blood Drawing Time`)))

# Add additional variables to 2023 Thiamine data -------------------------------
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

# Clean PIST genetic assignment data -------------------------------------------
pist_2023_genetics_clean <- pist_2023_genetics |> 
  mutate(new_ID = paste0("PIST23_", sprintf("%03d", ID))) |> 
  select(new_ID, Country, `Broad Group`, `Fine Group`, collection_date) |> 
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

# Clean 2023 PIST length weight data -------------------------------------------
pist_2023_lw_clean <- pist_2023_lw |> 
  mutate(id = paste0("PIST23_", sprintf("%03d", ID))) |> 
  rename(length_mideye_fork_mm = `METF Length`,
         fish_mass_kg = `Weight(kg)`,
         girth_mm = `Girth(mm)`) |> 
  select(id, length_mideye_fork_mm, fish_mass_kg, girth_mm)

# Join to genetics data 

pist_2023_genetics_clean <- pist_2023_genetics_clean |> 
  left_join(pist_2023_lw_clean)

# 2023 Rapids (Ich) ------------------------------------------------------------
rapids_2023_ich_clean <- rapids_2023_ich |> 
  filter(!is.na(`Egg Vial 1 Num`)) |> 
  mutate(
    thiamine_ID = case_when(
      Site == "Rapids" ~ paste0("RARA23_", sprintf("%03d", `Egg Vial 1 Num`)),
      Site == "Fort Yukon" ~ paste0("FOYU23_", sprintf("%03d", `Egg Vial 1 Num`))
      )
    ) |> 
  select(
    thiamine_ID, collection_date, `Latitude (N)`, `Longitude (W)`, `Length Mideye Fork (mm)`,
    `Length Nose Fork (mm)`, `Ich Visual Positive`, `Ich Visual Count`
  ) |> 
  mutate(
    ich_visual_positive = case_when(
      `Ich Visual Positive` == "No" ~ FALSE,
      `Ich Visual Positive` == "Yes" ~ TRUE,
      .default = NA
      )
  ) |> 
  select(
    -`Ich Visual Positive`
  ) |> 
  rename(
    lat_n = `Latitude (N)`,
    long_w = `Longitude (W)`,
    length_mideye_fork_mm = `Length Mideye Fork (mm)`,
    length_nose_fork_mm = `Length Nose Fork (mm)`,
    ich_visual_count = `Ich Visual Count`
  )

# Join PIST genetic assignment data to whole dataset ---------------------------
clean_2023_df_join1 <- clean_2023_df |> 
  left_join(pist_2023_genetics_clean, join_by(SIN == id))

# Add fine_group designations for other reaches of the river
# site == "WH", fine_group == "Canada"
# site == "FOYU" | site == "RARA", fine_group == "Canada"
# site == "Chena" | site == "Salcha", fine_group == "Tanana"
clean_2023_df_join1 <- clean_2023_df_join1 |>
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

# Join Rapids 2023 data to whole dataset ---------------------------------------
clean_2023_df_join2 <- clean_2023_df_join1 |> 
  left_join(rapids_2023_ich_clean, join_by(SIN == thiamine_ID)) |> 
  mutate(
    collection_date = coalesce(collection_date.x, collection_date.y),
    length_mideye_fork_mm = coalesce(length_mideye_fork_mm.x, length_mideye_fork_mm.y)
    ) |> 
  select(-collection_date.x, -collection_date.y,
         -length_mideye_fork_mm.x, -length_mideye_fork_mm.y)

# Write clean data -------------------------------------------------------------
write_rds(clean_2023_df_join2, "data/processed/clean_2023_data.rds")

write_csv(clean_2023_df_join2, "data/processed/clean_2023_yukon_data.csv")

# Write a clean copy of only Canada fish ---------------------------------------
clean_2023_canada_df <- clean_2023_df_join2 |> 
  filter(fine_group == "Canada")

write_rds(clean_2023_canada_df, "data/processed/clean_2023_canada_data.rds")
