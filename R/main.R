if (!require("pacman")) install.packages("pacman")

pkgs =
  c(
    "here",
    "lubridate",
    "sf",
    "taxize",
    "tidyverse"
  )

pacman::p_load(pkgs, character.only = T)

dir.create(here("data"))

# Download rodent trapping data
combined_data <- readRDS(gzcon(url("https://github.com/DidDrog11/rodent_trapping/raw/main/data/data_for_export/combined_data.rds"))) %>%
  write_rds(here("data", "combined_data.rds"))

# Download ELISA results
ELISA_data <- readRDS(gzcon(url("https://github.com/DidDrog11/SL_lassa_ELISA/raw/main/output/ELISA_output.rds"))) %>%
  write_rds(here("data", "ELISA_data.rds"))

# Match species names to NCBI taxa ID
ncbi_id <- tibble(species_name =  unique(str_to_sentence(str_replace_all(combined_data$rodent_data$clean_names, "_", " ")))) %>%
  mutate(ncbi_tax_id = unname(c(get_ids(species_name, db = "ncbi")$ncbi)))

# Match coordinates to trap data and extract latitude and longitude
rodent_coords <- tibble(rodent_uid = combined_data$rodent_data$rodent_uid) %>%
  left_join(combined_data$trap_data %>%
              select(rodent_uid) %>%
              mutate(Latitude = st_coordinates(geometry)[, 1],
                     Longitude = st_coordinates(geometry)[, 2])) %>%
  select(-geometry)

# Date of collection
rodent_dates <- tibble(rodent_uid = combined_data$rodent_data$rodent_uid) %>%
  left_join(combined_data$trap_data %>%
              tibble() %>%
              select(rodent_uid, date_set)) %>%
  mutate(`Collection day` = day(date_set),
         `Collection month` = month(date_set),
         `Collection year` = year(date_set)) %>%
  select(-date_set)

# Lassa mammarenavirus NCBI taxa id
LASV <- unname(c(get_ids("Lassa mammarenavirus", db = "ncbi")$ncbi))

# Sample detection results and processing
elisa_results <- ELISA_data$ELISA_enriched %>%
  select(rodent_uid, interpretation, result) %>%
  mutate(`Collection method or tissue` = "blood",
         `Detection method` = case_when(!is.na(interpretation) ~ "ELISA - antibody",
                                        TRUE ~ NA),
         `Detection target` = case_when(!is.na(interpretation) ~ "Lassa mammarenavirus",
                                        TRUE ~ NA),
         `Detection target NCBI tax ID` = case_when(!is.na(interpretation) ~ LASV,
                                                    TRUE ~ NA),
         `Detection outcome` = str_to_lower(interpretation),
         `Detection measurement` = result,
         `Detection measurement units` = case_when(!is.na(interpretation) ~ "Optical density 450 - 630",
                                                   TRUE ~ NA),
         `Pathogen` = case_when(!is.na(interpretation) ~ "Lassa mammarenavirus",
                                TRUE ~ NA),
         `Pathogen NCBI tax id` = case_when(!is.na(interpretation) ~ LASV,
                                TRUE ~ NA)) %>%
  select(-interpretation, -result)

# Find the names of the columns required
pharos_template <- read_csv(here("data", "pharos_template.csv"))

# Produce the final dataset
pharos_data <- tibble(
  `Sample ID` =  combined_data$rodent_data$rodent_uid,
  `Organism ID` = combined_data$rodent_data$rodent_uid,
  `Host species` = str_to_sentence(str_replace_all(combined_data$rodent_data$clean_names, "_", " "))) %>%
  left_join(ncbi_id, by = c("Host species" = "species_name")) %>%
  rename(`Host species NCBI tax ID` = ncbi_tax_id) %>%
  left_join(rodent_coords, by = c("Organism ID" = "rodent_uid")) %>%
  mutate(`Spatial uncertainty` = 10) %>%
  left_join(rodent_dates, by = c("Organism ID" = "rodent_uid")) %>%
  left_join(elisa_results, by = c("Organism ID" = "rodent_uid")) %>%
  mutate(`Organism sex` = combined_data$rodent_data$sex,
         `Dead or alive` = "alive",
         `Health notes` = "Length refers to head-body length of trapped rodents and shrews",
         # convert not known lifestage to unknown
         `Life stage` = combined_data$rodent_data$age_group,
         `Life stage` = case_when(`Life stage` == "not_known" ~ "unknown",
                                  TRUE ~ `Life stage`),
         `Age` = NA,
         # convert mass from grams to kg
         `Mass` = combined_data$rodent_data$weight/1000,
         # convert length from mm to m
         `Length` = combined_data$rodent_data$head_body/1000)

write_csv(pharos_data, here("data", "rodent_eastern_sl_2023-06-22.csv"))
