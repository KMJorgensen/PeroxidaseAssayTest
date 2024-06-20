# Cleaning data: Test Mn-dependency 

# Load packages ####
library(tidyverse) # ver. 2.0.0
library(readxl)

# Load data ####

# Dataset was produced 2023-07-28
# The Excel sheet contains lots of extra information about the read-method. Only read rows with absorbance data (rows 44:104) 
H2O2.raw <- read_excel("Data/Raw/MnP_enzyme_assay_20230728_121526.xlsx", range = "A44:CU104")
#Fix column names
colnames(H2O2.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))

# Load samples and treatments
H2O2.samples <- read_excel("Data/Raw/Template_samples_20230728_121526.xlsx", range = "A1:I9")
H2O2.treats <- read_excel("Data/Raw/Template_samples_20230728_121526.xlsx", range = "A12:I20")

# Wrangle data ####

# Wrangle sample and treatment data into long format with "well_ID" as ID-marker

H2O2.samples <-  H2O2.samples |> 
  pivot_longer(cols = 2:9, names_to = "column", values_to = "MnP_conc") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

H2O2.treats <- H2O2.treats |> 
  pivot_longer(cols = 2:9, names_to = "column", values_to = "H2O2_conc") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

H2O2.samples <- full_join(H2O2.samples, H2O2.treats, by = "well_ID")
# Dataframe ready to be merged with absorbance data

# Wrangle dataset with absorbance reads into long format
# Round Time [s] to 0 digits
H2O2.raw <- H2O2.raw |> 
  mutate(time_s = floor(time_s))

H2O2.raw.long <- H2O2.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")

H2O2.clean <- full_join(H2O2.raw.long, H2O2.samples, by = "well_ID") 

H2O2.clean <- H2O2.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = H2O2_conc) |> 
  select(-c(cycle_nr, temp, time_s))  #drops uninformative columns

# Export cleaned dataframe
write.csv(H2O2.clean, "Data/Clean/clean_H2O2_concentration.csv", row.names = FALSE)


