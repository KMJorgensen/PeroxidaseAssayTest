# Cleaning data: Masceration gradient, needle cultures


# Load packages ####
library(tidyverse) # ver. 2.0.0
library(readxl)

# Load data ####

# Dataset was produced 2023-11-28
# Analyses were made on eight plates, meaning that datasets will be cleaned separately but merged into one final clean dataset.


# Load data from all different substrate assays

MBTH_1.raw <- read_excel("Data/Raw/MnP_enzyme_assay_20231128_112423.xlsx", range = "A44:CU104")
MBTH_2.raw <- read_excel("Data/Raw/MnP_enzyme_assay_20231128_115603.xlsx", range = "A44:CU104")
ABTS_1.raw <- read_excel("Data/Raw/ABTS_peroxidase_assay_20231128_134719.xlsx", range = "A38:CU98")
ABTS_2.raw <- read_excel("Data/Raw/ABTS_peroxidase_assay_20231128_141906.xlsx", range = "A38:CU98")
LDOPA_1.raw <- read_excel("Data/Raw/L_DOPA_peroxidase_assay_20231128_144948.xlsx", range = "A38:CU98")
LDOPA_2.raw <- read_excel("Data/Raw/L_DOPA_peroxidase_assay_20231128_152219.xlsx", range = "A38:CU98")
DMP_1.raw <- read_excel("Data/Raw/DMP_peroxidase_assay_20231128_122728.xlsx", range = "A38:CU98")
DMP_2.raw <- read_excel("Data/Raw/DMP_peroxidase_assay_20231128_131524.xlsx", range = "A38:CU98")

#Fix column names
colnames(MBTH_1.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))
colnames(MBTH_2.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))
colnames(ABTS_1.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))
colnames(ABTS_2.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))
colnames(LDOPA_1.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))
colnames(LDOPA_2.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))
colnames(DMP_1.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))
colnames(DMP_2.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))


# Load samples and treatments, the same template is used for all substrates
samples1 <- read_excel("Data/Raw/Template_samples_20231128_plate1.xlsx", range = "A1:M9")
samples2 <- read_excel("Data/Raw/Template_samples_20231128_plate2.xlsx", range = "A1:M9")
samples2_MBTH <- read_excel("Data/Raw/Template_samples_20231128_plate2_MBTH.xlsx", range = "A1:M9")
treats1 <- read_excel("Data/Raw/Template_samples_20231128_plate1.xlsx", range = "A12:M20")
treats2 <- read_excel("Data/Raw/Template_samples_20231128_plate2.xlsx", range = "A12:M20")
treats2_MBTH <- read_excel("Data/Raw/Template_samples_20231128_plate2_MBTH.xlsx", range = "A12:M20")

# Wrangle data ####

# Wrangle sample and treatment data into long format with "well_ID" as ID-marker

samples1 <-  samples1 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "sample") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

samples2 <-  samples2 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "sample") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

samples2_MBTH <-  samples2_MBTH |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "sample") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

treats1 <- treats1 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "treatment") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

treats2 <- treats2 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "treatment") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

treats2_MBTH <- treats2_MBTH |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "treatment") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

# Join sample and treatment data
samples1 <- full_join(samples1, treats1, by = "well_ID")
samples2 <- full_join(samples2, treats2, by = "well_ID")
samples2_MBTH <- full_join(samples2_MBTH, treats2_MBTH, by = "well_ID")

# Dataframes ready to be merged with absorbance data

####
# Wrangle dataset with absorbance reads into long format
# Round Time [s] to 0 digits
MBTH_1.raw <- MBTH_1.raw |> 
  mutate(time_s = floor(time_s))
MBTH_2.raw <- MBTH_2.raw |> 
  mutate(time_s = floor(time_s))
ABTS_1.raw <- ABTS_1.raw |> 
  mutate(time_s = floor(time_s))
ABTS_2.raw <- ABTS_2.raw |> 
  mutate(time_s = floor(time_s))
LDOPA_1.raw <- LDOPA_1.raw |> 
  mutate(time_s = floor(time_s))
LDOPA_2.raw <- LDOPA_2.raw |> 
  mutate(time_s = floor(time_s))
DMP_1.raw <- DMP_1.raw |> 
  mutate(time_s = floor(time_s))
DMP_2.raw <- DMP_2.raw |> 
  mutate(time_s = floor(time_s))

# Transpose to long format

MBTH_1.raw.long <- MBTH_1.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")
MBTH_2.raw.long <- MBTH_2.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")

ABTS_1.raw.long <- ABTS_1.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")
ABTS_2.raw.long <- ABTS_2.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")

LDOPA_1.raw.long <- LDOPA_1.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")
LDOPA_2.raw.long <- LDOPA_2.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")

DMP_1.raw.long <- DMP_1.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")
DMP_2.raw.long <- DMP_2.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")

MBTH_1.clean <- full_join(MBTH_1.raw.long, samples1, by = "well_ID") 
MBTH_2.clean <- full_join(MBTH_2.raw.long, samples2_MBTH, by = "well_ID") 
# Merge datasets from both plates
MBTH.clean <- rbind(MBTH_1.clean, MBTH_2.clean)

ABTS_1.clean <- full_join(ABTS_1.raw.long, samples1, by = "well_ID") 
ABTS_2.clean <- full_join(ABTS_2.raw.long, samples2, by = "well_ID") 
# Merge datasets from both plates
ABTS.clean <- rbind(ABTS_1.clean, ABTS_2.clean)

LDOPA_1.clean <- full_join(LDOPA_1.raw.long, samples1, by = "well_ID") 
LDOPA_2.clean <- full_join(LDOPA_2.raw.long, samples2, by = "well_ID") 
# Merge datasets from both plates
LDOPA.clean <- rbind(LDOPA_1.clean, LDOPA_2.clean)

DMP_1.clean <- full_join(DMP_1.raw.long, samples1, by = "well_ID")
DMP_2.clean <- full_join(DMP_2.raw.long, samples2, by = "well_ID") 
# Merge datasets from both plates
DMP.clean <- rbind(DMP_1.clean, DMP_2.clean)

# Final polishing of datasets

MBTH.clean <- MBTH.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = sample) |> 
  select(-c(well_ID,cycle_nr, temp, time_s))  #drops uninformative columns

ABTS.clean <- ABTS.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = sample) |> 
  select(-c(well_ID,cycle_nr, temp, time_s))  #drops uninformative columns

LDOPA.clean <- LDOPA.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = sample) |> 
  select(-c(well_ID,cycle_nr, temp, time_s))  #drops uninformative columns

DMP.clean <- DMP.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = sample) |> 
  select(-c(well_ID,cycle_nr, temp, time_s))  #drops uninformative columns


# Export cleaned dataframes
write.csv(MBTH.clean, "Data/Clean/clean_MBTH_needle_masceration.csv", row.names = FALSE)
write.csv(ABTS.clean, "Data/Clean/clean_ABTS_needle_masceration.csv", row.names = FALSE)
write.csv(LDOPA.clean, "Data/Clean/clean_LDOPA_needle_masceration.csv", row.names = FALSE)
write.csv(DMP.clean, "Data/Clean/clean_DMP_needle_masceration.csv", row.names = FALSE)
