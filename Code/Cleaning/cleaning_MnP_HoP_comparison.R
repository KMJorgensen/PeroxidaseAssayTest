# Cleaning MnP and HoP comparison [2023-12-05]

library(tidyverse)
library(readxl)


# Load all datasets 

mbth.raw <- read_excel("Data/Raw/MnP_enzyme_assay_20231205_105857.xlsx", range = "A44:CU104")
#Fix column names
colnames(mbth.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))

abts.raw <- read_excel("Data/Raw/ABTS_peroxidase_assay_20231205_113002.xlsx", range = "A38:CU98")
colnames(abts.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))

ldopa.raw <- read_excel("Data/Raw/L_DOPA_peroxidase_assay_20231205_123154.xlsx", range = "A38:CU98")
colnames(ldopa.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))

dmp.raw <- read_excel("Data/Raw/DMP_peroxidase_assay_20240516_151046.xlsx", range = "A38:CU98")
colnames(dmp.raw)[1:3] <- (c("cycle_nr", "time_s", "temp"))

# Load samples and treatments
# Plate layout 1 is for MBTH and ABTS assay
treats1 <- read_excel("Data/Raw/Template_samples_20231205_1.xlsx", range = "A1:M9")
samples1 <- read_excel("Data/Raw/Template_samples_20231205_1.xlsx", range = "A12:M20")

# Plate layout 2 is for L-DOPA assay
treats2 <- read_excel("Data/Raw/Template_samples_20231205_2.xlsx", range = "A1:M9")
samples2 <- read_excel("Data/Raw/Template_samples_20231205_2.xlsx", range = "A12:M20")

# Plate layout 3 is for DMP assay
treats3 <- read_excel("Data/Raw/Template_samples_20240516_151046.xlsx", range = "A1:M9")
samples3 <- read_excel("Data/Raw/Template_samples_20240516_151046.xlsx", range = "A12:M20")

# Wrangle data ####

# Wrangle sample and treatment data into long format with "well_ID" as ID-marker

samples1 <-  samples1 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "treatment") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

treats1 <- treats1 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "enzyme") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

samples1 <- full_join(samples1, treats1, by = "well_ID")


samples2 <-  samples2 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "treatment") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

treats2 <- treats2 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "enzyme") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

samples2 <- full_join(samples2, treats2, by = "well_ID")

samples3 <-  samples3 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "treatment") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

treats3 <- treats3 |> 
  pivot_longer(cols = 2:13, names_to = "column", values_to = "enzyme") |> 
  unite(well_ID, c(row, column), sep = "", remove = TRUE)

samples3 <- full_join(samples3, treats3, by = "well_ID")
# Dataframe ready to be merged with absorbance data

# Wrangle dataset with absorbance reads into long format
# Round Time [s] to 0 digits
mbth.raw <- mbth.raw |> 
  mutate(time_s = floor(time_s))
abts.raw <- abts.raw |> 
  mutate(time_s = floor(time_s))
ldopa.raw <- ldopa.raw |> 
  mutate(time_s = floor(time_s))
dmp.raw <- dmp.raw |> 
  mutate(time_s = floor(time_s))

mbth.raw.long <- mbth.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")
abts.raw.long <- abts.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")
ldopa.raw.long <- ldopa.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")
dmp.raw.long <- dmp.raw |> 
  pivot_longer(cols = 4:99, names_to = "well_ID", values_to = "absorbance")

mbth.clean <- full_join(mbth.raw.long, samples1, by = "well_ID") 
abts.clean <- full_join(abts.raw.long, samples1, by = "well_ID") 
ldopa.clean <- full_join(ldopa.raw.long, samples2, by = "well_ID") 
dmp.clean <- full_join(dmp.raw.long, samples3, by = "well_ID") 


mbth.clean <- mbth.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = enzyme) |> 
  select(-c(cycle_nr, temp, time_s)) |> #drops uninformative columns
  mutate(substrate = "MBTH/DMAB")
  
abts.clean <- abts.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = enzyme) |> 
  select(-c(cycle_nr, temp, time_s)) |> #drops uninformative columns
  mutate(substrate = "ABTS")  

ldopa.clean <- ldopa.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = enzyme) |> 
  select(-c(cycle_nr, temp, time_s)) |> #drops uninformative columns
  mutate(substrate = "L-DOPA")  

dmp.clean <- dmp.clean |> 
  drop_na() |> # removes readings of empty wells
  mutate(time_min = time_s/60) |> #change time to min instead of s
  relocate(time_min, .before = absorbance ) |> 
  relocate(absorbance, .after = enzyme) |> 
  select(-c(cycle_nr, temp, time_s)) |> #drops uninformative columns
  mutate(substrate = "DMP")  


# Merge all datasets
dat.clean <- rbind(mbth.clean, abts.clean, ldopa.clean, dmp.clean )

# Export cleaned dataframe
write.csv(dat.clean, "Data/Clean/clean_MnP_HoP_comparison.csv", row.names = FALSE)

