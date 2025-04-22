
#B. Bock
#22 April 2025
#Dark Web Paper
#All data prep is in this file, to make the analysis file cleaner.

library(readxl)
library(dplyr)
library(tidyr)
library(nanodRop)
library(ggplot2)

setwd("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web")


# Biomass -----------------------------------------------------------------

#Biomass + metadata, make sure it's in your working directory
ds0 <- read_excel("Datasets/DWN_Data_2023.xlsx", sheet=2)


#Calculate fresh:dry root weight ratio to estimate dry weights of the root subsamples used for scoring.
ds <- ds0%>%
  mutate(A_roots_fd = if_else(
    is.na(A_Root_Subsample_Wet_Weight_g),
    A_Root_Wet_Weight_g / A_Root_Dry_Weight_g,
    (A_Root_Wet_Weight_g - A_Root_Subsample_Wet_Weight_g) / A_Root_Dry_Weight_g),
    B_roots_fd = if_else(
      is.na(B_Root_Subsample_Wet_Weight_g),
      B_Root_Wet_Weight_g / B_Root_Dry_Weight_g,
      (B_Root_Wet_Weight_g - B_Root_Subsample_Wet_Weight_g) / B_Root_Dry_Weight_g
    ))

fd_ratio <- mean(rbind(ds$A_roots_fd, ds$B_roots_fd), na.rm = T)
#Mean 

#Wide dataset up until here.


ds <- ds %>%
  mutate(
    A_subsample_estd_dry_g = A_Root_Subsample_Wet_Weight_g / fd_ratio,
    B_subsample_estd_dry_g = B_Root_Subsample_Wet_Weight_g / fd_ratio,
    A_Root_Dry_Weight_g = A_Root_Dry_Weight_g + coalesce(A_subsample_estd_dry_g, 0),
    B_Root_Dry_Weight_g = B_Root_Dry_Weight_g + coalesce(B_subsample_estd_dry_g, 0)
  ) %>%
  pivot_longer(
    cols = matches("^[AB]_"),  # all columns starting with A_ or B_
    names_to = c("Chamber", ".value"),
    names_pattern = "^([AB])_(.*)"
  )%>%
  mutate(Experiment_Round = as.character(Experiment_Round))%>%
  filter(!is.na(Type_Barrier))

ds_long <- ds %>%
  pivot_longer(cols = c("Shoot_Dry_Weight_g", "Root_Dry_Weight_g"), names_to = "Compartment", values_to = "Dry_Weight_g", names_pattern = "(.*)_Dry_Weight_g")%>%
  mutate(Compartment_abbr = case_when(
    Compartment == "Root" ~ "R",
    Compartment == "Shoot" ~ "S",
    TRUE ~ Compartment  # Default case if other values exist
  )) %>%
  mutate(Sample_Name = paste0(Experiment_Round, ".", Box_Nr, Chamber, ".", Compartment_abbr))%>%
  relocate(Sample_Name) %>%
  mutate(Chamber = case_match(Chamber,
                              "A" ~ "Donor",
                              "B" ~ "Receiver",
                              .default = Chamber))%>%
  select(!contains(c("Root", "Shoot", "subsample")))

write.csv(ds, "Datasets/biomass_ds_wide.csv")
write.csv(ds_long, "Datasets/biomass_ds_long.csv")

dist <- combo_ds2 %>%
  group_by(Experiment_Round, Type_Barrier)%>%
  summarize(n = n()/2)
dist
#write.csv(dist, "Datasets/dist_boxes.csv")




# Spectral ----------------------------------------------------------------

setwd("Datasets/NanoDrop Readings/Spectrophotometer_Reads")

key <- read_excel("Key.xlsx")%>%
  mutate(samps = as.character(Tube)) 

key2 <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/biomass_ds_long.csv") #Has the rest of the data. Not sure if should choose wide or long here.

#Supernatant, read on same day as centrifuged
ds1 <- read_nanodrop("UV-Vis 5_21_2024 5_40_43 PM 1-10.tsv")%>%
  mutate(n = "sn1")%>%
  mutate(samps = case_when(samps == "6" ~ NA, #First reading of 6 had an error
                           samps == "6.2" ~ "6", #6.2 is the second reading of 6 and does not appear to have the same issue
                           .default = samps))%>%
  filter(!is.na(samps))

ds2 <- read_nanodrop("UV-Vis 5_23_2024 2_47_05 PM.tsv")%>% 
  mutate(n = "sn2") #Hadn't included this one for some reason

ds3 <- read_nanodrop("UV-Vis 5_29_2024 6_21_54 PM.tsv")%>% 
  mutate(n = "sn3")

ds4 <- read_nanodrop("UV-Vis 5_31_2024 6_11_04 PM.tsv")%>%
  mutate(n = "sn4")

# ds5 <- read_nanodrop("UV-Vis 1_13_2025 1_16_18 PM.tsv")%>%
#   mutate(n = "sn5")%>%
#   filter(!samps %in% c("h1", "h2", "s1", "s2", "s3", "s4", "147.2", "147.19999999999999")) #Remove the standards I used. Can refer back to later if wanted.
#Commenting out because I redid these samples in ds6 (below)

ds6 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/Dark_Web_Samples/UV-Vis 1_16_2025 12_49_48 PM.tsv")%>% 
  mutate(n = "sn7",
         samps = case_when(samps == "149" ~ NA, #First and second reading of 149 had an error
                           samps == "149.2" ~ NA,
                           samps == "156.3" ~ NA,
                           samps == "150" ~ NA,
                           samps == "149.3" ~ "149", #6.2 is the second reading of 6 and does not appear to have the same issue
                           samps == "150.3" ~ NA,
                           .default = samps))%>%
  filter(!is.na(samps))

# Start analyzing spectral data -------------------------------------------



peaks <- ds1 %>%
  full_join(ds3, by = colnames(ds1))%>%
  full_join(ds4, by = colnames(ds1))%>%
  full_join(ds6, by = colnames(ds1)) %>%
  full_join(key, by = "samps") %>%
  mutate(Chamber = case_when(grepl("A", Sample_Name) ~ "Donor",
                             grepl("B", Sample_Name) ~ "Receiver", 
                             .default = NA
  ),
  Compartment = case_when(grepl("R", Sample_Name) ~ "Root",
                          grepl("S", Sample_Name) ~ "Shoot", 
                          .default = NA
  )) %>%  
  full_join(key2, by = c("Sample_Name", "Chamber", "Compartment"))%>% 
  dplyr::select(c("samps", "abs", "Sample_Name", "Chamber", "Compartment", "Box_Nr", "Type_Barrier", "Experiment_Round", "Weight_G", "uL_H2O", "waves", "Dry_Weight_g"))%>%
  filter(!is.na(Compartment))%>%
  filter(!is.na(Type_Barrier))%>%
  mutate(Experiment_Round = recode(Experiment_Round,
                                   `2` = "Preliminary",
                                   `4` = "Main",
                                   `5` = "Follow-Up"),
         Experiment_Round = factor(Experiment_Round,
                                   levels = c("Preliminary", "Main", "Follow-Up"))
  )

peaks %>%
  group_by(Experiment_Round, Type_Barrier)%>%
  summarize(n=n())

#The "average root" and "average shoot" are readings using the same method as in this experiment, but the plants were not treated with dye. All other conditions are the same.
avg_root <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/avg_root.csv")
avg_shoot <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/avg_shoot.csv")


peaks <- peaks %>%
  filter(!is.na(abs))%>%
  mutate(abs_adj = case_when(
    Compartment == "Shoot" ~ (abs - avg_shoot$abs),
    .default = abs))

#Create a dataset of all the NanoDrop data
write.csv(peaks, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/DW_NanoDrop_all.csv")


abs545 <- peaks %>%
  filter(waves== 545) %>%
  #rename("waves545" = abs)%>%
  mutate(uL_H2O = as.numeric(uL_H2O),
         Weight_G = as.numeric(Weight_G))


abs545 %>%
  group_by(Experiment_Round, Type_Barrier, Box_Nr)%>%
  summarize(n=n())

#Make a dataset of all the absorbances at 545 nm
write.csv(abs545, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/abs_545.csv")

shoots <- abs545 %>%
  filter(Compartment == "Shoot")

roots <- abs545 %>%
  filter(Compartment == "Root")

shoots %>%
  ggplot()+
  geom_histogram(aes(x = Weight_G), bins = 10, stat = "count")


shoots %>%
  group_by(Experiment_Round, Type_Barrier)%>%
  summarize(n=n())

shoots %>%
  ggplot()+
  geom_histogram(aes(x = uL_H2O))

#mods and mods2
