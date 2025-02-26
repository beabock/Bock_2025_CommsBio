#B. Bock
#May-June 2024
#This code creates a singular dataset for the controls I used to create the model to predict amount of dye based on the spectral readings of the supernatant of the ground plant samples in water

library(readr)
library(readxl)
library(stringr)
library(ggplot2)
library(dplyr)
library(nanodRop)
library(rempsyc)
library(ggpmisc) #Add lines and equations to ggplot
library(tidyr) #separate function

theme_set(theme_bw())

controls1 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_Data_Controls_5_2024/UV-Vis 5_7_2024 12_48_38 PM.tsv") %>%
  mutate(samps = case_match(samps, 
                            "Sample 1" ~ "1f"))
  
controls2 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_Controls_5-8-24/UV-Vis 5_8_2024 12_53_35 PM.tsv")

test1 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/NanodropOne_bb_test2_4-25-24/bea UV-Vis 4_25_2024 6_16_34 PM.tsv")

test2 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/NanodropOne_bb_test2_4-25-24/bea #2 UV-Vis 4_25_2024 7_56_07 PM.tsv")

ds <- rbind(controls1, controls2, test1, test2)

ds <- ds %>% 
  mutate(waves = as.numeric(waves),
         abs = as.numeric(abs),
         samps = case_match(samps,
                            "1x 1" ~ "10f",
                            "1x 2" ~ "100f",
                            "1x 3" ~ "1000f", .default = samps))

#Looking at just the blanks
blanks <- ds %>%
  filter(str_detect(samps, 'f'))%>%
  mutate(samps = factor(samps, levels = c("1f", "2f 1", "4f 1", "8f 1", "10f","16f 1", "32f 1", "64f 1", "100f", "128f 1", "256f 1", "1000f")))

#1f dilution is weird - too concentrated I think.


# Better controls ---------------------------------------------------------


n1 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_5-14-24/UV-Vis 5_14_2024 5_03_11 PM.tsv")%>% #Bare plant samples, need to remove a couple where the readings were flawed
  mutate(n = "n1")%>%
  mutate(samps2= case_when(grepl("6", samps) ~ "5R",
                           grepl("7", samps) ~ "5R",
                           grepl("8", samps) ~ "5R",
                           grepl("9", samps) ~ "2S",
                           grepl("10", samps) ~ "idk",
                           grepl("11", samps)~ "D4S",
                           grepl("12", samps) ~ "D5S",
                           grepl("13", samps) ~ "D5R",
                           grepl("14", samps) ~ "D4R",
                           grepl("15", samps) ~ "S1",
                           grepl("16", samps) ~ "R4",
                           .default = samps),
         Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", 
                                 .default = NA
         ))%>%
  filter(samps != "15 2" & samps != "9 2") %>% #Removing these because their readings were erroneous
  mutate(Dye = case_when(
    grepl("^(6|7|8|9|1[0-6])", samps) ~ 0,  # Matches 6-16
    grepl("^(1[7-9]|2[0-7])", samps) ~ 1,   # Matches 17-27
    .default = NA                          # Default to NA for any other cases
  ))
  

blanks <-blanks%>% #Keep the blanks for now - note that they are dissimilar from the other datasets in that they have no plant material in them.
  mutate(Dye = case_when(
    samps == "1f" ~ 1,
    samps ==  "2f 1" ~ 1/2,
    samps ==  "4f 1" ~ 1/4,
    samps ==  "8f 1" ~ 1/8,
    samps ==  "16f 1" ~ 1/16,
    samps == "32f 1" ~ 1/32,
    samps ==  "64f 1" ~ 1/64,
    samps ==  "128f 1" ~ 1/128,
    samps == "256f 1" ~ 1/256,
    samps == "1000f" ~ 1/1000,
    .default = NA                          # Default to NA for any other cases
  ))


n1 %>%
  group_by(Dye)%>%
  summarize(n=n())


n2 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_5-14-24/UV-Vis 5_14_2024 5_27_12 PM.tsv")%>% 
  mutate(n = "n2")%>%
  mutate(samps = case_when(samps == "5/14/2024 5:31 PM" ~ "24", .default = samps))%>% #Sample name is missing from this one so it was listed as a date instead
  mutate(samps2= case_when(grepl("17", samps) ~ "5R",
                           grepl("18", samps) ~ "5R",
                           grepl("19", samps) ~ "5R",
                           grepl("20", samps) ~ "2S",
                           grepl("21", samps) ~ "idk",
                           grepl("22", samps)~ "D4S",
                           grepl("23", samps) ~ "D5S",
                           grepl("24", samps) ~ "D5R",
                           grepl("25", samps) ~ "D4R",
                           grepl("26", samps) ~ "S1",
                           grepl("27", samps) ~ "R4",
                           .default = samps),
         Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", .default = samps2))%>%
  filter(samps != "15 2" & samps != "9 2" & samps != "")%>% #Remove erroneous readings
  mutate(Dye = case_when(
    grepl("^(6|7|8|9|1[0-6])", samps) ~ 0,  # Matches 6-16
    grepl("^(1[7-9]|2[0-7])", samps) ~ 1/32,   # Matches 17-27
    TRUE ~ NA_real_                          # Default to NA for any other cases
  ))

n2 %>%
  group_by(samps2, Dye)%>%
  summarize(n=n())

n3 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_5-14-24/UV-Vis 5_15_2024 2_08_18 PM.tsv")%>% 
  mutate(n = "n3")%>%
  mutate(samps2= case_when(grepl("39", samps) ~ "5R",
                           grepl("40", samps) ~ "5R",
                           grepl("41", samps) ~ "5R",
                           grepl("42", samps) ~ "2S",
                           grepl("43", samps) ~ "idk",
                           grepl("44", samps)~ "D4S",
                           grepl("45", samps) ~ "D5S",
                           grepl("46", samps) ~ "D5R",
                           grepl("47", samps) ~ "D4R",
                           grepl("48", samps) ~ "S1",
                           grepl("49", samps) ~ "R4",
                           .default = samps),
         Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", .default = samps2))%>%
  mutate(Dye = case_when(
    grepl("^(39|4[0-9])", samps) ~ 1/16,  # Matches 39-49
    TRUE ~ NA_real_                          # Default to NA for any other cases
  ))

n3 %>%
  group_by(samps2, Dye)%>%
  summarize(n=n())

n4 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_5-14-24/UV-Vis 5_16_2024 12_41_46 PM.tsv")%>% 
  mutate(n = "n4")%>%
  mutate(samps2= case_when(grepl("50", samps) ~ "5R",
                           grepl("51", samps) ~ "5R",
                           grepl("52", samps) ~ "5R",
                           grepl("53", samps) ~ "2S",
                           grepl("54", samps) ~ "idk",
                           grepl("55", samps)~ "D4S",
                           grepl("56", samps) ~ "D5S",
                           grepl("57", samps) ~ "D5R",
                           grepl("58", samps) ~ "D4R",
                           grepl("59", samps) ~ "S1",
                           grepl("60", samps) ~ "R4",
                           grepl("61", samps) ~ "5R",
                           grepl("62", samps) ~ "5R",
                           grepl("63", samps) ~ "5R",
                           grepl("64", samps) ~ "2S",
                           grepl("65", samps) ~ "idk",
                           grepl("66", samps)~ "D4S",
                           grepl("67", samps) ~ "D5S",
                           grepl("68", samps) ~ "D5R",
                           grepl("69", samps) ~ "D4R",
                           grepl("70", samps) ~ "S1",
                           grepl("71", samps) ~ "R4",
                           .default = samps),
         Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", .default = samps2))%>%
  mutate(Dye = case_when(
    grepl("^(5[0-9]|60)", samps) ~ 1/8,  
    grepl("^(6[1-9]|7[0-1])", samps) ~ 1,
    TRUE ~ NA_real_                          # Default to NA for any other cases
  ))

n4 %>%
  group_by(samps2, Dye)%>%
  summarize(n=n())

n5 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_5-14-24/UV-Vis 5_17_2024 2_45_58 PM.tsv")%>% 
  mutate(n = "n5")%>%
  mutate(samps2= case_when(grepl("72", samps) ~ "5R",
                           grepl("73", samps) ~ "5R",
                           grepl("74", samps) ~ "5R",
                           grepl("75", samps) ~ "2S",
                           grepl("76", samps) ~ "idk",
                           grepl("77", samps)~ "D4S",
                           grepl("78", samps) ~ "D5S",
                           grepl("79", samps) ~ "D5R",
                           grepl("80", samps) ~ "D4R",
                           grepl("81", samps) ~ "S1",
                           grepl("82", samps) ~ "R4",
                           grepl("83", samps) ~ "5R",
                           grepl("84", samps) ~ "5R",
                           grepl("85", samps) ~ "5R",
                           grepl("86", samps) ~ "2S",
                           grepl("87", samps) ~ "idk",
                           grepl("88", samps)~ "D4S",
                           grepl("89", samps) ~ "D5S",
                           grepl("90", samps) ~ "D5R",
                           grepl("91", samps) ~ "D4R",
                           grepl("92", samps) ~ "S1",
                           grepl("93", samps) ~ "R4",
                           grepl("94", samps) ~ "5R",
                           grepl("95", samps) ~ "5R",
                           grepl("96", samps) ~ "5R",
                           grepl("97", samps) ~ "2S",
                           grepl("98", samps) ~ "idk",
                           grepl("99", samps)~ "D4S",
                           grepl("100", samps) ~ "D5S",
                           grepl("101", samps) ~ "D5R",
                           grepl("102", samps) ~ "D4R",
                           grepl("103", samps) ~ "S1",
                           grepl("104", samps) ~ "R4",
                           .default = samps),
         Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", .default = samps2))%>%
  mutate(Dye = case_when(
    grepl("^(7[2-9]|8[0-2])", samps) ~ 1/2,  
    grepl("^(8[3-9]|9[0-3])", samps) ~ 3/4,
    grepl("^(9[4-9]|10[0-4])", samps) ~ 9/10,
    TRUE ~ NA_real_                          # Default to NA for any other cases
  ))

n5 %>%
  group_by(samps2, Dye)%>%
  summarize(n=n())

#The pipettor seemed off in volume expelled when setting up n5, and the results of n5 do seem different enough from n6 (same samples repeated with different pipettor) that n5 will be removed from the analyses.

n6 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_5-14-24/UV-Vis 5_17_2024 5_49_34 PM.tsv")%>% 
  mutate(n = "n6")%>%
  mutate(samps2= case_when(grepl("105", samps) ~ "5R",
                           grepl("106", samps) ~ "5R",
                           grepl("107", samps) ~ "5R",
                           grepl("108", samps) ~ "2S",
                           grepl("109", samps) ~ "idk",
                           grepl("110", samps)~ "D4S",
                           grepl("111", samps) ~ "D5S",
                           grepl("112", samps) ~ "D5R",
                           grepl("113", samps) ~ "D4R",
                           grepl("114", samps) ~ "S1",
                           grepl("115", samps) ~ "R4",
                           grepl("116", samps) ~ "5R",
                           grepl("117", samps) ~ "5R",
                           grepl("118", samps) ~ "5R",
                           grepl("119", samps) ~ "2S",
                           grepl("120", samps) ~ "idk",
                           grepl("121", samps)~ "D4S",
                           grepl("122", samps) ~ "D5S",
                           grepl("123", samps) ~ "D5R",
                           grepl("124", samps) ~ "D4R",
                           grepl("125", samps) ~ "S1",
                           grepl("126", samps) ~ "R4",
                           grepl("127", samps) ~ "5R",
                           grepl("128", samps) ~ "5R",
                           grepl("129", samps) ~ "5R",
                           grepl("130", samps) ~ "2S",
                           grepl("131", samps) ~ "idk",
                           grepl("132", samps)~ "D4S",
                           grepl("133", samps) ~ "D5S",
                           grepl("134", samps) ~ "D5R",
                           grepl("135", samps) ~ "D4R",
                           grepl("136", samps) ~ "S1",
                           grepl("137", samps) ~ "R4",
                           .default = samps),
         Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", .default = samps2))%>%
   filter(samps != "110")%>% #samp 110 had an error in the reading
  mutate(Dye = case_when(
    grepl("^(10[5-9]|11[0-5])", samps) ~ 1/2,  
    grepl("^(11[6-9]|12[0-6])", samps) ~ 3/4,
    grepl("^(12[7-9]|13[0-7])", samps) ~ 9/10,
    TRUE ~ NA_real_                          # Default to NA for any other cases
  ))

n6 %>%
  group_by(samps2, Dye)%>%
  summarize(n=n())

#Exlcuding n5 
dataset <- rbind(n1, n2, n3, n4, n6)%>%
  mutate(plant_g = case_when(samps2 == "5R" ~ 0.0015, #g
                             samps2 == "2S" ~ 0.0006,
                             samps2 == "D4S" ~ 0.0003,
                             samps2 == "D5S" ~ 0.0007,
                             samps2 == "D5R" ~ 0.0007,
                             samps2 == "D4R" ~ 0.0008,
                             samps2 == "S1" ~ 0.0005,
                             samps2 == "R4" ~ 0.0013,
                             .default = Dye),
         uL = 30) #uL H2O added to each sample


key <- read_excel("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/Dark_Web_Samples/Key.xlsx")

n7 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/Dark_Web_Samples/UV-Vis 1_16_2025 12_49_48 PM.tsv")%>% 
  mutate(n = "n7")%>%
  filter(!samps %in% c(149, 150, 156, 149.2, 149.3, 150.3, 156.3, 157.2))%>% #These are actual samples, not controls
  left_join(key, by = c("samps" = "Tube"))%>%
  mutate(uL = as.numeric(uL_H2O),
         Dye = Dye_Dilution,
    samps2= case_when(grepl("151", Sample_Name) ~ "S1",
                           grepl("152", Sample_Name) ~ "S2",
                           grepl("155", Sample_Name) ~ "DS5"),
    plant_g = case_when(grepl("151", Sample_Name) ~ 0.0011,
                            grepl("152", Sample_Name) ~ 0.0009,
                            grepl("155", Sample_Name) ~ 0.0009))%>%
  mutate(Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", .default = samps2))


dataset <- full_join(dataset, n7)

#Do avg shoot here

#Average plant without dye

dataset <- dataset %>%
  mutate(dye_ug = 
           Dye #Dilution of dye
         * 0.2 #.2 mg/mL of acid fuchsin in soln
         * 2/3 #further dilution of dye in samples run
         * 0.001 #mL of solution added to nanodrop
         * 1000 #Convert from mg to ug
  ) #to get ug of dye

avg_shoot <- dataset %>%
  filter(Dye == 0 & Compartment == "Shoots") %>% #Compartment won't ne root for many of the 0 samples
  group_by(waves) %>%
  summarize(abs = mean(abs, na.rm = T))

write.csv(avg_shoot, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/avg_shoot.csv")

avg_root <- dataset %>%
  filter(Dye == 0 & Compartment == "Roots") %>% #Compartment won't ne root for many of the 0 samples
  group_by(waves) %>%
  summarize(abs = mean(abs, na.rm = T))

write.csv(avg_root, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/avg_root.csv")

# zeros <- dataset %>%
#   filter(dye_ug == 0)


#write.csv(zeros, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/zero_dye_samples.csv")

dataset <- dataset %>%
  filter(!samps %in% c("9 1","12 2", "9 3", "11 1", "15 3", "20", "158", "185", "192", "67", "132", "136", "48", "59", "53", "66", "15 1", "23", "42"))#%>%
 # filter(Compartment == "Shoots")

#Let's look just at the wavelengths where the blanks peaked
peaks <- dataset %>%
  filter(waves== 600 |
           waves == 548 |
           waves == 545 |
           waves == 537 | 
           waves == 492 |
           waves == 443 |
           waves == 386 |
           waves == 300 |
           waves == 290 |
           waves == 260 |
           waves == 210 | 
           waves == 190 )


#Shoots start to shwo the dye signal with the 1x dye added

#0.2 mg per mL of acid fuchsin in solution

#1/16 dilution means it has 0.0125 mg/ml dye

#2 uL dye plus 1 microliter plant = 2/3 dilution of above * 0.001 to get mL
#1 uL of that 

dataset %>%
  group_by(Dye)%>%
  summarize(n = n())


dataset %>%
  group_by(dye_ug)%>%
  summarize(n = n())

dataset %>%
  ggplot(aes(x = dye_ug))+
  geom_histogram(bins = 10)+
  facet_grid(~Compartment)

test <- dataset %>%
  filter(Compartment == "Shoots")%>%
  filter(waves ==545)

shapiro.test(test$dye_ug)
#Shoots data are normally distd

test <- dataset %>%
  filter(Compartment == "Roots")%>%
  filter(waves ==545)

shapiro.test(test$dye_ug)
#Roots not so much.

# dataset <- dataset %>%
#   group_by(samps)%>%
#   mutate(abs_adj = case_when(Compartment == "Shoots" ~ (abs - avg_shoot$abs),
#                              .default = abs))
#Ended up not doing this.


blanks <- blanks %>%
  mutate(dye_ug = round(
           Dye #Dilution of dye
         * 0.2 #.2 mg/mL of acid fuchsin in soln
         * 2/3 #further dilution of dye in samples run
         * 0.001 #mL of solution added to nanodrop
         * 1000,
         digits = 4) #Convert from mg to ug
  ) #to get ug of dye


#Plot for supplementary showing why we picked 545 nm
blanks %>%
  filter(samps != "1f")%>% #1f is too concentrated, getting erroneous readings
  ggplot(aes(x = waves, y= abs))+
  geom_line(aes(color = factor(dye_ug)))+
  # geom_smooth(color = "#00AFBB")+
  labs(x = "Wavelength (nm)", y = "Absorbance", color = "uG Dye")+
  geom_vline(xintercept = 545)+
  annotate(geom="text", x=700, y=25, label="Vertical Line at 545nm",
           color="black")+
  theme(
        text = element_text(size = 18))





shoots <- dataset %>%
  filter(waves == 545)%>% #training data shouldn't have abs > 5 bc the data do not go above 5
  filter(Compartment == "Shoots",
         abs < 5
  )%>%
  mutate(g_ml = plant_g/uL*1000,
         abs_g_ml = abs/g_ml)


roots <- dataset %>%
  filter(waves == 545)%>% #training data shouldn't have abs > 5 bc the data do not go above 5
  filter(Compartment == "Roots",
         abs < 5
  )%>%
  mutate(g_ml = plant_g/uL*1000,
         abs_g_ml = abs/g_ml)

  write.csv(roots, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/roots_mod.csv")
  write.csv(shoots, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/shoots_mod.csv")

# roots <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Datasets/roots_mod.csv")
# shoots <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Datasets/shoots_mod.csv")

avg_root_abs<- dataset %>%
  filter(Dye == 0)%>%
  filter(Compartment == "Roots",
         waves > 450,
         waves < 650)%>%
    mutate(g_ml = plant_g/uL*1000,
           abs_g_ml = abs/g_ml)%>%
   ggplot(aes(y = abs, x = waves, color = samps))+
  geom_line()+
  theme(legend.position = "none", legend.title.position = "none")+
  labs(x = "Wavelength (nm)", y = "Absorbance")+
  geom_vline(xintercept= 545, linetype = "dashed")

ggsave("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/avg_root_dye.png", avg_root_abs, dpi = 1000, width = 3, height = 2.5)


# Create model ------------------------------------------------------------

shoots %>%
  ggplot(aes(x = dye_ug))+
  geom_histogram()

#Manually need to increase sampling of dyed samples. More in the 0.07 realm (10 or so) and between 0.02 and 0.1. Around 0.05 is really missing.
shoots %>%
  ggplot(aes(y = dye_ug))+
  geom_point(aes(x = abs), color = "blue")+
#  geom_point(aes(x = abs_adj), color = "red")+
  geom_smooth(aes(x = abs), color = "blue", method = "lm")#+
 # geom_smooth(aes(x = abs_adj), color = "red", method = "lm")
#Maybe dont need to adjust for average root/shoots.


mod1 <- lm(dye_ug ~ abs, shoots) #Best
summary(mod1)

mod2 <- lm(dye_ug ~ abs+plant_g, shoots) #Maybe better to use, slightly higher r2
summary(mod2)

mod5 <- lm(dye_ug ~ abs*plant_g, shoots) #Lower F stat but higher p-val. need to think on this.
summary(mod5)

mod3 <- lm(dye_ug ~ abs + g_ml, shoots) #not great
summary(mod3)

mod4 <-lm(dye_ug ~ abs_g_ml, shoots) #Worst
summary(mod4)


#Roots

mod4 <- lm(dye_ug ~ abs, roots) 
summary(mod4) #Not as good as shoots...









# Suspended samples, tested to see what it looks like. Not helpful. -------------------------------------------------------

#s datasets are suspended; n are centrifuged and the supernatant is read
#If "e" is in a sample, it means I suspect the previous read of that sample was in error
s1 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/BB_NanoDrop_5-14-24/UV-Vis 5_23_2024 4_53_54 PM.tsv")%>%
  mutate(n = "s1")%>% #suspended 1
  mutate(samps2= case_when(grepl("c6", samps) ~ "5R",
                           grepl("c7", samps) ~ "5R",
                           grepl("c8", samps) ~ "5R",
                           grepl("c9", samps) ~ "2S",
                           grepl("c10", samps) ~ "idk",
                           grepl("c11", samps)~ "D4S",
                           grepl("c12", samps) ~ "D5S",
                           grepl("c13", samps) ~ "D5R",
                           grepl("c14", samps) ~ "D4R",
                           grepl("c15", samps) ~ "S1",
                           grepl("c16", samps) ~ "R4",
                           .default = samps),
         Compartment = case_when(grepl("R", samps2) ~ "Roots",
                                 grepl("S", samps2) ~ "Shoots", 
                                 .default = NA
         ))%>%
  mutate(Dye = 0)

#Suspended real samples
key <- read_excel("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/Dark_Web_Samples/Key.xlsx")%>%
  mutate(samps = as.character(Tube))

#Suspended real samples
#Iteration 2 is supernatant readings,
#Iteration 3 is in suspension
s2 <- read_nanodrop("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/NanoDrop Readings/Dark_Web_Samples/UV-Vis 5_23_2024 2_47_05 PM.tsv")%>%
  separate(samps, c("samps", "Iteration"))%>%
  full_join(key, by = "samps")%>%
  mutate(n = "s1",
         Compartment = case_when(grepl("R", Sample_Name) ~ "Roots",
                                 grepl("S", Sample_Name) ~ "Shoots", 
                                 .default = NA
         )) %>% #suspended 1
  filter(Iteration == "2")


#Same thing to suspended dataset
s1 <- s1 %>%
  mutate(plant_g = case_when(samps2 == "5R" ~ 0.0015/30, #g/uL
                             samps2 == "2S" ~ 0.0006/30,
                             samps2 == "D4S" ~ 0.0003/30,
                             samps2 == "D5S" ~ 0.0007/30,
                             samps2 == "D5R" ~ 0.0007/30,
                             samps2 == "D4R" ~ 0.0008/30,
                             samps2 == "S1" ~ 0.0005/30,
                             samps2 == "R4" ~ 0.0013/30,
                             .default = 0),
         abs_to_plant = abs/plant_g)

avg_root_s <- s1 %>%
  filter(Dye == 0 & Compartment == "Roots") %>%
  group_by(waves) %>%
  summarize(abs = mean(abs, na.rm = T))

avg_shoot_s <- s1 %>%
  filter(Dye == 0 & Compartment == "Shoots") %>%
  group_by(waves) %>%
  summarize(abs = mean(abs, na.rm = T))


s2 <- s2 %>%
  mutate(abs_adj = case_when(Compartment == "Roots" ~ (abs - avg_root_s$abs),
                             Compartment == "Shoots" ~ (abs - avg_shoot_s$abs),
                             .default = abs),
         abs_adj_plus_plant = abs_adj/Weight_G)



s2 %>%
  ggplot(aes(x = waves, y = abs, color = samps))+
  geom_point()

dataset %>%
  filter(Compartment != "idk" | Compartment =="Roots" | Dye != 0) %>%
  ggplot(aes(x=waves, y = abs))+
  geom_line(aes(color = factor(Dye)))+
  stat_summary(geom = "line", fun = "mean")#+
facet_grid(~Compartment)

ggplot()+
  geom_point(data = s2, aes(x = waves, y = abs, color = samps))+
  stat_summary(data = s1, aes(x = waves, y = abs), geom = "line", fun = "mean", linetype = "dashed")

# Misc code that isn't useful anymore -----------------------------------------------------------
#Note: This next section of data is a bit weird and maybe should not be used

# pattern <- "\\s+(.*)$" #This one seems to work
# 
# plants <- ds%>%
#          mutate(id = trimws(str_extract(samps, pattern)),
#                 control = case_when(grepl("f", samps) ~ "Standard",
#                                     grepl("p", samps) ~ "Plant", .default = id))%>%
#   mutate(id2 = case_match(id, 
#                              "1" ~ "Shoots1",
#                            "2" ~ "Shoots2",
#                            "3" ~ "Shoots3+dilutepink",
#                            "4" ~ "Roots4",
#                            "5" ~ "Roots5",
#                            "6" ~ "Shoots2 w/ 1/128",
#                            "7" ~ "Roots5 w/ 1/128",
#                            "8" ~ "Shoots2 w/ 1/32",
#                            "9" ~ "Roots5 w/  1/32",
#                            "10" ~ "Shoots2 w/ 1/16",
#                            "11" ~ "Roots5 w/ 1/16",
#                            "12" ~ "Shoots2 w/ 1/8",
#                            "13" ~ "Roots5 w/ 1/8",
#                            "14" ~ "Shoots2 w/ 1/2",
#                            "15" ~ "Roots5 w/ 1/2",
#                            .default = id),
#                 Compartment = case_when(grepl("Shoots", id2) ~ "Shoots",
#                                         grepl("Roots", id2) ~ "Roots", .default = id),
#                 Dye = case_when(grepl("1/128", id2) ~ "1/128",
#                                 grepl("1/32", id2) ~ "1/32",
#                                 grepl("1/16", id2) ~ "1/16",
#                                 grepl("1/8", id2) ~ "1/8",
#                                 grepl("1/2", id2) ~ "1/2",
#                                 control == "Standard" ~ "Standard",
#                                 .default = "0"))%>%
#   filter(control != "Standard") #Need to do this better



# Trying to account for dilutions of samples ------------------------------
#Note: May want to skip to next section

#Also need to account for the dilutions of plants

# 
# x <- .0013/15 #starting dilution of plant, g/uL
# x2 <- x - 1.5 + 1*(1/128)
# x3 <- x - 1.5 + 1.5*(1/32)
# 
# V_start <- 15
# V_used <- 1.5
# C_start <- 8.6
# V_dilution <- 1
# V_total<- V_start-V_used + V_dilution
# 
# 
# C_1_128 <- ((V_start-V_used)*C_start+V_dilution)/V_total
# 
# V_start_1_32 <- V_total
# V_used <- 1.5
# C_start_1_32 <- C_1_128
# V_dilution_1_32 <- 1.5
# V_total_1_32<- V_start_1_32-V_used + V_dilution_1_32
# 
# 
# C_end_1_32 <- ((V_start_1_32-V_used)*C_start_1_32+V_dilution_1_32)/V_total_1_32
# 
# 
# V_start_1_16 <- V_total_1_32
# C_start_1_16 <- C_end_1_32
# V_dilution_1_16 <- 1.5
# V_total_1_16<- V_start_1_16-V_used + V_dilution_1_16
# 
# 
# C_end_1_16 <- ((V_start_1_16-V_used)*C_start_1_16+V_dilution_1_16)/V_total_1_16
# 
# 
# 
# V_start_1_8 <- V_total_1_16
# C_start_1_8 <- C_end_1_16
# V_dilution_1_8 <- 1.5
# V_total_1_8<- V_start_1_8-V_used + V_dilution_1_8
# 
# C_end_1_8 <- ((V_start_1_8-V_used)*C_start_1_8+V_dilution_1_8)/V_total_1_8
# 
# 
# 
# V_start_1_2 <- V_total_1_8
# C_start_1_2 <- C_end_1_8
# V_dilution_1_2 <- 1.5
# V_total_1_2<- V_start_1_2-V_used + V_dilution_1_2
# 
# C_end_1_2 <- ((V_start_1_2-V_used)*C_start_1_2+V_dilution_1_2)/V_total_1_2
# 
# dilution_factor <- 0.91 #Or 1.098
# 
# ds_adj <- plants %>%
#   mutate(abs_adj = case_when(Dye == "1/128" ~ abs*(C_start/C_1_128), 
#                              Dye == "1/32" ~ abs*C_start/C_end_1_32,
#                              Dye == "1/16" ~ abs*C_start/C_end_1_16,
#                              Dye == "1/8" ~ abs*C_start/C_end_1_8,
#                              Dye == "1/2" ~ abs*C_start/C_end_1_2,.default = abs))
# 
# 
# 
# #Average root without dye
# avg_root <- ds_adj %>%
#   filter(Dye == "0" & Compartment == "Roots") %>%
#   group_by(waves) %>%
#   summarize(abs = mean(abs, na.rm = T))
# 
# avg_root %>%
#   ggplot(aes(x=waves, y = abs))+
#   geom_point()
# 
# #Average shoot without dye
# avg_shoot <- ds_adj %>%
#   filter(Dye == "0" & Compartment == "Shoots") %>%
#   group_by(waves) %>%
#   summarize(abs = mean(abs, na.rm = T))
# 
# avg_shoot %>%
#   ggplot(aes(x=waves, y = abs))+
#   geom_point()
# 
# #Now subtract the averages from the readings
# 
# ds_adj_3 <- ds_adj %>%
#    mutate(abs_adj_2 = case_when(Compartment == "Roots" ~ (abs_adj - avg_root$abs),
#                                Compartment == "Shoots" ~ (abs_adj - avg_shoot$abs),
#                                .default = abs_adj))%>%
#   mutate(Dye = factor(Dye, levels = c("0", "1/128", "1/32", "1/16", "1/8", "1/2")))

#Now plot, it's messy

#not sure my mathematical adjustments work - assuming that absorbance scales with dilution, which it may not. Going to redo test re-dyed samples. 



#Below I'm checking the color-blind-friendliness of the colors I use for most plots. They pass for the four main types of color-blindess.
# 
# library(dichromat)
# red_green_colors <- c("#FC4E07", "#00AFBB")
# 
# # convert to the three dichromacy approximations
# protan <- dichromat(red_green_colors, type = "protan")
# deutan <- dichromat(red_green_colors, type = "deutan")
# tritan <- dichromat(red_green_colors, type = "tritan")
# 
# # plot for comparison
# layout(matrix(1:4, nrow = 4)); par(mar = rep(1, 4))
# recolorize::plotColorPalette(red_green_colors, main = "Trichromacy")
# recolorize::plotColorPalette(protan, main = "Protanopia")
# recolorize::plotColorPalette(deutan, main = "Deutanopia")
# recolorize::plotColorPalette(tritan, main = "Tritanopia")

cites_cmn <- read.table("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Nature 2024/Datasets/cmn_papers_wos.txt", header=TRUE, sep="\t")

# Convert Publication Years to numeric
cites_cmn$Publication.Years <- as.numeric(cites_cmn$Publication.Years)

# Create bins for every 5 years
bin <- 10

max_year <- max(cites_cmn$Publication.Years)
# Adjust the min_year to ensure it is included in the bins
min_year <- min(cites_cmn$Publication.Years)

# Adjust the min_year to ensure it is included in the bins
adjusted_min_year <- min_year - (bin - (max_year - min_year) %% bin)

# Create the breaks
breaks <- seq(max_year, adjusted_min_year, by = -bin)

# Cut the data into bins
cites_cmn$Year.Bin <- cut(cites_cmn$Publication.Years, breaks=breaks, right=TRUE, include.lowest=TRUE, dig.lab = 50)


# Summarize the Record Counts by 5-year bins
binned_data <- aggregate(Record.Count ~ Year.Bin, cites_cmn, sum) 

theme_set(theme_bw())

ggplot(binned_data, aes(x = Year.Bin, y = Record.Count))+
  geom_col()+
 # theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Years (10-Year Bins)",
       y = "Number of Records (Web of Science)")

binned_data <- binned_data %>%
  mutate(
    Percent_Diff = (Record.Count - lag(Record.Count)) / lag(Record.Count) * 100
  )
#>90% increase between 2004-2014 and 2014-2024