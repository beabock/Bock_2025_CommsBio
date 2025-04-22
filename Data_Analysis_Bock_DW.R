#B. Bock
#21 June 2024
#This file has all data import, data analysis, and plots for the Dark Web experiment/paper. 
#4/17/25: Edits back from editor, making some figures for the main text regarding both rounds of experimentation

library(readr) #read_csv
library(readxl)
library(stringr) #read strings
library(ggplot2)
library(dplyr)
library(rempsyc)
library(tidyr) #separate function
library(EnvStats) #For adding sample sizes to plots
library(emmeans)
library(purrr)
#library(multcomp)

theme_set(theme_bw())

custom_col <- c("#D15A62", "#5AA7D1")

setwd("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web")



# Intro barplot -----------------------------------------------------------
#From Web of Science on 2/11/2025 by BB
#Searching just "Common Mycorrhizal Network"
# cmns <- read.delim("Comms Bio 2025/Datasets/pubs_CMN.txt", sep = "\t", header = TRUE)%>%
#   rename("CMNs" = Record.Count)
# 
# #Searching for "mycorrhizae OR mycorrhizal"
# myc <- read.delim("Comms Bio 2025/Datasets/pubs_mycorrhizae_mycorrhizal.txt", sep = "\t", header = TRUE)%>%
#   rename("Mycorrhizae_Mycorrhizal" = Record.Count)
# 
# both <- full_join(cmns, myc, by = "Publication.Years")%>%
#   select(!starts_with("X"))%>%
#   mutate(CMNs_Relativized = CMNs/Mycorrhizae_Mycorrhizal)%>%
#   replace(is.na(.), 0) 
# 
# both %>%
#   filter(Publication.Years != 2025,
#          Publication.Years >= 1950)%>%
#   pivot_longer(cols = -Publication.Years, names_to = "Type", values_to = "Records") %>%
#   ggplot(aes(x = Publication.Years, y = Records)) +
#   geom_col() +
#   geom_smooth(se = F)+
#   facet_wrap(~Type, scales = "free_y")

# Dataset Cleanup -----------------------------------------------------------------


#write.csv(combo_ds2, "Datasets/combo_ds2.csv")
#Uncomment the above line to export the dataset with all data combined


# Biomass Plot ------------------------------------------------------------


#Main plot

combo_ds2 <- read.csv("Datasets/biomass_ds_long.csv")

desired_order <- c("Experimental", "Impermeable", "Sterile")

plot <- combo_ds2 %>%
  filter(Experiment_Round != "Preliminary")%>%
  filter(!Type_Barrier %in% c("Barrierless", "Diffusion1", "Diffusion2", "Diffusion", "Diffusion3", "Oyster"))%>% 
  group_by(Experiment_Round, Type_Barrier, Box_Nr, Chamber)%>%
  summarize(Plant_Weight_g = Dry_Weight_g[Compartment == "Root"]+
              Dry_Weight_g[Compartment == "Shoot"])%>%
  mutate(Type_Barrier = factor(Type_Barrier, levels = desired_order),
         Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")))%>%
  as.data.frame()


#Normality
plot %>%
  ggplot(aes(x=Plant_Weight_g))+
  geom_histogram(bins = 5)
#Somewhat normal.


nested_models <- plot %>%
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Plant_Weight_g ~Type_Barrier:Chamber, data = .x)),
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
    emmeans_summary = map(emmeans, ~ summary(.x))
  )

all_emmeans <- nested_models %>%
  select(Experiment_Round, emmeans) %>%
  mutate(
    emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans)),
    contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x))))
  )

all_emmeans_long <- all_emmeans %>%
  unnest(emmeans_df)%>%
  select(-emmeans,-contrast_df)

all_contrasts_long <- all_emmeans %>%
  unnest(contrast_df)%>%
  select(-emmeans,-emmeans_df)

write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_biomass.csv", row.names = FALSE)
write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_biomass.csv", row.names = FALSE)


n_labels <- plot %>%
  group_by(Experiment_Round, Type_Barrier) %>%
  summarise(n_half = n() / 2, .groups = "drop") %>%
  mutate(label = paste0("n = ", n_half))

significant_contrasts <- all_contrasts_long %>%
  filter(p.value < 0.05) %>%
  mutate(
    Type_Barrier = factor(Type_Barrier, levels = c("Experimental", "Impermeable", "Sterile")),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*"
    )
  )


ggplot(all_emmeans_long, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = plot,
             aes(x = Type_Barrier, y = Plant_Weight_g, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 4, stroke = 1,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -0.05, label = label),
            color = "black", 
            size = 2.5, inherit.aes = FALSE)+
  geom_text(
    data = significant_contrasts,
    aes(x = Type_Barrier, y = 0.155, label = label),  # see note below
    color = "black",
    size = 4,
    position = position_dodge(width = 0.5),
    inherit.aes = FALSE
  )+
  scale_color_manual(values = custom_col, name = "Chamber")+
  labs(x = "Treatment",
       y = "Plant Biomass (g)"
  )+
  scale_x_discrete(labels = c(
    "Experimental" = "Permeable",
    "Impermeable" = "Impermeable",
    "Sterile" = "Axenic"
  )) 

ggsave(
  "Comms Bio 2025/Pub_Figures/biomass_main.png",
  width = 6.5, height = 5, units = "in", dpi = 600
)





#This one is for the supplementary.

desired_order <- c("Experimental", "Impermeable", "Sterile","Oyster", "Barrierless", "Diffusion1", "Diffusion2", "Diffusion3")

plot <- combo_ds2 %>%
  filter(!is.na(Type_Barrier))%>%
  filter(Experiment_Round != "Preliminary")%>%
  group_by(Experiment_Round, Type_Barrier, Box_Nr, Chamber)%>%
  summarize(Plant_Weight_g = Dry_Weight_g[Compartment == "Root"]+
              Dry_Weight_g[Compartment == "Shoot"])%>%
  mutate(Type_Barrier = factor(Type_Barrier, levels = desired_order),
         Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")))%>%
  as.data.frame()


#Normality
plot %>%
  ggplot(aes(x=Plant_Weight_g))+
  geom_histogram(bins = 5)
#Somewhat normal.


nested_models <- plot %>%
  filter(!is.na(Type_Barrier))%>%
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Plant_Weight_g ~Type_Barrier:Chamber, data = .x)),
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
    emmeans_summary = map(emmeans, ~ summary(.x))
  )

all_emmeans <- nested_models %>%
  select(Experiment_Round, emmeans) %>%
  mutate(
    emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans)),
    contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x))))
  )

all_emmeans_long <- all_emmeans %>%
  unnest(emmeans_df)%>%
  select(-emmeans,-contrast_df)

all_contrasts_long <- all_emmeans %>%
  unnest(contrast_df)%>%
  select(-emmeans,-emmeans_df)

write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_biomass_supp.csv", row.names = FALSE)
write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_biomass_supp.csv", row.names = FALSE)


n_labels <- plot %>%
  group_by(Experiment_Round, Type_Barrier) %>%
  summarise(n_half = n() / 2, .groups = "drop") %>%
  mutate(label = paste0("n = ", n_half))

significant_contrasts <- all_contrasts_long %>%
  filter(p.value < 0.05) %>%
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*"
    )
  )
  


ggplot(all_emmeans_long, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = plot,
             aes(x = Type_Barrier, y = Plant_Weight_g, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 4, stroke = 1,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -0.03, label = label),
            color = "black", 
            size = 2.5, inherit.aes = FALSE)+
  # geom_text(
  #   data = significant_contrasts,
  #   aes(x = Type_Barrier, y = 0.225, label = label),  # see note below
  #   color = "black",
  #   size = 4,
  #   position = position_dodge(width = 0.5),
  #   inherit.aes = FALSE
  # )+
  scale_color_manual(values = custom_col, name = "Chamber")+
  labs(x = "Treatment",
       y = "Plant Biomass (g)"
  )+
  scale_x_discrete(labels = c(
    "Experimental" = "Permeable",
    "Impermeable" = "Impermeable",
    "Sterile" = "Axenic"
  )) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "Comms Bio 2025/Pub_Figures/biomass_supp.png",
  width = 6.5, height = 5, units = "in", dpi = 600
)


#




# Spectral Data -----------------------------------------------------------


#Import the datasets needed to build the model we will use later. These were obtained by taking spec readings of the supernatant of samples with a known amount of dye in them. 

setwd("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/")

abs545 <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/abs_545.csv")%>%
  mutate(Experiment_Round = factor(Experiment_Round,
                            levels = c("Preliminary", "Main", "Follow-Up"))
  )

shoots <- abs545 %>%
  filter(Compartment == "Shoot")

roots <- abs545 %>%
  filter(Compartment == "Root")

root_control <- read.csv("Datasets/roots_mod.csv")%>%
  dplyr::select(!Weight_G)%>%
  rename(Weight_G = plant_g)

shoot_control <- read.csv("Datasets/shoots_mod.csv") %>%
  dplyr::select(!Weight_G)%>%
  rename(Weight_G = plant_g)

#Shoots model 

mod_s3 <- lm(data=shoot_control, dye_ug ~ abs)
summary(mod_s3)

mod_r3 <- lm(data = root_control, dye_ug ~ abs)
summary(mod_r3)


preds <- predict(mod_s3, newdata = shoots, se.fit = TRUE)


desired_order <- c("Experimental", "Impermeable", "Sterile","Oyster", "Barrierless", "Diffusion1", "Diffusion2", "Diffusion3")

shoots <- shoots %>%
  mutate(
    preds_mod_s3 = preds$fit,
    se_mod_s3 = preds$se.fit,
    Dry_Weight_ug = Dry_Weight_g * 1e6,
    dye_ug = preds_mod_s3 * Dry_Weight_ug,
    dye_ug_lower = (preds_mod_s3 - 1.96 * se_mod_s3) * Dry_Weight_ug,
    dye_ug_upper = (preds_mod_s3 + 1.96 * se_mod_s3) * Dry_Weight_ug
  )%>%
  filter(Experiment_Round != "Preliminary")%>% #Throws an error later if I don't do this
  mutate(Type_Barrier = factor(Type_Barrier, levels = desired_order))


preds_r <- predict(mod_r3, newdata = roots, se.fit = TRUE)

roots <- roots %>%
  mutate(
    preds_mod_r3 = preds_r$fit,
    se_mod_r3 = preds_r$se.fit,
    Dry_Weight_ug = Dry_Weight_g * 1e6,
    dye_ug = preds_mod_r3 * Dry_Weight_ug,
    dye_ug_lower = (preds_mod_r3 - 1.96 * se_mod_r3) * Dry_Weight_ug,
    dye_ug_upper = (preds_mod_r3 + 1.96 * se_mod_r3) * Dry_Weight_ug
  )%>%
  mutate(Type_Barrier = factor(Type_Barrier, levels = desired_order))%>%
  filter(Experiment_Round != "Preliminary") #Throws an error later if I don't do this


peaks <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/DW_NanoDrop_all.csv")



#Goes into supplementary
avg_abs <- peaks %>%
  mutate(Chamber = case_when(grepl("A", Sample_Name) ~ "Donor",
                             grepl("B", Sample_Name) ~ "Receiver", 
                             .default = NA
  ))%>%
#  filter(Type_Barrier != "Diffusion1" & Type_Barrier != "Diffusion2" & Type_Barrier != "Diffusion3")%>%
  ggplot(aes(x= waves, y =abs, color = Type_Barrier))+
  scale_color_discrete(labels = c("Experimental" = "Permeable", "Impermeable" = "Impermeable", "Sterile" = "Axenic"))+
  stat_summary(fun.y=mean, geom="line")+
  geom_vline(xintercept = 545)+
  guides(color = guide_legend(title = "Barrier Type")) +
  labs(x = "Wavelength (nm)", y = "Avg Absorbance")+
  facet_grid(Compartment ~ Chamber)

#450x500 or so

#avg_abs


#Data included in the dye plot, only comparing receiver plants to each other
plot <- shoots %>%
  filter(Experiment_Round != "Preliminary")%>% 
  filter(!Type_Barrier %in% c("Diffusion1", "Barrierless", "Oyster", "Diffusion2", "Diffusion3"))

plot_roots <- roots %>%
  filter(Experiment_Round != "Preliminary")%>% 
  filter(!Type_Barrier %in% c("Diffusion1", "Barrierless", "Oyster", "Diffusion2", "Diffusion3"))


#Testing if any differences among treatment in amount of dyed shoot material weighed out for this analysis. Answer: no detectable difference in Round 4. If including round 5, may need to include that in the model. 

anova(lm(Weight_G~Type_Barrier, plot)) #p > 0.05
iqr(plot$Weight_G) #0.0003425

anova(lm(Weight_G~Type_Barrier, plot_roots)) #p > 0.05
iqr(plot_roots$Weight_G) #0.000405, update in ms


 plot %>%
   group_by(Experiment_Round, Type_Barrier, Box_Nr)%>%
   summarize(n=n())

 
 plot_roots %>%
   group_by(Experiment_Round, Type_Barrier)%>%
   summarize(n=n()) #Not missing any samples.
 
#Main dye plot for text, shoots
 
 nested_models <- plot %>%
   filter(Experiment_Round != "Preliminary")%>% #Only one barrier type for prelim
   group_by(Experiment_Round) %>%
   nest() %>%
   mutate(
     model = map(data, ~ lm(preds_mod_s3 * Dry_Weight_ug ~ Type_Barrier:Chamber, data = .x)),
     emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
     emmeans_summary = map(emmeans, ~ summary(.x))
   )
 
all_emmeans <- nested_models %>%
   select(Experiment_Round, emmeans) %>%
   mutate(
     emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans)),
     contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x))))
   )
 
 all_emmeans_long <- all_emmeans %>%
   unnest(emmeans_df)%>%
   select(-emmeans,-contrast_df)
 
 all_contrasts_long <- all_emmeans %>%
   unnest(contrast_df)%>%
   select(-emmeans,-emmeans_df)
 
 write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_dye_shoots.csv", row.names = FALSE)
 write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_dye_shoots.csv", row.names = FALSE)
 
 
 n_labels <- plot %>%
   group_by(Experiment_Round, Type_Barrier) %>%
   summarise(n_half = n() / 2, .groups = "drop") %>%
   mutate(label = paste0("n = ", n_half))
 
 significant_contrasts <- all_contrasts_long %>%
   filter(p.value < 0.05) %>%
   mutate(
     label = case_when(
       p.value < 0.001 ~ "***",
       p.value < 0.01 ~ "**",
       p.value < 0.05 ~ "*"
     )
   )
 
 #Best dye plot
 
 ggplot(all_emmeans_long, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
   geom_hline(yintercept = 0, linetype = "dashed")+
   geom_point(data = plot,
               aes(x = Type_Barrier, y = preds_mod_s3*Dry_Weight_ug, color = Chamber), shape = 4,
               width = 0.15, alpha = .4, size = 4, stroke = 1,
              position = position_dodge(width = 0.5)) +
   geom_point(position = position_dodge(width = 0.5), size = 3) +
   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                 position = position_dodge(width = 0.5),
                 width = 0.2) +
   facet_grid(Experiment_Round~., scales = "free_y") +
   geom_text(data = n_labels,
             aes(x = Type_Barrier, y = -2000, label = label),
             color = "black", 
             size = 2.5, inherit.aes = FALSE)+
   geom_text(
     data = significant_contrasts,
     aes(x = Type_Barrier, y = 6500, label = label),  # see note below
     color = "black",
     size = 4,
     inherit.aes = FALSE
   )+
   scale_color_manual(values = custom_col, name = "Chamber")+
   labs(x = "Treatment",
     y = "Dye Content in Leaves (μg)"
   )+
   scale_x_discrete(labels = c(
     "Experimental" = "Permeable",
     "Impermeable" = "Impermeable",
     "Sterile" = "Axenic"
   )) 
 
 ggsave(
   "Comms Bio 2025/Pub_Figures/dye_leaves_main.png",
   width = 6.5, height = 5, units = "in", dpi = 600
 )
 
#Adjust the code here for saving the tables
 

 

#Now let's do the same for the supplementary, so not excluding any treatments

 plot <- shoots %>%
   filter(Experiment_Round != "Preliminary")

nested_models <-  shoots%>%
  filter(Experiment_Round != "Preliminary")%>% #Still can't include this because no contrats available
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(preds_mod_s3 * Dry_Weight_ug ~ Type_Barrier:Chamber, data = .x)),
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
    emmeans_summary = map(emmeans, ~ summary(.x))
  )


all_emmeans <- nested_models %>%
  select(Experiment_Round, emmeans) %>%
  mutate(
    emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans)),
    contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x))))
  )

all_emmeans_long <- all_emmeans %>%
  unnest(emmeans_df)%>%
  select(-emmeans,-contrast_df)

all_contrasts_long <- all_emmeans %>%
  unnest(contrast_df)%>%
  select(-emmeans,-emmeans_df)

write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_dye_shoots_supp.csv", row.names = FALSE)
write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_dye_shoots_supp.csv", row.names = FALSE)


n_labels <- shoots %>%
  group_by(Experiment_Round, Type_Barrier) %>%
  summarise(n_half = n() / 2, .groups = "drop") %>%
  mutate(label = paste0("n = ", n_half))

significant_contrasts <- all_contrasts_long %>%
  filter(p.value < 0.05) %>%
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*"
    )
  )

#Best dye plot

ggplot(all_emmeans_long, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = plot,
             aes(x = Type_Barrier, y = preds_mod_s3*Dry_Weight_ug, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 4, stroke = 1,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -2500, label = label),
            color = "black", 
            size = 2.5, inherit.aes = FALSE)+
  # geom_text(
  #   data = significant_contrasts,
  #   aes(x = Type_Barrier, y = 8000, label = label),  # see note below
  #   color = "black",
  #   size = 4,
  #   inherit.aes = FALSE
  # )+
  scale_color_manual(values = custom_col, name = "Chamber")+
  labs(x = "Treatment",
       y = "Dye Content in Leaves (μg)"
  )+
  scale_x_discrete(labels = c(
    "Experimental" = "Permeable",
    "Impermeable" = "Impermeable",
    "Sterile" = "Axenic"
  )) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "Comms Bio 2025/Pub_Figures/dye_leaves_supp.png",
  width = 6.5, height = 5, units = "in", dpi = 600
)

#Supp roots fig now

plot_roots <- roots %>%
  filter(Experiment_Round != "Preliminary")

nested_models <-  roots%>%
  filter(Experiment_Round != "Preliminary")%>% #Still can't include this because no contrats available
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(preds_mod_r3 * Dry_Weight_ug ~ Type_Barrier:Chamber, data = .x)),
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
    emmeans_summary = map(emmeans, ~ summary(.x))
  )


all_emmeans <- nested_models %>%
  select(Experiment_Round, emmeans) %>%
  mutate(
    emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans)),
    contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x))))
  )

all_emmeans_long <- all_emmeans %>%
  unnest(emmeans_df)%>%
  select(-emmeans,-contrast_df)

all_contrasts_long <- all_emmeans %>%
  unnest(contrast_df)%>%
  select(-emmeans,-emmeans_df)

write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_dye_roots_supp.csv", row.names = FALSE)
write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_dye_roots_supp.csv", row.names = FALSE)


n_labels <- roots %>%
  group_by(Experiment_Round, Type_Barrier) %>%
  summarise(n_half = n() / 2, .groups = "drop") %>%
  mutate(label = paste0("n = ", n_half))

significant_contrasts <- all_contrasts_long %>%
  filter(p.value < 0.05) %>%
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*"
    )
  )

#Best dye plot

ggplot(all_emmeans_long, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = plot_roots,
             aes(x = Type_Barrier, y = preds_mod_r3*Dry_Weight_ug, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 4, stroke = 1,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -500, label = label),
            color = "black", 
            size = 2.5, inherit.aes = FALSE)+
  # geom_text(
  #   data = significant_contrasts,
  #   aes(x = Type_Barrier, y = 1000, label = label),  # see note below
  #   color = "black",
  #   size = 4,
  #   inherit.aes = FALSE
  # )+
  scale_color_manual(values = custom_col, name = "Chamber")+
  labs(x = "Treatment",
       y = "Dye Content in Roots (μg)"
  )+
  scale_x_discrete(labels = c(
    "Experimental" = "Permeable",
    "Impermeable" = "Impermeable",
    "Sterile" = "Axenic"
  )) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave(
  "Comms Bio 2025/Pub_Figures/dye_roots_supp.png",
  width = 6.5, height = 5, units = "in", dpi = 600
)

