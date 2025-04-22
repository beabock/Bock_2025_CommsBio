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

#Supplementary biomass plot, with all data in it, separated out by aboveground and belowground parts of plants. 

#Note to self: Should I plot the emmeans like i do with the dye below?

combo_ds2 <- read.csv(ds_biomass_long)

supp_biomass <- combo_ds2 %>%
  mutate(Compartment = factor(Compartment, levels = c("Shoot", "Root"), 
                              labels = c("Aboveground", "Belowground")))%>%
  ggplot(aes(y = Dry_Weight_g, x= Type_Barrier, fill = Chamber))+
  geom_boxplot()+
  facet_grid(Compartment~Experiment_Round, scales = "free", space = "free")+
  stat_n_text()+
  geom_point(aes(fill = Chamber),stroke = 0.5, position = position_dodge(width = 0.4), size = 2, alpha = 0.7,shape = 21, color = "black")+
  scale_color_manual(values = custom_col, name = "Chamber")+
  scale_fill_manual(values = custom_col, name = "Chamber")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_discrete(labels = c("Experimental" = "Permeable", "Impermeable" = "Impermeable", "Sterile" = "Axenic", "Diffusion1" = "Diffusion"))+
  labs(
    x = "Treatment", 
    y = "Dry Mass (g)")+
  theme(text = element_text(size = 18))+
  stat_n_text()

supp_biomass

ggsave("Comms Bio 2025/Pub_Figures/supp_biomass.png", supp_biomass, dpi = 1000, width = 13, height = 6)

#Adding aboveground and belowground masses in this ds
plot <- combo_ds2 %>%
  filter(!Type_Barrier %in% c("Barrierless", "Diffusion1", "Diffusion2", "Diffusion", "Diffusion3", "Oyster"))%>% #Barrierless and Diffusion1 were treatments used as controls, but we did not make hypotheses based on these treatments, so they will go into the supplementary.
 # filter(Experiment_Round == 4)%>% #Round 4 is the experiment in the main text.
  group_by(Experiment_Round, Type_Barrier, Box_Nr, Chamber)%>%
  summarize(Plant_Weight_g = Dry_Weight_g[Compartment == "Root"]+
              Dry_Weight_g[Compartment == "Shoot"])%>%
  as.data.frame()


#Normality
plot %>%
  ggplot(aes(x=Plant_Weight_g))+
  geom_histogram(bins = 5)
#Somewhat normal.


#Stats for biomass plot in main text
mod <- lm(Plant_Weight_g ~ Chamber * Type_Barrier + Experiment_Round, data = plot)

emmeans_results <- emmeans(mod, pairwise ~ Chamber | Type_Barrier)

emmeans_results

emmeans_table <- as_tibble(summary(emmeans_results$emmeans)) 
contrasts_table <- as_tibble(summary(emmeans_results$contrasts))  

# Write to CSV files
write.csv(emmeans_table, "Comms Bio 2025/Pub_Figures/emmeans_table_biomass.csv", row.names = FALSE)
write.csv(contrasts_table, "Comms Bio 2025/Pub_Figures/contrasts_table_biomass.csv", row.names = FALSE)


significant_contrasts <- contrasts_table %>%
  filter(p.value < 0.05) %>%
  mutate(
    Type_Barrier = factor(Type_Barrier, levels = c("Experimental", "Impermeable", "Sterile")),
    label = ifelse(p.value < 0.01, "*", "*") # Significance labels
  )

#This one is for the supplementary, lumped by plant, not separated out by aboveground and belowground.
# plot_supp <- combo_ds2 %>%
#   group_by(Experiment_Round, Type_Barrier, Box_Nr, Chamber)%>%
#   summarize(Plant_Weight_g = Dry_Weight_g[Compartment == "Root"]+
#               Dry_Weight_g[Compartment == "Shoot"]) 


biomass <- plot %>%
  mutate(Type_Barrier = factor(Type_Barrier, levels = c("Experimental", "Impermeable", "Sterile")))%>%
  ggplot(aes(
    x= Type_Barrier,
    y = Plant_Weight_g,
    fill = Chamber
  ))+
  geom_boxplot()+
  geom_point(aes(fill = Chamber),stroke = 0.5, position = position_dodge(width = 0.4), size = 2, alpha = 0.7,shape = 21, color = "black")+
  scale_color_manual(values = custom_col, name = "Chamber")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_fill_manual(values = c("#D15A62", "#5AA7D1"), name = "Chamber")+
  scale_x_discrete(labels = c("Experimental" = "Permeable", "Impermeable" = "Impermeable", "Sterile" = "Axenic", "Diffusion1" = "Diffusion"))+
  labs(
    x = "Treatment", 
    y = "Total Plant Dry Mass (g)")+
  theme(#legend.title=element_blank(),
        #axis.title.x=element_blank(),
        text = element_text(size = 11))+
   stat_n_text()+
   facet_grid(~Experiment_Round)
  
biomass

ggsave("Comms Bio 2025/Pub_Figures/biomass.png", biomass, dpi = 600, width = 5, height = 4)




  #Not filtering for aboveground vs. belowground here
# biomass_supp2 <- plot_supp %>%
#   ggplot(aes(
#     x= Type_Barrier,
#     y = Plant_Weight_g,
#     fill = Chamber
#   ))+
#   geom_boxplot()+
#   scale_fill_manual(values = c("#D15A62", "#5AA7D1"))+
#   scale_x_discrete(labels = c("Experimental" = "Permeable", "Impermeable" = "Impermeable", "Sterile" = "Axenic"))+
#   labs(
#       x = "Treatment", 
#     y = "Plant Dry Mass (g)")+
#   theme(legend.title=element_blank(),
#         axis.title.x = element_blank(),
#         text = element_text(size = 18))+
#   stat_n_text()
  
# biomass_supp2





# Spectral Data -----------------------------------------------------------


#Import the datasets needed to build the model we will use later. These were obtained by taking spec readings of the supernatant of samples with a known amount of dye in them. 

abs545 <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/abs_545.csv")

shoots <- abs545 %>%
  filter(Compartment == "Shoot")

roots <- abs545 %>%
  filter(Compartment == "Root")

root_control <- read.csv("Comms Bio 2025/Datasets/roots_mod.csv")%>%
  dplyr::select(!Weight_G)%>%
  rename(Weight_G = plant_g)

shoot_control <- read.csv("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Datasets/shoots_mod.csv") %>%
  dplyr::select(!Weight_G)%>%
  rename(Weight_G = plant_g)

#Shoots model 

mod_s3 <- lm(data=shoot_control, dye_ug ~ abs)
summary(mod_s3)

mod_r3 <- lm(data = root_control, dye_ug ~ abs)
summary(mod_r3)

setwd("Comms Bio 2025")



preds <- predict(mod_s3, newdata = shoots, se.fit = TRUE)

shoots <- shoots %>%
  mutate(
    preds_mod_s3 = preds$fit,
    se_mod_s3 = preds$se.fit,
    Dry_Weight_ug = Dry_Weight_g * 1e6,
    dye_ug = preds_mod_s3 * Dry_Weight_ug,
    dye_ug_lower = (preds_mod_s3 - 1.96 * se_mod_s3) * Dry_Weight_ug,
    dye_ug_upper = (preds_mod_s3 + 1.96 * se_mod_s3) * Dry_Weight_ug
  )

roots <- roots %>%
  mutate(preds_mod_r3 = predict(mod_r3, newdata = roots))




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

avg_abs




#Data included in the dye plot, only comparing receiver plants to each other
plot <- shoots %>%
 # filter(Experiment_Round == 4)%>% 
  filter(!Type_Barrier %in% c("Diffusion1", "Barrierless", "Oyster", "Diffusion2", "Diffusion3"))

plot_roots <- roots %>%
 # filter(Experiment_Round == 4)%>% 
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
 

 
 nested_models <- plot %>%
   group_by(Experiment_Round) %>%
   nest() %>%
   mutate(
     model = map(data, ~ lm(preds_mod_s3 * Dry_Weight_ug ~ Type_Barrier:Chamber, data = .x)),
     emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
     emmeans_summary = map(emmeans, ~ summary(.x))
   )
 
 all_emmeans <- nested_models %>%
   select(Experiment_Round, emmeans) %>%
   mutate(emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans))) %>%
   unnest(emmeans_df)
 
 all_emmeans
 #Remember to save the above to include output in manuscript
 
 n_labels <- plot %>%
   group_by(Experiment_Round, Type_Barrier) %>%
   summarise(n_half = n() / 2, .groups = "drop") %>%
   mutate(label = paste0("n = ", n_half))
 
 
 #Best dye plot
 
 ggplot(all_emmeans, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
   geom_hline(yintercept = 0, linetype = "dashed")+
   geom_point(data = plot,
               aes(x = Type_Barrier, y = preds_mod_s3*Dry_Weight_ug, color = Chamber), shape = 4,
               width = 0.15, alpha = .4, size = 4, stroke = 1,
              position = position_dodge(width = 0.5)) +
   geom_point(position = position_dodge(width = 0.5), size = 3) +
   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                 position = position_dodge(width = 0.5),
                 width = 0.2) +
   facet_wrap(~ Experiment_Round) +
   geom_text(data = n_labels,
             aes(x = Type_Barrier, y = 0, label = label),
             color = "black",  # or a neutral color since it's not per chamber
             vjust = 9, size = 3, inherit.aes = FALSE)+
   scale_color_manual(values = custom_col, name = "Chamber")+
   labs(
     y = "Dye Content in Leaves (μg)",
     title = "EMMeans by Chamber and Type_Barrier, Faceted by Experiment Round"
   )+
   scale_x_discrete(labels = c(
     "Experimental" = "Permeable",
     "Impermeable" = "Impermeable",
     "Sterile" = "Axenic"
   )) 
 
 
#Adjust the code here for saving the tables
 
#Now same thing but With the roots
 nested_models <- plot_roots %>%
   group_by(Experiment_Round) %>%
   nest() %>%
   mutate(
     model = map(data, ~ lm(preds_mod_r3 * Dry_Weight_ug ~ Type_Barrier:Chamber, data = .x)),
     emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
     emmeans_summary = map(emmeans, ~ summary(.x))
   )
 
 all_emmeans <- nested_models %>%
   select(Experiment_Round, emmeans) %>%
   mutate(emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans))) %>%
   unnest(emmeans_df)
 
 all_emmeans
 
 
 n_labels <- plot_roots %>%
   group_by(Experiment_Round, Type_Barrier) %>%
   summarise(n_half = n() / 2, .groups = "drop") %>%
   mutate(label = paste0("n = ", n_half))
 
 
 #Root dye plot
 
 ggplot(all_emmeans, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
   geom_hline(yintercept = 0, linetype = "dashed")+
   geom_point(data = plot_roots,
              aes(x = Type_Barrier, y = preds_mod_r3*Dry_Weight_ug, color = Chamber), shape = 4,
              width = 0.15, alpha = .4, size = 4, stroke = 1,
              position = position_dodge(width = 0.5)) +
   geom_point(position = position_dodge(width = 0.5), size = 3) +
   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                 position = position_dodge(width = 0.5),
                 width = 0.2) +
   facet_wrap(~ Experiment_Round) +
   geom_text(data = n_labels,
             aes(x = Type_Barrier, y = 0, label = label),
             color = "black",  # or a neutral color since it's not per chamber
             vjust = 9, size = 3, inherit.aes = FALSE)+
   scale_color_manual(values = custom_col, name = "Chamber")+
   labs(
     y = "Dye Content in Roots (μg)",
     title = "EMMeans by Chamber and Type_Barrier, Faceted by Experiment Round"
   )+
   scale_x_discrete(labels = c(
     "Experimental" = "Permeable",
     "Impermeable" = "Impermeable",
     "Sterile" = "Axenic"
   ))  
 #Roots are super variable, not much to pull from this.
 
 # Write to CSV files
 write.csv(emmeans_table, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_dye.csv", row.names = FALSE)
 write.csv(contrasts_table, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_dye.csv", row.names = FALSE)
 
 write.csv(emmeans_table_r, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_dye_roots.csv", row.names = FALSE)
 write.csv(contrasts_table_r, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_dye_roots.csv", row.names = FALSE)

 
 # Add significant contrasts
 significant_contrasts_dye <- contrasts_table %>%
   filter(p.value < 0.05) %>%
   mutate(label = case_when(
     p.value < 0.001 ~ "*",
     p.value < 0.01 ~ "*",
     p.value < 0.05 ~ "*"
   ))
 #Maybe come back and add the above onto the plots

 

#850 x 500
 
 

setwd("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures")

ggsave("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/dye.png", dye, dpi = 1000, width = 5.5, height = 3.5)

ggsave("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/dye_roots.png", dye_r, dpi = 1000, width = 5.5, height = 3.5)



nested_models <-  shoots%>%
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(preds_mod_s3 * Dry_Weight_ug ~ Type_Barrier:Chamber, data = .x)),
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
    emmeans_summary = map(emmeans, ~ summary(.x))
  )


n_labels <- shoots %>%
  group_by(Experiment_Round, Type_Barrier) %>%
  summarise(n_half = n() / 2, .groups = "drop") %>%
  mutate(label = paste0("n = ", n_half))

all_emmeans <- nested_models %>%
  select(Experiment_Round, emmeans) %>%
  mutate(emmeans_df = map(emmeans, ~ as_tibble(summary(.x)$emmeans))) %>%
  unnest(emmeans_df)

dye_supp <- ggplot(all_emmeans, aes(x = Type_Barrier, y = emmean, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = shoots,
             aes(x = Type_Barrier, y = preds_mod_s3*Dry_Weight_ug, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 4, stroke = 1,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_wrap(~ Experiment_Round, labeller = labeller(Experiment_Round = function(x) paste("Experiment Round:", x)))+
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = 0, label = label),
            color = "black",  # or a neutral color since it's not per chamber
            vjust = 10.5, size = 3, inherit.aes = FALSE)+
  scale_color_manual(values = custom_col, name = "Chamber")+
  labs(
    y = "Dye Content in Leaves (μg)",
    title = "EMMeans by Chamber and Type_Barrier, Faceted by Experiment Round"
  )+
  scale_x_discrete(labels = c(
    "Experimental" = "Permeable",
    "Impermeable" = "Impermeable",
    "Sterile" = "Axenic"
  )) 

dye_supp



dye_supp_r <- roots %>%
  mutate(Dry_Weight_ug = Dry_Weight_g * 1e6)%>%
  ggplot(aes(y=preds_mod_r3*Dry_Weight_g, x = Type_Barrier, fill = Chamber))+
  geom_hline(yintercept = 0,  linetype = "dashed")+
  geom_boxplot()+
  geom_point(aes(fill = Chamber),stroke = 0.5, position = position_dodge(width = 0.4), size = 2, alpha = 0.7,shape = 21, color = "black")+
  scale_fill_manual(labels = c("Donor", "Receiver"), values = c("#D15A62", "#5AA7D1"),)+
  labs(y = "Dye in roots (ug)", x = "Treatment")+
  scale_x_discrete(labels = c("Experimental" = "Permeable", "Impermeable" = "Impermeable", "Sterile" = "Axenic"))+
  theme(legend.title=element_blank(),
        #axis.title.x=element_blank(),
        text = element_text(size = 18))+
  stat_n_text()


dye_supp_r

ggsave("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/dye_supp.png", dye_supp, dpi = 1000, width = 14, height = 8)

ggsave("C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/dye_supp_r.png", dye_supp_r, dpi = 1000, width = 14, height = 8)

