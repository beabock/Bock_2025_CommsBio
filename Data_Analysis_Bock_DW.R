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
library(broom)
library(lmtest)
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

desired_order <- c("Permeable", "Impermeable", "Axenic")

plot <- combo_ds2 %>%
  filter(Experiment_Round != "Preliminary")%>%
  filter(!Type_Barrier %in% c("Barrierless", "Diffusion1", "Diffusion2", "Diffusion", "Diffusion3", "Oyster"))%>% 
  mutate(
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  )%>%
  group_by(Experiment_Round, Type_Barrier, Box_Nr, Chamber)%>%
  summarize(Plant_Weight_g = Dry_Weight_g[Compartment == "Root"]+
              Dry_Weight_g[Compartment == "Shoot"])%>%
  mutate(Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")))%>%
  as.data.frame()


#Normality
plot %>%
  ggplot(aes(x=Plant_Weight_g))+
  geom_histogram(bins = 5)
#Somewhat normal.


nested_models_log <- plot %>%
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(Plant_Weight_g ~ Type_Barrier:Chamber, family = quasipoisson(link = "log"), data = .x)),
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
    
    # Back-transformed emmeans and contrasts
    emmeans_backtrans = map(emmeans, ~ summary(.x, type = "response")),
    emmeans_summary_df = map(emmeans_backtrans, ~ as_tibble(.x$emmeans)),
    contrast_summary_df = map(emmeans_backtrans, ~ as_tibble(.x$contrasts)),
    
    # Optional: diagnostics
    model_aug = map(model, ~ broom::augment(.x)),
    residuals = map(model_aug, ~ .x$.resid),
    fitted = map(model_aug, ~ .x$.fitted),
    
    shapiro_test = map(residuals, ~ shapiro.test(.x)),
    
    # Breusch-Pagan test for heteroscedasticity
    bp_test = map(model, ~ bptest(.x)),
    
    # QQ plot for residuals
    qq_plot = map2(residuals, Experiment_Round, ~
                     ggplot(data.frame(resid = .x), aes(sample = resid)) +
                     stat_qq() +
                     stat_qq_line() +
                     ggtitle(paste("QQ Plot (log) -", .y))
    ),
    resid_fit_plot = map2(residuals, fitted, ~
                            ggplot(data.frame(fitted = .y, resid = .x), aes(x = fitted, y = resid)) +
                            geom_point(alpha = 0.6) +
                            geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                            ggtitle("Residuals vs Fitted") +
                            theme_minimal()
    )
  )


nested_models_log %>%
  pull(qq_plot) %>%
  walk(print)

nested_models_log %>%
  pull(resid_fit_plot) %>%
  walk(print)
#Looks fine.

nested_models_log$bp_test

nested_models_log$shapiro_test
#Doesnt pass, but quasipoisson doesnt require normality so we're okay.

all_emmeans <- nested_models_log %>%
  select(Experiment_Round, emmeans) %>%
  mutate(
    emmeans_df = map(emmeans, ~ as_tibble(summary(.x, type = "response")$emmeans)),
    contrast_df = map(emmeans, ~ as_tibble(summary(.x, type = "response")$contrasts))
  )

all_emmeans_long <- all_emmeans %>%
  dplyr::select(Experiment_Round, emmeans_df) %>%
  unnest(emmeans_df) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  ) %>%
  select(-df)

all_contrasts_long <- all_emmeans %>%
  dplyr::select(Experiment_Round, contrast_df) %>%
  unnest(contrast_df) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
    p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "f", digits = 3)),
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  ) %>%
  select(-df)

write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_biomass.csv", row.names = FALSE)
write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_biomass.csv", row.names = FALSE)


n_labels <- plot %>%
  group_by(Experiment_Round, Type_Barrier) %>%
  summarise(n_half = n() / 2, .groups = "drop") %>%
  mutate(label = paste0("n = ", n_half))

significant_contrasts <- all_contrasts_long %>%
  filter(p.value < 0.05) %>%
  mutate(
    Type_Barrier = factor(Type_Barrier, levels = c("Permeable", "Impermeable", "Axenic")),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*"
    )
  )


ggplot(all_emmeans_long, aes(x = Type_Barrier, y = rate, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = plot,
             aes(x = Type_Barrier, y = Plant_Weight_g, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 7.5, stroke = 1.2,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -0.05, label = label),
            color = "black", 
           inherit.aes = FALSE)+
  geom_text(
    data = significant_contrasts,
    aes(x = Type_Barrier, y = 0.155, label = label),  # see note below
    color = "black",
    size = 6,
    position = position_dodge(width = 0.5),
    inherit.aes = FALSE
  )+
  scale_color_manual(values = custom_col, name = "Chamber")+
  labs(x = "Treatment",
       y = "Plant Dry Biomass (g)"
  )+
  scale_x_discrete(labels = c(
    "Experimental" = "Permeable",
    "Impermeable" = "Impermeable",
    "Sterile" = "Axenic"
  )) +
  theme(
    text = element_text(size = 16))

ggsave(
  "Comms Bio 2025/Pub_Figures/biomass_main.png",
  width = 6.5, height = 5, units = "in", dpi = 600
)





#This one is for the supplementary.

desired_order <- c("Permeable", "Impermeable", "Axenic","Oyster", "Barrierless", "Diffusion1", "Diffusion2", "Diffusion3")

plot <- combo_ds2 %>%
  mutate(
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  )%>%
  filter(!is.na(Type_Barrier))%>%
  filter(Experiment_Round != "Preliminary")%>%
  group_by(Experiment_Round, Type_Barrier, Box_Nr, Chamber)%>%
  summarize(Plant_Weight_g = Dry_Weight_g[Compartment == "Root"]+
              Dry_Weight_g[Compartment == "Shoot"])%>%
  mutate(Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")))%>%
  as.data.frame()


#Normality
plot %>%
  ggplot(aes(x=Plant_Weight_g))+
  geom_histogram(bins = 5)
#Somewhat normal.


nested_models_log <- plot %>%
  filter(!is.na(Type_Barrier)) %>%
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    # Log-transformed model
    model = map(data, ~ glm(Plant_Weight_g ~ Type_Barrier:Chamber, family = quasipoisson(link = "log"), data = .x)),
    
    # emmeans still works normally
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
    
    # Back-transformed emmeans and contrasts
    emmeans_backtrans = map(emmeans, ~ summary(.x, type = "response")),
    emmeans_summary_df = map(emmeans_backtrans, ~ as_tibble(.x$emmeans)),
    contrast_summary_df = map(emmeans_backtrans, ~ as_tibble(.x$contrasts)),
    
    # Residual diagnostics
    model_aug = map(model, ~ broom::augment(.x)),
    residuals = map(model_aug, ~ .x$.resid),
    fitted = map(model_aug, ~ .x$.fitted),
    
    # Shapiro-Wilk test for normality
    shapiro_test = map(residuals, ~ shapiro.test(.x)),
    
    bp_test = map(model, ~ bptest(.x)),
    # QQ plot for residuals
    qq_plot = map2(residuals, Experiment_Round, ~
                     ggplot(data.frame(resid = .x), aes(sample = resid)) +
                     stat_qq() +
                     stat_qq_line() +
                     ggtitle(paste("QQ Plot (log) -", .y))
    ),
    
    # Residuals vs Fitted plot
    resid_fit_plot = map2(residuals, fitted, ~
                            ggplot(data.frame(fitted = .y, resid = .x), aes(x = fitted, y = resid)) +
                            geom_point(alpha = 0.6) +
                            geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                            ggtitle("Residuals vs Fitted")
    )
  )


nested_models_log %>%
  pull(qq_plot) %>%
  walk(print)

nested_models_log %>%
  pull(resid_fit_plot) %>%
  walk(print)

nested_models_log %>%
  transmute(Experiment_Round,
            shapiro_p = map_dbl(shapiro_test, ~ .x$p.value),
            bp_p = map_dbl(bp_test, ~ .x$p.value))
#Passes.

all_emmeans <- nested_models_log %>%
  select(Experiment_Round, emmeans) %>%
  mutate(
    emmeans_df = map(emmeans, ~ as_tibble(summary(.x, type = "response")$emmeans)),
    contrast_df = map(emmeans, ~ as_tibble(summary(.x, type = "response")$contrasts))
  )

all_emmeans_long <- all_emmeans %>%
  dplyr::select(Experiment_Round, emmeans_df) %>%
  unnest(emmeans_df) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  ) %>%
  select(-df)

all_contrasts_long <- all_emmeans %>%
  dplyr::select(Experiment_Round, contrast_df) %>%
  unnest(contrast_df) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
    p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "f", digits = 3)),
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  ) %>%
  select(-df)

#write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_biomass_supp.csv", row.names = FALSE)
#write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_biomass_supp.csv", row.names = FALSE)


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
  


ggplot(all_emmeans_long, aes(x = Type_Barrier, y = rate, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = plot,
             aes(x = Type_Barrier, y = Plant_Weight_g, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 7.5, stroke = 1.2,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -0.03, label = label),
            color = "black", 
            inherit.aes = FALSE)+
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
  )+
  theme(
    text = element_text(size = 16))

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

hist(shoot_control$dye_ug)

length(shoot_control$samps)

mod_s3 <- lm(data=shoot_control, dye_ug ~ abs)
summary(mod_s3)

bptest(mod_s3) #Fine.

shapiro.test(mod_s3$residuals) #Does not pass.

hist(mod_s3$residuals)
ggplot(data.frame(residuals = residuals(mod_s3)), aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line()
#Actually fine. Shapiro-wilk p = 0.017, but the plot looks good.



ggsave("Comms Bio 2025/Pub_Figures/qqplot_shoot_controlmodel.png", width = 3, height = 3, dpi = 600)


mod_r3 <- lm(data = root_control, dye_ug ~ abs)
summary(mod_r3)

ggplot(data.frame(residuals = residuals(mod_r3)), aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line()
#Not awesome.

ggsave("Comms Bio 2025/Pub_Figures/qqplot_root_controlmodel.png", width = 3, height = 3, dpi = 600)


desired_order <- c("Permeable", "Impermeable", "Axenic","Oyster", "Barrierless", "Diffusion1", "Diffusion2", "Diffusion3")

preds <- predict(mod_s3, newdata = shoots, se.fit = TRUE)

#Something is getting messed up here. Come back to this

# Step 4: Continue with the pipeline for calculating dye_ug and confidence intervals
shoots <-  shoots %>% mutate(
  preds_mod_s3 = preds$fit,
  se_mod_s3 = preds$se.fit,
  Dry_Weight_ug = Dry_Weight_g * 1e6,
  dye_ug = preds_mod_s3 * Dry_Weight_ug,
  dye_ug_lower = (preds_mod_s3 - 1.96 * se_mod_s3) * Dry_Weight_ug,
  dye_ug_upper = (preds_mod_s3 + 1.96 * se_mod_s3) * Dry_Weight_ug
)%>%
  filter(Experiment_Round != "Preliminary")%>% #Throws an error later if I don't do this
  mutate(
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  )


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
  filter(!Type_Barrier %in% c("Diffusion1", "Barrierless", "Oyster", "Diffusion2", "Diffusion3"))%>%
  mutate(
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  )

plot_roots <- roots %>%
  filter(Experiment_Round != "Preliminary")%>% 
  filter(!Type_Barrier %in% c("Diffusion1", "Barrierless", "Oyster", "Diffusion2", "Diffusion3"))%>%
  mutate(
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  )


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
 
 ggplot(plot, aes(x = preds_mod_s3))+
   geom_histogram()+
   facet_grid(~Experiment_Round)
 

 nested_models <- plot %>%
   group_by(Experiment_Round) %>%
   nest() %>%
   mutate(
    model = map(data, ~ glm(preds_mod_s3 * Dry_Weight_ug ~ Type_Barrier:Chamber,
                                                 family = quasipoisson(link = "identity"), data = .x)),
     emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier)),
     emmeans_summary = map(emmeans, ~ summary(.x, type = "response")),
     
     # Extract model diagnostics
     model_aug = map(model, ~ augment(.x)),
     residuals = map(model_aug, ~ .x$.resid),
     fitted = map(model_aug, ~ .x$.fitted),
     
     shapiro_test = map(residuals, ~ shapiro.test(.x)),
     
     # Breusch-Pagan test for heteroscedasticity
     bp_test = map(model, ~ bptest(.x)),
     
     # QQ plot for residuals
     qq_plot = map2(residuals, Experiment_Round, ~
                      ggplot(data.frame(resid = .x), aes(sample = resid)) +
                      stat_qq() +
                      stat_qq_line() +
                      ggtitle(paste("QQ Plot -", .y))
     ),
     resid_fit_plot = map2(residuals, fitted, ~
                             ggplot(data.frame(fitted = .y, resid = .x), aes(x = fitted, y = resid)) +
                             geom_point(alpha = 0.6) +
                             geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                             ggtitle("Residuals vs Fitted") +
                             theme_minimal()
     )
   )
 
 
 nested_models %>%
   pull(qq_plot) %>%
   walk(print)
 
 nested_models %>%
   pull(resid_fit_plot) %>%
   walk(print)
   
 nested_models$shapiro_test
 
 nested_models$bp_test
 #Passes both.
 
 

 all_emmeans <- nested_models %>%
  # select(Experiment_Round, emmeans) %>%
   mutate(
     # Extracting EMMs and contrasts for Gamma GLM with back-transformation to response scale
     emmeans_df = map(emmeans, ~ as_tibble(summary(.x, type = "response")$emmeans)),
     contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x, type = "response"))))
   )
 
 
 all_emmeans_long <- all_emmeans %>%
   dplyr::select(Experiment_Round, emmeans_df) %>%
   unnest(emmeans_df) %>%
   mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
   mutate(
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
     Type_Barrier = case_when(
       Type_Barrier == "Experimental" ~ "Permeable",
       Type_Barrier == "Sterile" ~ "Axenic",
       TRUE ~ as.character(Type_Barrier)
     ),
     Type_Barrier = factor(Type_Barrier, levels = desired_order),
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")
     )) %>%
   select(-df)
 
 all_contrasts_long <- all_emmeans %>%
   dplyr::select(Experiment_Round, contrast_df) %>%
   unnest(contrast_df) %>%
   mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
   mutate(
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
     p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "f", digits = 3)),
     Type_Barrier = case_when(
       Type_Barrier == "Experimental" ~ "Permeable",
       Type_Barrier == "Sterile" ~ "Axenic",
       TRUE ~ as.character(Type_Barrier)
     ),
     Type_Barrier = factor(Type_Barrier, levels = desired_order),
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")
   )) %>%
   select(-df)
 
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
               width = 0.15, alpha = .4, size = 7.5, stroke = 1.2,
              position = position_dodge(width = 0.5)) +
   geom_point(position = position_dodge(width = 0.5), size = 5) +
   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 position = position_dodge(width = 0.5),
                 width = 0.2) +
   facet_grid(Experiment_Round~., scales = "free_y") +
   geom_text(data = n_labels,
             aes(x = Type_Barrier, y = -2000, label = label),
             color = "black",
             inherit.aes = FALSE)+
   geom_text(
     data = significant_contrasts,
     aes(x = Type_Barrier, y = 6500, label = label),  # see note below
     color = "black",
     size = 6,
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
   )) +
   theme(
     text = element_text(size = 16))
 
 ggsave(
   "Comms Bio 2025/Pub_Figures/dye_leaves_main.png",
   width = 6.5, height = 5, units = "in", dpi = 600
 )
 


 

#Now let's do the same for the supplementary, so not excluding any treatments

 plot <- shoots %>%
   filter(Experiment_Round != "Preliminary")
 
 nested_models <- plot %>%
   group_by(Experiment_Round) %>%
   nest() %>%
   mutate(
     # GLM model with quasipoisson distribution
     model = map(data, ~ glm((preds_mod_s3 * Dry_Weight_ug) ~ Type_Barrier:Chamber,
                             family = quasipoisson(link = "identity"), data = .x)),
     
     # emmeans with bootstrapped confidence intervals
     emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier, 
                                    adjust = "none", boot = TRUE, 
                                    reps = 1000)),  # Specify the number of bootstrap iterations
     
     # Summary of emmeans with back-transformation
     emmeans_summary = map(emmeans, ~ summary(.x, type = "response")),
     
     # Extract model diagnostics
     model_aug = map(model, ~ augment(.x)),
     residuals = map(model_aug, ~ .x$.resid),
     fitted = map(model_aug, ~ .x$.fitted),
     
     # Shapiro-Wilk test for normality of residuals
     shapiro_test = map(residuals, ~ shapiro.test(.x)),
     
     # Breusch-Pagan test for heteroscedasticity
     bp_test = map(model, ~ bptest(.x)),
     
     # QQ plot for residuals
     qq_plot = map2(residuals, Experiment_Round, ~
                      ggplot(data.frame(resid = .x), aes(sample = resid)) +
                      stat_qq() +
                      stat_qq_line() +
                      ggtitle(paste("QQ Plot -", .y))
     ),
     
     # Residuals vs Fitted plot for heteroscedasticity
     resid_fit_plot = map2(residuals, fitted, ~
                             ggplot(data.frame(fitted = .y, resid = .x), aes(x = fitted, y = resid)) +
                             geom_point(alpha = 0.6) +
                             geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                             ggtitle("Residuals vs Fitted") +
                             theme_minimal()
     )
   )
 nested_models$shapiro_test 
 nested_models$bp_test
 
 #Both not perfect. But also not trying to extract meaningful values for analysis. 
 
 nested_models %>%
   pull(qq_plot) %>%
   walk(print)
 
 nested_models %>%
   pull(resid_fit_plot) %>%
   walk(print)
 #Not great.
 
 #Bootstrapping because of the violations.
 
 all_emmeans <- nested_models %>%
  # select(Experiment_Round, emmeans) %>%
   mutate(
     # Extracting EMMs and contrasts for Gamma GLM with back-transformation to response scale
     emmeans_df = map(emmeans, ~ as_tibble(summary(.x, type = "response")$emmeans)),
     contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x, type = "response"))))
   )

 
 all_emmeans_long <- all_emmeans %>%
   dplyr::select(Experiment_Round, emmeans_df) %>%
   unnest(emmeans_df) %>%
   mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
   mutate(
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
     Type_Barrier = case_when(
       Type_Barrier == "Experimental" ~ "Permeable",
       Type_Barrier == "Sterile" ~ "Axenic",
       TRUE ~ as.character(Type_Barrier)
     ),
     Type_Barrier = factor(Type_Barrier, levels = desired_order),
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")
     )) %>%
   select(-df)
 
 all_contrasts_long <- all_emmeans %>%
   dplyr::select(Experiment_Round, contrast_df) %>%
   unnest(contrast_df) %>%
   mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
   mutate(
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
     p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "f", digits = 3)),
     Type_Barrier = case_when(
       Type_Barrier == "Experimental" ~ "Permeable",
       Type_Barrier == "Sterile" ~ "Axenic",
       TRUE ~ as.character(Type_Barrier)
     ),
     Type_Barrier = factor(Type_Barrier, levels = desired_order),
     Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")
     )) %>%
   select(-df)
 

#write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_dye_shoots_supp.csv", row.names = FALSE)
#write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_dye_shoots_supp.csv", row.names = FALSE)


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

#Dye shoots supp

ggplot(all_emmeans_long, aes(x = Type_Barrier, y = response, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = plot,
             aes(x = Type_Barrier, y = preds_mod_s3*Dry_Weight_ug, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 7.5, stroke = 1.2,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -2500, label = label),
            color = "black", 
            inherit.aes = FALSE)+
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
  )+
  theme(
    text = element_text(size = 16))

ggsave(
  "Comms Bio 2025/Pub_Figures/dye_leaves_supp.png",
  width = 6.5, height = 5, units = "in", dpi = 600
)

#Supp roots fig now

preds <- predict(mod_r3, newdata = roots, se.fit = TRUE)


#Something is getting messed up here. Come back to this

# Step 4: Continue with the pipeline for calculating dye_ug and confidence intervals
roots <-  roots %>% mutate(
  preds_mod_r3 = preds$fit,
  se_mod_r3 = preds$se.fit,
  Dry_Weight_ug = Dry_Weight_g * 1e6,
  dye_ug = preds_mod_r3 * Dry_Weight_ug,
  dye_ug_lower = (preds_mod_r3 - 1.96 * se_mod_r3) * Dry_Weight_ug,
  dye_ug_upper = (preds_mod_r3 + 1.96 * se_mod_r3) * Dry_Weight_ug
)%>%
  filter(Experiment_Round != "Preliminary")%>% #Throws an error later if I don't do this
  mutate(
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order)
  )


nested_models <- roots%>%
  filter(preds_mod_r3*Dry_Weight_ug>=0)%>% #Mod cant deal with negs
  group_by(Experiment_Round) %>%
  nest() %>%
  mutate(
    # GLM model with quasipoisson distribution
    model = map(data, ~ glm((preds_mod_r3 * Dry_Weight_ug) ~ Type_Barrier:Chamber,
                            family = quasipoisson(link = "identity"), data = .x)),
    
    # emmeans with bootstrapped confidence intervals
    emmeans = map(model, ~ emmeans(.x, pairwise ~ Chamber | Type_Barrier, 
                                   adjust = "none", boot = TRUE, 
                                   reps = 1000)),  # Specify the number of bootstrap iterations
    
    # Summary of emmeans with back-transformation
    emmeans_summary = map(emmeans, ~ summary(.x, type = "response")),
    
    # Extract model diagnostics
    model_aug = map(model, ~ augment(.x)),
    residuals = map(model_aug, ~ .x$.resid),
    fitted = map(model_aug, ~ .x$.fitted),
    
    # Shapiro-Wilk test for normality of residuals
    shapiro_test = map(residuals, ~ shapiro.test(.x)),
    
    # Breusch-Pagan test for heteroscedasticity
    bp_test = map(model, ~ bptest(.x)),
    
    # QQ plot for residuals
    qq_plot = map2(residuals, Experiment_Round, ~
                     ggplot(data.frame(resid = .x), aes(sample = resid)) +
                     stat_qq() +
                     stat_qq_line() +
                     ggtitle(paste("QQ Plot -", .y))
    ),
    
    # Residuals vs Fitted plot for heteroscedasticity
    resid_fit_plot = map2(residuals, fitted, ~
                            ggplot(data.frame(fitted = .y, resid = .x), aes(x = fitted, y = resid)) +
                            geom_point(alpha = 0.6) +
                            geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                            ggtitle("Residuals vs Fitted") +
                            theme_minimal()
    )
  )

nested_models %>%
  pull(qq_plot) %>%
  walk(print)

nested_models %>%
  pull(resid_fit_plot) %>%
  walk(print)

nested_models$shapiro_test
nested_models$bp_test

#Bootstrapped bc of slight violation

all_emmeans <- nested_models %>%
  dplyr::select(Experiment_Round, emmeans) %>%
  mutate(
    # Extracting EMMs and contrasts for Gamma GLM with back-transformation to response scale
    emmeans_df = map(emmeans, ~ as_tibble(summary(.x, type = "response")$emmeans)),
    contrast_df = map(emmeans, ~ as_tibble(summary(contrast(.x, type = "response"))))
  )


all_emmeans_long <- all_emmeans %>%
  dplyr::select(Experiment_Round, emmeans_df) %>%
  unnest(emmeans_df) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order),
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")
    )) %>%
  select(-df)

all_contrasts_long <- all_emmeans %>%
  dplyr::select(Experiment_Round, contrast_df) %>%
  unnest(contrast_df) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")),
    p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "f", digits = 3)),
    Type_Barrier = case_when(
      Type_Barrier == "Experimental" ~ "Permeable",
      Type_Barrier == "Sterile" ~ "Axenic",
      TRUE ~ as.character(Type_Barrier)
    ),
    Type_Barrier = factor(Type_Barrier, levels = desired_order),
    Experiment_Round = factor(Experiment_Round, levels = c("Main", "Follow-Up")
    )) %>%
  select(-df)

#write.csv(all_emmeans_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/emmeans_table_dye_roots_supp.csv", row.names = FALSE)
#write.csv(all_contrasts_long, "C:/Users/beabo/OneDrive/Documents/NAU/Dark Web/Comms Bio 2025/Pub_Figures/contrasts_table_dye_roots_supp.csv", row.names = FALSE)


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



ggplot(all_emmeans_long, aes(x = Type_Barrier, y = response, color = Chamber)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = roots,
             aes(x = Type_Barrier, y = preds_mod_r3*Dry_Weight_ug, color = Chamber), shape = 4,
             width = 0.15, alpha = .4, size = 7.5, stroke = 1.2,
             position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_grid(Experiment_Round~., scales = "free_y") +
  geom_text(data = n_labels,
            aes(x = Type_Barrier, y = -500, label = label),
            color = "black", 
            inherit.aes = FALSE)+
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
  )+
  theme(
    text = element_text(size = 16))


ggsave(
  "Comms Bio 2025/Pub_Figures/dye_roots_supp.png",
  width = 6.5, height = 5, units = "in", dpi = 600
)

