# Experiment #2: Starvation resistance analysis
# Caylee Chan
# Created: 13 Jan 2025
# Updated: 23 Jan 2026
# Notes: 


# Libraries:
library(dplyr)
library(ggplot2)
library(tidyr)
library(car)
library(emmeans)
library(MuMIn)
library(coxme)
library(DHARMa)
library(survival)
library(survminer)
library(patchwork)
library(extrafont)

# Set wd
setwd("C:/Users/cayle.LAPTOP-QMLQMN4A/OneDrive - University of Illinois - Urbana/green_stormwater_infrastructure_Aurora/green_stormwater_infrastructure_Aurora")

# Set theme
survtheme <- theme(
  strip.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = 'black', linewidth = 0.3), 
  axis.title.x = element_text(size = 7, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.title.y = element_text(size = 7, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.text.x = element_text(size = 7, color="black", family = "HelveticaNeueforSAS"), 
  axis.text.y = element_text(size = 7, color="black", family = "HelveticaNeueforSAS"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 7), 
  legend.position = "top",
  plot.title = element_text(size = 7, face = "bold", family = "HelveticaNeueforSAS")
)

figtheme <- theme(
  strip.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = 'black', linewidth = 0.3), 
  axis.title.x = element_text(size = 7, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.title.y = element_text(size = 7, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.text.x = element_text(size = 7, color="black", family = "HelveticaNeueforSAS"), 
  axis.text.y = element_text(size = 7, color="black", family = "HelveticaNeueforSAS"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  legend.text = element_text(size = 7, family = "HelveticaNeueforSAS"), 
  legend.position = "top",
  legend.title = element_blank())

# Set N size size
N_size_size <- 2

#### Data organizing and sample sizes ----

# Raw data
starv_raw <- read.csv("data/Rain gardens/pupal sampling/starvation_resistance/Aurora_2015_starvation_resistance_20241007.csv")

## Data organizing
starv_raw <- starv_raw %>%
  mutate(CB_Class = factor(CB_Class)) %>%
  mutate(CBID = factor(CBID)) %>%
  mutate(Epiweek = factor(Epiweek))

# Females
females <- starv_raw %>%
  subset(Species == "Cx. pipiens" & Sex == "F") %>% # Select only Cx. pipiens and females
  drop_na(Survival_days, Winglength_mm) # Drop obs with missing survival days and wing length

females$status <- ifelse(!is.na(females$Survival_days), 1, 0) # status = 1 died by completion of study; status = 0 didn't die by completion of the study (censored)

# Female sample sizes and summary stats
females %>%
  group_by(CB_Class) %>%
  summarise(
    total_obs = n(),
    avg_survivaldays = mean(Survival_days),
    se_survivaldays = sd(Survival_days) / sqrt(total_obs)
  )

# Males
males <- starv_raw %>%
  subset(Species == "Cx. pipiens" & Sex == "M") %>% # Select only Cx. pipiens and males
  drop_na(Survival_days, Winglength_mm) # Drop obs with missing survival days

males$status <- ifelse(!is.na(males$Survival_days), 1, 0) # status = 1 died by completion of study; status = 0 didn't die by completion of the study (censored)

# Male sample sizes and summary stats
males %>%
  group_by(CB_Class) %>%
  summarise(
    total_obs = n(),
    avg_survivaldays = mean(Survival_days),
    se_survivaldays = sd(Survival_days) / sqrt(total_obs)
  )

#### Cox proportional hazards model (FEMALES) ----

# Response variable: Survival_days

# Predictor variables:
# CB_Class 
# Winglength_mm

# Random effects:
# CBID
# Epiweek

# Check data types
str(females$Survival_days)
str(females$CB_Class)
str(females$Winglength_mm)
str(females$CBID)
str(females$Epiweek)
str(females$status)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models - females

# Full model

full_model_F <- coxme(Surv(Survival_days, status) ~ CB_Class + Winglength_mm + (1|CBID) + (1|Epiweek),
                      data = females)

# Reduced model #1: dropped winglength variable
red_model1_F <- coxme(Surv(Survival_days, status) ~ CB_Class + (1|CBID) + (1|Epiweek),
                      data = females)

# Reduced model #2: dropped catch basin class variable
red_model2_F <- coxme(Surv(Survival_days, status) ~ Winglength_mm + (1|CBID) + (1|Epiweek),
                      data = females)

# Compare models - females
model.sel(full_model_F, red_model1_F, red_model2_F,
          rank = AICc)

# Summary
summary(full_model_F)

exp(confint(full_model_F, level = 0.95))

# Anova
Anova(full_model_F, type = 3)

# Emmeans
emobject_females <- emmeans(full_model_F, specs = pairwise ~ CB_Class, type = "response", adjust = "tukey")
contrast(emobject_females, "revpairwise", type = "response", adjust = "tukey")


#### Female plot ----
femaleSurvCurve <- survfit(Surv(Survival_days, status) ~ CB_Class, data = females)

femaleSurvCurvePlot <- ggsurvplot(femaleSurvCurve,
                                  break.x.by = 0.5,
                                  axes.offset = FALSE,
                                  conf.int = TRUE,
                                  xlim = c(0,7),
                                  risk.table = FALSE,
                                  size = 0.5,
                                  legend.labs = c("Conventional\nCatch Basin", "Rain Garden Overflow\nCatch Basin"),
                                  palette = c("grey27", "green4"),
                                  ggtheme = survtheme) 

femaleSurvCurvePlot <- femaleSurvCurvePlot$plot + 
  xlab("Time (Days)") +
  ylab("Survival Probability") +
  scale_y_continuous(limits = c(0,1.01), breaks = seq(0,1,0.25)) +
  annotate("text", x = 1, y = 0.41, label = "Females", size = 3, family = "HelveticaNeueforSAS") +
  annotate("rect", xmin = 0.61, xmax = 1.39, ymin = 0.28, ymax = 0.32, fill = "grey27", alpha = 0.3) + # CCB
  annotate("rect", xmin = 0.61, xmax = 1.39, ymin = 0.23, ymax = 0.27, fill = "green4", alpha = 0.3) + # RGO
  annotate("text", x = 1, y = 0.3, label = expression(italic(n) == 188), size = N_size_size, family = "HelveticaNeueforSAS") + # CCB
  annotate("text", x = 1, y = 0.25, label = expression(italic(n) == 157), size = N_size_size, family = "HelveticaNeueforSAS") + # RGO
  theme(
    legend.position = "bottom"
  )


femaleSurvCurvePlot

#### Cox proportional hazards model (MALES) ----

# Response variable: Survival_days

# Predictor variables:
# CB_Class
# Winglength_mm

# Random effects:
# CBID
# Epiweek

# Check data types
str(males$Survival_days)
str(males$CB_Class)
str(males$Winglength_mm)
str(males$CBID)
str(males$Epiweek)
str(males$status)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models - males

# Full model

full_model_M <- coxme(Surv(Survival_days, status) ~ CB_Class + Winglength_mm + (1|CBID) + (1|Epiweek),
                      data = males)

# Reduced model #1: dropped winglength variable
red_model1_M <- coxme(Surv(Survival_days, status) ~ CB_Class + (1|CBID) + (1|Epiweek),
                      data = males)

# Reduced model #2: dropped catch basin class variable
red_model2_M <- coxme(Surv(Survival_days, status) ~ Winglength_mm + (1|CBID) + (1|Epiweek),
                      data = males)

# Compare models - males
model.sel(full_model_M, red_model1_M, red_model2_M,
          rank = AICc)

# Selected model: full model

# Summary
summary(full_model_M)

exp(confint(full_model_M, level = 0.95))

# Anova
Anova(full_model_M, type = 3)

# Emmeans
emobject_males <- emmeans(full_model_M, specs = pairwise ~ CB_Class, type = "response", adjust = "tukey")
contrast(emobject_males, "revpairwise", type = "response", adjust = "tukey")


#### Male plot ----

maleSurvCurve <- survfit(Surv(Survival_days, status) ~ CB_Class, data = males)

maleSurvCurvePlot <- ggsurvplot(maleSurvCurve,
                                break.x.by = 0.5,
                                axes.offset = FALSE,
                                conf.int = TRUE,
                                xlim = c(0,7),
                                risk.table = FALSE,
                                size = 0.5,
                                legend.labs = c("Conventional\nCatch Basin", "Rain Garden Overflow\nCatch Basin"),
                                palette = c("grey27", "green4"),
                                ggtheme = survtheme) 

maleSurvCurvePlot <- maleSurvCurvePlot$plot + 
  xlab("Time (Days)") +
  ylab("Survival Probability") +
  scale_y_continuous(limits = c(0,1.01), breaks = seq(0,1,0.25)) +
  annotate("text", x = 1, y = 0.41, label = "Males", size = 3, family = "HelveticaNeueforSAS") + 
  annotate("rect", xmin = 0.61, xmax = 1.39, ymin = 0.28, ymax = 0.32, fill = "grey27", alpha = 0.3) + # CCB
  annotate("rect", xmin = 0.61, xmax = 1.39, ymin = 0.23, ymax = 0.27, fill = "green4", alpha = 0.3) + # RGO
  annotate("text", x = 1, y = 0.3, label = expression(italic(n) == 185), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = 0.25, label = expression(italic(n) == 142), size = N_size_size, family = "HelveticaNeueforSAS") +
  theme(
    legend.position = "bottom"
  )

maleSurvCurvePlot


#### Combination figures ----

## Winglength
# Raw winglength data
raw_wing <- read.csv("data/Rain gardens/pupal sampling/starvation_resistance/Aurora_2015_starvation_resistance_20241007.csv")

# Subset Cx. pipiens and remove NA wing lengths
wing_clean_all <- raw_wing %>%
  filter(!is.na(Winglength_mm)) %>%
  subset(Species == "Cx. pipiens") 

## Wing length sample sizes

# Females - CCBs
females_CCBs_wing <- wing_clean_all %>%
  subset(Sex == "F" & CB_Class == "C")
nrow(females_CCBs_wing) # Number of obs females CCBs

# Females - RGOs
females_RGOs_wing <- wing_clean_all %>%
  subset(Sex == "F" & CB_Class == "RGO")
nrow(females_RGOs_wing) # Number of obs females RGOs

# Males - CCBs
males_CCBs_wing <- wing_clean_all %>%
  subset(Sex == "M" & CB_Class == "C")
nrow(males_CCBs_wing) # Number of obs males CCBs

# Males - RGOs
males_RGOs_wing <- wing_clean_all %>%
  subset(Sex == "M" & CB_Class == "RGO")
nrow(males_RGOs_wing) # Number of obs males RGOs

# Number of CBs by type and sampling periods
n_distinct(wing_clean_all$CBID[wing_clean_all$CB_Class == "C"])
n_distinct(wing_clean_all$CBID[wing_clean_all$CB_Class == "RGO"])
n_distinct(wing_clean_all$DateColl)

# Wing length figure
winglength <- ggplot(wing_clean_all, aes(x = Sex, y = Winglength_mm, fill = CB_Class, color = CB_Class)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  scale_fill_manual(values = c("grey27", "green4"), labels = c("Conventional\nCatch Basin", "Rain Garden Overflow\nCatch Basin")) +
  scale_color_manual(values = c("black", "black", "black"), labels = c("Conventional\nCatch Basin", "Rain Garden Overflow\nCatch Basin")) +
  scale_x_discrete(labels = c("Females", "Males")) +
  scale_y_continuous(expand = c(0,0), limits = c(1.6,4.75), breaks = seq(1.5, 4.5, 0.5)) +
  labs(x = "Sex", y = "Wing Length (mm)") +
  annotate("text", x = 0.79, y = 1.8, label = expression(italic(n) == 189), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.23, y = 1.8, label = expression(italic(n) == 157), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.79, y = 1.8, label = expression(italic(n) == 185), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.23, y = 1.8, label = expression(italic(n) == 142), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme +
  theme(
    legend.position = "bottom"
  )

winglength

## Diapause
# Raw diapause data
diapause <- read.csv("data/Rain gardens/pupal sampling/diapause/Aurora_2015_diapause_data_20241203_w_watertemp.csv")

# Diapause data organizing
diapause <- diapause %>%
  mutate(CB_CLASS = factor(CB_CLASS)) %>%
  mutate(CB_ID = factor(CB_ID)) 

diapause %>%
  group_by(CB_CLASS) %>%
  summarise(
    total_obs = n(),
    total_indiapause = sum(Diapause == 1),
    percent_indiapause = (total_indiapause / total_obs) * 100
  )

diapause$Diapause_factor <- factor(diapause$Diapause, levels = c(1, 0), labels = c("Diapausing Female", "Non-diapausing Female"))
str(diapause$Diapause_factor)

diapause_viz <- ggplot(diapause, aes(x = CB_CLASS, fill = Diapause_factor)) +
  geom_bar(position = "fill", color = "black") +
  labs(x = "Catch Basin Type", y = "Relative Frequency") +
  scale_fill_manual(values = c("green4", "grey57")) + 
  scale_x_discrete(labels = c("Conventional\nCatch Basin", "Rain Garden Overflow\nCatch Basin")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05)) +
  annotate("text", x = 1, y = 1.02, label = expression(italic(n) == 77), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = 1.02, label = expression(italic(n) == 67), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme + theme(legend.key.size = unit(4, "mm")) +
  theme(
    legend.position = "top"
  )

diapause_viz

Fig6_top <- ((femaleSurvCurvePlot | maleSurvCurvePlot)) +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", plot.tag.position = c(0.02, 1.01), plot.tag = element_text(size = 9, face = "bold", family = "HelveticaNeueforSAS"))

Fig6 <- Fig6_top / (winglength | diapause_viz) +
  #plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", plot.tag.position = c(0.02, 1.01), plot.tag = element_text(size = 9, face = "bold", family = "HelveticaNeueforSAS"))

Fig6

ggsave(
  plot = Fig6,
  filename = "figures/manuscript/Figure6.jpeg",
  device ="jpeg",
  units = "mm",
  height = 185, width = 180, dpi = 300
)


