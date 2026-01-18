# Aurora story map figs 
# Caylee Chan
# Created: 14 Jan 2026
# Updated: 18 Jan 2026
# Notes: 

# Libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(patchwork)
library(extrafont)
library(ggtext)

# Set wd
setwd("C:/Users/cayle.LAPTOP-QMLQMN4A/OneDrive - University of Illinois - Urbana/green_stormwater_infrastructure_Aurora/green_stormwater_infrastructure_Aurora")

#### Specify plotting elements ----

# Set theme
figtheme <- theme(
  strip.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = 'black', linewidth = 0.3), 
  axis.title.x = element_text(size = 18, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.title.y = element_text(size = 18, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.text.x = element_text(size = 16, color="black", family = "HelveticaNeueforSAS"), 
  axis.text.y = element_text(size = 12, color="black", family = "HelveticaNeueforSAS"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  legend.text = element_text(size = 14, family = "HelveticaNeueforSAS"), 
  legend.position = "top",
  legend.title = element_text(size = 14, color="black", family = "HelveticaNeueforSAS"))

#### ICB experiment figures -----------------------------------------------------------------------------------------

#### Clean raw ICB data ----

# Raw data
icb_raw <- read.csv("data/Infiltration catch basins/catch basin sampling/Aurora_ICB_CBsamples_2013to2015.csv")

# Organize raw data
icb_clean <- icb_raw %>% # Select CBs of interest
  subset((Replace_class == "C-C" & CB_type == "C" & Area == "Control" & Year %in% c(2013, 2014, 2015)) | # Control C-C CBs from all years
           (Replace_class == "C-C" & CB_type == "C" & Area == "Treatment" & Year %in% c(2013, 2014, 2015)) | # Treatment C-C CBs from all years
           (Replace_class == "C-I" & CB_type == "C" & Area == "Treatment" & Year == 2013) | # Treatment C-I CBs that were C in 2013
           (Replace_class == "C-I" & CB_type == "I" & Area == "Treatment" & Year %in% c(2014, 2015))) %>% # Treatment C-I CBs that were I in 2014 & 2015
  mutate(AreaCB_Rclass = paste(Area, Replace_class, sep = "-")) %>% # Create new variable combining Area and Replace_class
  mutate(AreaCB_Rclass = factor(AreaCB_Rclass)) %>% # Make CB type a factor
  mutate(EpiweekYear = paste(Epiweek, Year, sep = "-")) %>% # Create new variable combining Epiweek and Year
  mutate(EpiweekYear = ifelse(EpiweekYear == "24.5-2015", "25-2015", EpiweekYear)) %>% # Epiweek 24.5 in 2015 typo - should be Epiweek 25
  mutate(EpiweekYear = factor(EpiweekYear)) %>% # Make epiweek-year a factor
  mutate(Time = ifelse(Year == 2013, "Pre ICB Construction", "Post ICB Construction")) %>% # Create pre and post construction time blocks
  mutate(Time = factor(Time, levels = c("Pre ICB Construction", "Post ICB Construction"))) %>% # Make time a factor
  mutate(CB_ID = factor(CB_ID)) %>% # Make CB ID a factor
  mutate(Surf_Detrit = factor(Surf_Detrit, levels = c("L", "M", "H"))) %>% # Make surface detritus level a factor
  subset(ICB_wk_type == "BACI") # Select BACI (two-week interval) obs

#### Standing water presence ----

# Data for standing water analyses
standingwater <- icb_clean %>%
  subset(Inspected == 1) %>% # Select catch basins that were inspected for standing water presence
  drop_na(Dry) # Drop obs where standing water presence was NA

## Mixed effects logistic regression model
# Response variable: Dry (0 = water present; 1 = no water present)

# Predictor variables:
# NSA_TMIN_C
# NS_PRCP_mm_sum_Lag4 
# AreaCB_Rclass
# Time
# AreaCB_Rclass * Time

# Random effects
# CB_ID
# EpiweekYear

# Check data types
str(standingwater$Dry)
str(standingwater$NSA_TMIN_C)
str(standingwater$NS_PRCP_mm_sum_Lag4)
str(standingwater$AreaCB_Rclass)
str(standingwater$Time)
str(standingwater$CB_ID)
str(standingwater$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

# SELECTED MODEL: Reduced model 1: dropped temp var
red_model1_water_pres <- glmmTMB(Dry ~ NS_PRCP_mm_sum_Lag4 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                 data = standingwater,
                                 family = binomial(link= "logit"))

# Emmeans
em_water_pres <- emmeans(red_model1_water_pres, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
em_water_pres
contrast(em_water_pres, "revpairwise", by = "Time", type = "response", adjust = "sidak")

## Set plotting elements

# Legend
legend_colors <- c("grey27", "grey60", "green4")
legend_labs <- c("**C-C** Catch Basin<br>(Control Area)", "**C-C** Catch Basin<br>(Treatment Area)", "**C-I** Catch Basin<br>(Treatment Area)")

# Set emmip arguments
linearg <- list(size = 1.5, linetype = "solid")
dotarg <- list(size = 3, shape = "circle")
CIarg <- list(lwd = 0.75, alpha = 0.3)

# Plot
ICB_holding_water_plot <- emmip(object = em_water_pres, 
                                formula = AreaCB_Rclass ~ Time, 
                                CIs = TRUE,
                                engine = "ggplot",
                                linearg = linearg,
                                dotarg = dotarg,
                                CIarg = CIarg) +
  scale_color_manual(values = legend_colors, labels = legend_labs) +
  #scale_linetype_manual(values = c("solid", "solid", "solid"), labels = legend_labs, guide = "none", size = 1) +
  #scale_size_manual(values = c(2,2,2), labels = legend_labs, guide = "none") +
  #aes(linetype = AreaCB_Rclass, size = AreaCB_Rclass, color = AreaCB_Rclass) +
  labs(color = "") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.015,1.1)) +
  scale_x_discrete(labels = c("Pre","Post")) +
  ylab("Estimated Probability of a Catch <br>Basin <span style='color:#FF0000;'>Not</span> Holding Water") +
  xlab("Construction Period") +
  #xlab("") +
  figtheme + theme(legend.text = element_markdown(),
                   axis.title.y = element_markdown())

ICB_holding_water_plot

ggsave(
  plot = ICB_holding_water_plot,
  filename = "figures/ICB_exp_holdingwater.jpeg",
  device = "jpeg",
  units = "in",
  height = 6, width = 8, dpi = 300
)

