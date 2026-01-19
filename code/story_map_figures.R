# Aurora story map figs 
# Caylee Chan
# Created: 14 Jan 2026
# Updated: 19 Jan 2026
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
library(ggpubr)

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
  filename = "figures/story map/ICB_exp_holdingwater.jpeg",
  device = "jpeg",
  units = "in",
  height = 6, width = 8, dpi = 300
)


#### RG experiment figures -----------------------------------------------------------------------------------------

#### Clean raw ICB data ----

# Raw data 
rg_raw <- read.csv("data/Rain gardens/catch basin sampling/Aurora_RG_CBsamples_2013to2015.csv")

# Organize raw data
rg_clean <- rg_raw %>%
  mutate(Inter_type = ifelse(grepl("C", Inter_ID), "Control", "Rain Garden")) %>% # Create inter_type var
  mutate(interTypecbType = paste(Inter_type, CB_CLASS, sep = "-")) %>% # Create interTypecbType var combining the intersection type and CB type
  mutate(EpiweekYear = paste(Epiweek, Year, sep = "-")) %>% # Create new variable combining Epiweek and Year
  mutate(EpiweekYear = factor(EpiweekYear)) %>% # Change data type
  subset(CB_CLASS == "C" | CB_CLASS == "RGO" | CB_CLASS == "RGI") %>% # Subset the following CB types: C, RGO, and RGI
  subset(Area == "SE" & Year == 2015 & Epiweek %in% c(23, 25, 27, 29, 31, 33, 35)) %>% # Subset CBs in the SE area, from 2015, and in BACI Epiweeks (correspond to Epiweeks 23, 25, 27, 29, 31, 33, 35)
  subset(Inter_ID != "R16") %>% # Remove RG16 (labeled as RG intersection but no actual RG)
  mutate(interTypecbType = factor(interTypecbType)) %>% # Change data type
  mutate(Inter_ID = factor(Inter_ID)) %>% # Change data type
  mutate(CB_ID_RG = factor(CB_ID_RG)) %>% # Change data type
  mutate(Surf_Detrit = factor(Surf_Detrit, levels = c("L", "M", "H"))) # Change data type

#### Standing water presence in all catch basins  ----

## Data organizing
standingwater_all <- rg_clean %>%
  subset(Inspected == 1) %>% # Subset observations where catch basin was inspected for the presence of standing water
  drop_na(Dry) # Drop observations where Dry is NA

# Factor for proportional stacked bar chart
standingwater_all$Dry_factor <- factor(standingwater_all$Dry)

standing_water_interTypecbType <- ggplot(standingwater_all, aes(x = interTypecbType, fill = Dry_factor)) +
  geom_bar(position = "fill", color = "black") +
  labs(x = "Catch Basin Type", y = "Relative Frequency") +
  scale_fill_manual(values = c("dodgerblue1","wheat3"), labels = c("Holding\nWater", "Not Holding\nWater")) +
  scale_x_discrete(labels = c("Control\nConventional", 
                              "Rain Garden\nConventional",
                              "Rain Garden\nInfiltration",
                              "Rain Garden\nOverflow")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = seq(0,1,0.25)) +
  #annotate("text", x = 1, y = 1.02, label = expression(italic(n) == 421), size = N_size_size, family = "HelveticaNeueforSAS") +
  #annotate("text", x = 2, y = 1.02, label = expression(italic(n) == 128), size = N_size_size, family = "HelveticaNeueforSAS") +
  #annotate("text", x = 3, y = 1.02, label = expression(italic(n) == 985), size = N_size_size, family = "HelveticaNeueforSAS") +
  #annotate("text", x = 4, y = 1.02, label = expression(italic(n) == 221), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = 1.07, label = "Control Intersection", size = 5, family = "HelveticaNeueforSAS") +
  annotate("text", x = 3, y = 1.07, label = "Rain Garden Intersection", size = 5, family = "HelveticaNeueforSAS") +
  geom_bracket(data = standingwater_all, xmin = 1.55, xmax = 4.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # RG Intersection grouping
  geom_bracket(data = standingwater_all, xmin = 0.55, xmax = 1.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # Control Intersection grouping
  figtheme + theme(legend.title = element_blank())#, legend.key.size = unit(4, "mm"))


standing_water_interTypecbType

ggsave(
  plot = standing_water_interTypecbType,
  filename = "figures/story map/RG_holdingwater.jpeg",
  device = "jpeg",
  units = "in",
  height = 6, width = 8, dpi = 300
)

#### Data organizing for analyses involving juvenile mosquitoes ----

wet_CBs <- rg_clean %>% 
  subset(Dry == 0 & Sampled == 1) %>% # Select catch basins that were holding water and were sampled for juvenile mosquitoes
  mutate(interTypecbType = as.character(interTypecbType)) %>% # Change data type for subsetting
  subset(interTypecbType != "Rain Garden-RGI") %>% # Drop single obs where a RGI was holding water and sampled for mosquitoes
  mutate(interTypecbType = factor(interTypecbType)) # Convert back to factor

# Create new column for positive for juvenile mosquitoes (1 = positive for juvenile mosquitoes; 0 = negative for juvenile mosquitoes)
# All columns are 0: juvenile_pos == 1
# All columns are NA: juvenile_pos == 1
# Some columns are NA but the rest are 0: juvenile_pos == 1
# Some columns are NA but at least one non-NA column isn’t 0: juvenile_pos == 0
# At least one column isn’t equal to 0: juvenile_pos == 0
juv_mos_cols <- c("L1L3_adj", "L4_Cx_total", "P_Cx_Cx_adj", "Other_adj")

wet_CBs$juvenile_pos <- apply(wet_CBs[, juv_mos_cols], 1, function(row) {
  
  # Captures all NULL values
  if (all(is.na(row))) {
    return(1)
  }
  
  # Captures all 0 values &  some NULL but all non-NULL are zero
  if (all(row[!is.na(row)] == 0)) {
    return(1)
  }
  
  # Captures at least one non-zero value & some NULL but at least one non-NULL !0
  return(0)
})

#### Juvenile mosquito presence in all catch basin analysis ----

## Data for juvenile mosquito presence analysis in all catch basins
mosquitoes_all <- wet_CBs

#### Visualize juvenile mosquito presence ----
mosquitoes_all$juvenile_pos_factor <- factor(mosquitoes_all$juvenile_pos)


juvenile_pos_interTypecbType <- ggplot(mosquitoes_all, aes(x = interTypecbType, fill = juvenile_pos_factor)) +
  geom_bar(position = "fill", color = "black") +
  labs(x = "Catch Basin Type", y = "Relative Frequency") +
  scale_fill_manual(values = c("firebrick4","wheat3"), labels = c("Holding\nJuvenile Mosquitoes", "Not Holding\nJuvenile Mosquitoes")) +
  scale_x_discrete(labels = c("Control\nConventional", 
                              "Rain Garden\nConventinoal",
                              "Rain Garden Overflow")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), breaks = seq(0,1,0.25)) +
  #annotate("text", x = 1, y = 1.02, label = expression(italic(n) == 109), size = N_size_size, family = "HelveticaNeueforSAS") +
  #annotate("text", x = 2, y = 1.02, label = expression(italic(n) == 22), size = N_size_size, family = "HelveticaNeueforSAS") +
  #annotate("text", x = 3, y = 1.02, label = expression(italic(n) == 83), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = 1.07, label = "Control Intersection", size = 5, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.5, y = 1.07, label = "Rain Garden Intersection", size = 5, family = "HelveticaNeueforSAS") +
  geom_bracket(data = mosquitoes_all, xmin = 1.55, xmax = 3.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # RG Intersection grouping
  geom_bracket(data = mosquitoes_all, xmin = 0.55, xmax = 1.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # Control Intersection grouping
  figtheme + theme(legend.title = element_blank())#, legend.key.size = unit(4, "mm"))

juvenile_pos_interTypecbType

ggsave(
  plot = juvenile_pos_interTypecbType,
  filename = "figures/story map/RG_holdingmosquitoes.jpeg",
  device = "jpeg",
  units = "in",
  height = 6, width = 8, dpi = 300
)
