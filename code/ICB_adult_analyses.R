# Experiment #1: Adult analyses
# Caylee Chan
# Created: 9 Jan 2025
# Updated: 1 June 2026
# Notes: test for working from E drive


# Libraries:
library(dplyr)
library(ggplot2)
library(tidyr)
library(car)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(extrafont)


# Set wd
setwd("C:/Users/cayle.LAPTOP-QMLQMN4A/OneDrive - University of Illinois - Urbana/green_stormwater_infrastructure_Aurora/green_stormwater_infrastructure_Aurora")


## Read in raw data
adult_raw <- read.csv("data/Infiltration catch basins/gravid trap sampling/Aurora_ICB_GravidTrap_2013to2015_byweek.csv")

sum(is.na(adult_raw$AVG_CULEX_F_INTACT)) # No NA data

# Data cleaning
adult_raw <- adult_raw %>%
  mutate(EpiweekYear = paste(Epiweek, Year, sep = "-")) %>% # Add EpiweekYear col
  mutate(EpiweekYear = factor(EpiweekYear)) %>% # Change data type
  mutate(Area = factor(Area)) %>% # Change data type
  mutate(Period = factor(Period, levels = c("Pre-ICB", "Post-ICB"))) %>% # Change data type
  mutate(TrapID = factor(TrapID)) # Change data type

## Selecting best temperature variable
adult_tempVars <- adult_raw %>% # Create dataframe with selected temp cols and AVG_CULEX_F_INTACT
  select(NSA_TMAX_C_mean:NSA_TAVG_C_mean, NSA_avgT_d0_d1_mean:NSA_minT_d0_d1_mean, AVG_CULEX_F_INTACT)

# Correlation btwn tempVars and AVG_CULEX_F_INTACT
corrs_adult_tempVars <- cor(x = adult_tempVars[, colnames(adult_tempVars) != "AVG_CULEX_F_INTACT"], 
                            y = adult_tempVars$AVG_CULEX_F_INTACT, 
                            method = "spearman")

corrs_adult_tempVars <- as.data.frame(corrs_adult_tempVars)
colnames(corrs_adult_tempVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_adult_tempVars$absval_r <- abs(corrs_adult_tempVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_adult_tempVars$absval_r, na.rm = TRUE)

# Extract the temp variable corresponding to the greatest correlation 
corrs_adult_tempVars %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable
adult_prcpVars <- adult_raw %>% # Create dataframe with selected prcp cols and AVG_CULEX_F_INTACT
  select(NS_AVG_PRCP_mm_mean:NSA_PRCP_mm_mean, NS_PRCP_mm_sum_Lag1_mean:NS_PRCP_mm_sum_d14_d20_mean, AVG_CULEX_F_INTACT)

# Correlation btwn prcpVars and AVG_CULEX_F_INTACT
corrs_adult_prcpVars <- cor(x = adult_prcpVars[, colnames(adult_prcpVars) != "AVG_CULEX_F_INTACT"], 
                            y = adult_prcpVars$AVG_CULEX_F_INTACT,
                            method = "spearman")

corrs_adult_prcpVars <- as.data.frame(corrs_adult_prcpVars)
colnames(corrs_adult_prcpVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_adult_prcpVars$absval_r <- abs(corrs_adult_prcpVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_adult_prcpVars$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_adult_prcpVars %>%
  filter(absval_r == maxAbsval_r)


# Generalized linear mixed effects model
# Response: AVG_CULEX_F_INTACT

# Predictor variables:
# NS_PRCP_mm_sum_d1_d7_mean
# NSA_TMAX_C_mean
# Area (control, ICB, ICB buffer)
# Period
# Area * Period

# Random effects
# TrapID
# EpiweekYear

# Check data types
str(adult_raw$AVG_CULEX_F_INTACT)
str(adult_raw$NS_PRCP_mm_sum_d1_d7_mean)
str(adult_raw$NSA_TMAX_C_mean)
str(adult_raw$Area)
str(adult_raw$Period)
str(adult_raw$TrapID)
str(adult_raw$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models

# Full model
full_model_adult1 <- glmmTMB(AVG_CULEX_F_INTACT ~ NS_PRCP_mm_sum_d1_d7_mean + NSA_TMAX_C_mean + Area + Period + Area*Period + (1|TrapID) + (1|EpiweekYear),
                             data = adult_raw,
                             family = Gamma(link = "log"))

# Reduced model 1: removed temp var
red_model1_adult1 <- glmmTMB(AVG_CULEX_F_INTACT ~ NS_PRCP_mm_sum_d1_d7_mean + Area + Period + Area*Period + (1|TrapID) + (1|EpiweekYear),
                             data = adult_raw,
                             family = Gamma(link = "log"))

# Reduced model 2: removed prcp var
red_model2_adult1 <- glmmTMB(AVG_CULEX_F_INTACT ~ NSA_TMAX_C_mean + Area + Period + Area*Period + (1|TrapID) + (1|EpiweekYear),
                             data = adult_raw,
                             family = Gamma(link = "log"))

# Reduced model 3: removed temp and prcp vars
red_model3_adult1 <- glmmTMB(AVG_CULEX_F_INTACT ~ Area + Period + Area*Period + (1|TrapID) + (1|EpiweekYear),
                             data = adult_raw,
                             family = Gamma(link = "log"))

# Compare models
model.sel(full_model_adult1, red_model1_adult1, red_model2_adult1, red_model3_adult1,
          rank = AICc)

# Selected model: reduced model #2

# Summary
adult_summary <- summary(red_model2_adult1)
adult_summary

confint(red_model2_adult1, level = 0.95)

# Anova
Anova(red_model2_adult1, type = 3)

# Emmeans
emmeans(red_model2_adult1, ~ Area*Period, type = "response", adjust = "sidak")

# By Period by Area  class
emobject_adult1 <- emmeans(red_model2_adult1, ~ Area*Period, type = "response", adjust = "sidak")
contrast(emobject_adult1, "revpairwise", by = "Period", type = "response", adjust = "sidak")

# By Area by Period
emobject2_adult1 <- emmeans(red_model2_adult1, ~ Area*Period, type = "response", adjust = "sidak")
contrast(emobject2_adult1, "revpairwise", by = "Area", type = "response", adjust = "sidak")


# Model diagnostics

# Calculate residuals
simulationOutput_adult1 <- simulateResiduals(fittedModel = red_model2_adult1, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_adult1)
plotQQunif(simulationOutput_adult1) # Left plot
plotResiduals(simulationOutput_adult1) # Right plot


#### Proportional change figure ----

# Calculate mean AVG_CULEX_F_INTACT for each TrapID by Area and Period
adultData_intact <- adult_raw %>%
  group_by(Period, Area, TrapID) %>%
  summarise(
    mean_Cx_F_intact = mean(AVG_CULEX_F_INTACT))

# Separate the mean_Cx_F_intact into pre and post dataframes & update column name
adultData_intact_PRE <- subset(adultData_intact, Period == "Pre-ICB")
colnames(adultData_intact_PRE)[4] <- "mean_Cx_F_intact_PRE"

adultData_intact_POST <- subset(adultData_intact, Period == "Post-ICB")
colnames(adultData_intact_POST)[4] <- "mean_Cx_F_intact_POST"

# Combine pre and post mean_Cx_F_intact into one dataframe 
adultData2_intact <- merge(adultData_intact_PRE, 
                           adultData_intact_POST[, c("TrapID", "mean_Cx_F_intact_POST")], 
                           by = "TrapID")

# Proportional change in mean_Cx_F_intact by TrapID from pre to post
adultData2_intact$prop_change_mean_Cx_F_intact <- (adultData2_intact$mean_Cx_F_intact_POST - adultData2_intact$mean_Cx_F_intact_PRE) / adultData2_intact$mean_Cx_F_intact_PRE

## Calculating sample sizes and summary stats
adult_raw %>%
  group_by(Area, Period) %>%
  summarise(
    total_trapping_events = n(),
    avg_adult = mean(AVG_CULEX_F_INTACT),
    se_adult = sd(AVG_CULEX_F_INTACT) / sqrt(total_trapping_events)
  )


## Boxplot

# Change xaxis order
adultData2_intact$Area <- factor(adultData2_intact$Area, levels = c("control", "ICB buffer", "ICB"))

# Set theme
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
  legend.position = "top")

adult_prop_change_intact <- ggplot(adultData2_intact, aes(x = Area, y = prop_change_mean_Cx_F_intact, fill = Area, color = Area)) + 
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.75), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19) +  
  scale_x_discrete(labels = c("Control Area", "Buffer Area", "Treatment Area")) +
  scale_fill_manual(values = c("grey27", "mediumblue", "green4")) +
  scale_color_manual(values = c("black", "black", "black")) +
  xlab("Area") +
  ylab("Proportional Change in Average Female Culex spp.\nAbundance from Pre- to Post-Construction Period") +
  figtheme +
  theme(
    legend.position = "none"
  )

adult_prop_change_intact


ggsave(
  plot = adult_prop_change_intact,
  filename = "figures/manuscript/Figure4.jpeg",
  device = "jpeg", 
  units = "cm", 
  height = 8.5, width = 8.5, dpi = 300
)

ggsave(
  plot = adult_prop_change_intact,
  filename = "figures/manuscript/Figure4.pdf",
  device = "pdf", 
  units = "cm", 
  height = 8.5, width = 8.5, dpi = 300
)

#### Rainfall during pre- and post-construction periods ####
library(lubridate)

# Precipitation data
all_weather_13_15 <- read.csv("data/weather/Aurora_weather_2013_2015_all_months.csv")

str(all_weather_13_15$DATE)

# Data cleaning
weather_clean <- all_weather_13_15 %>%
  mutate(Date_new = as.Date(DATE, format = "%Y-%m-%d")) %>%
  mutate(year = year(Date_new)) %>%
  mutate(month = month(Date_new)) %>%
  mutate(day = day(Date_new)) %>%
  subset(month %in% c(4, 5, 6, 7, 8, 9)) %>%
  mutate(month = factor(month, levels = c(4, 5, 6, 7, 8, 9), labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep"))) %>%
  mutate(year = factor(year))

# Sum precip by month
Aurora_prcp_summary <- weather_clean %>%
  group_by(year, month) %>%
  summarise(
    total_prcp = sum(PRCP)
  )

# Read in Aurora prcip norms data (monthly)
Aurora_precip_norms <- read.csv("data/weather/Aurora_monthly_prcp_normals_1981_2010.csv")
str(Aurora_precip_norms)
precip_norms_clean <- Aurora_precip_norms %>%
  mutate(prcp_mm = MLY.PRCP.NORMAL * 25.4) %>%
  subset(DATE %in% c(4, 5, 6, 7, 8, 9)) %>%
  mutate(DATE = factor(DATE, levels = c(4, 5, 6, 7, 8, 9), labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep")))

# Merge dfs together
prcp_13_15_and_normals <- merge(x = Aurora_prcp_summary,
                                y = precip_norms_clean,
                                by.x = "month", 
                                by.y = "DATE")

str(prcp_13_15_and_normals)
prcp_13_15_and_normals$group2 <- "Rainfall normals\n(1981-2010)"

precip_13_15 <- ggplot(prcp_13_15_and_normals) +
  geom_bar(aes(x = month, y = total_prcp, fill = year), stat = "identity", position = position_dodge(0.75), width = 0.75, color = "black", linewidth = 0.25) +
  geom_line(aes(x = month, y = prcp_mm, group = group2, linetype = group2)) +
  geom_point(aes(x = month, y = prcp_mm, group = group2, shape = group2), fill = "white", color = "black", size = 2) +
  scale_shape_manual(values = c(21)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,275), breaks = seq(0, 300, 50)) +
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086")) +
  labs(x = "Month", y = "Precipitation (mm)") +
  figtheme +
  theme(
    legend.title = element_blank(),
    legend.key.width = unit(4, "mm"),
    legend.key.height = unit(4, "mm")
  ) 

precip_13_15

# Save to Aurora Box folder
ggsave(
  plot = precip_13_15,
  filename = "figures/manuscript/precipitation_2013_15.jpeg",
  device = "jpeg", 
  units = "cm", 
  height = 8.5, width = 8.5, dpi = 300
)

ggsave(
  plot = precip_13_15,
  filename = "figures/manuscript/precipitation_2013_15.pdf",
  device = "pdf", 
  units = "cm", 
  height = 8.5, width = 8.5, dpi = 300
)
