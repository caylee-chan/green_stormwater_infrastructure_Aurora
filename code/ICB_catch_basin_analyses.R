# Experiment #1: ICB analyses (all)
# Caylee Chan
# Created: 2 Jan 2025
# Updated: 22 Jan 2025
# Notes: 

# Libraries:
library(dplyr)
library(ggplot2)
library(tidyr)
library(car)
library(formula.tools)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(performance)
library(patchwork)
library(extrafont)


# Set wd
setwd("C:/Users/cayle.LAPTOP-QMLQMN4A/OneDrive - University of Illinois - Urbana/green_stormwater_infrastructure_Aurora/green_stormwater_infrastructure_Aurora")

#### Specify plotting elements ----

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
  legend.position = "top",
  legend.title = element_text(size = 7, color="black", family = "HelveticaNeueforSAS"))

# Set legend labels
legend_labs <- c("Control C-C", "Treatment C-C", "Treatment C-I")
legend_labs <- c("Control CCB", "Treatment CCB", "Treatment ICB")
legend_colors <- c("grey27", "grey60", "green4")

# Set N size size
N_size_size <- 2

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

#### Standing water presence in catch basins analysis ----

# Data for standing water analyses
standingwater <- icb_clean %>%
  subset(Inspected == 1) %>% # Select catch basins that were inspected for standing water presence
  drop_na(Dry) # Drop obs where standing water presence was NA

## Selecting best temperature variable
tempVars <- icb_raw %>% # Create dataframe with selected temp cols
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1)
tempVars <- colnames(tempVars) # Create a vector of the temp by extracting the col names

# Declare vectors
Predictor <- vector() # Holds the variable being used for prediction
Predictor_Coef <- vector() # Holds the coefficient of the predictor variable
Model_AIC <- vector() # Holds the AIC score for the model
Model_Formula <- vector() # Holds the model formula

# Specify data
for (i in tempVars) { # Specify prcpVars or tempVars
  formula <- as.formula(paste("Dry ~", i)) # Change response variable as needed!!!!
  print(formula)
  model <- glm(formula, data = standingwater, family = binomial(link = "logit"))
  Predictor[i] = i
  Predictor_Coef[i] = model$coefficients[2]
  Model_AIC[i] = model$aic
  Model_Formula[i] = as.character(model$formula)
}

# Create dataframe to compare models
compareModels_df <- data.frame(cbind(Predictor, Predictor_Coef, Model_AIC, Model_Formula))
nrow(compareModels_df) # Triple check all predictors accounted for

# Find the model with the smallest AIC score
smallest_AIC <- min(compareModels_df$Model_AIC) # Determine the smallest AIC score
compareModels_df[Model_AIC == smallest_AIC, ] # Select the row with the smallest AIC score


## Selecting best precipitation variable
prcpVars <- icb_raw %>%
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20)
prcpVars <- colnames(prcpVars)

# Declare vectors
Predictor <- vector() # Holds the variable being used for prediction
Predictor_Coef <- vector() # Holds the coefficient of the predictor variable
Model_AIC <- vector() # Holds the AIC score for the model
Model_Formula <- vector() # Holds the model formula

# Specify data
for (i in prcpVars) { # Specify prcpVars or tempVars
  formula <- as.formula(paste("Dry ~", i)) # Change response variable as needed!!!!
  print(formula)
  model <- glm(formula, data = standingwater, family = binomial(link = "logit"))
  Predictor[i] = i
  Predictor_Coef[i] = model$coefficients[2]
  Model_AIC[i] = model$aic
  Model_Formula[i] = as.character(model$formula)
}

# Create dataframe to compare models
compareModels_df <- data.frame(cbind(Predictor, Predictor_Coef, Model_AIC, Model_Formula))
nrow(compareModels_df) # Triple check all predictors accounted for

# Find the model with the smallest AIC score
smallest_AIC <- min(compareModels_df$Model_AIC) # Determine the smallest AIC score
compareModels_df[Model_AIC == smallest_AIC, ] # Select the row with the smallest AIC score


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
options(contrasts = c("contr.sum","contr.poly")) # Comparisons to the grand mean, NOT to only one baseline (reference) level

## Models

# Full model
full_model_water_pres <- glmmTMB(Dry ~ NSA_TMIN_C + NS_PRCP_mm_sum_Lag4 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                 data = standingwater,
                                 family = binomial(link= "logit"))


# Reduced model 1: dropped temp var
red_model1_water_pres <- glmmTMB(Dry ~ NS_PRCP_mm_sum_Lag4 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                 data = standingwater,
                                 family = binomial(link= "logit"))


# Reduced model 2: dropped prcp var
red_model2_water_pres <- glmmTMB(Dry ~ NSA_TMIN_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                 data = standingwater,
                                 family = binomial(link= "logit"))


# Reduced model 3: dropped temp and prcp vars
red_model3_water_pres <- glmmTMB(Dry ~ AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                 data = standingwater,
                                 family = binomial(link= "logit"))


# Model selection
model.sel(full_model_water_pres, red_model1_water_pres, red_model2_water_pres, red_model3_water_pres,
          rank = AICc)

# Selected model: reduced model #1

# Summary
summary_holdingwater <- summary(red_model1_water_pres)
summary_holdingwater
exp(summary_holdingwater$coefficients$cond)
exp(confint(red_model1_water_pres, level = 0.95))

# Anova
Anova(red_model1_water_pres, type = 3)

# Emmeans
emmeans(red_model1_water_pres, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")

# By Year by Time
emobject_water_pres <- emmeans(red_model1_water_pres, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_water_pres, "revpairwise", by = "Time", type = "response", adjust = "sidak")

# By AreaCB_Rclass by Time
emobject_water_pres2 <- emmeans(red_model1_water_pres, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_water_pres2, "revpairwise", by = "AreaCB_Rclass", type = "response", adjust = "sidak")


## Model diagnostics

# Calculate residuals
simulationOutput_waterPres <- simulateResiduals(fittedModel = red_model1_water_pres, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_waterPres)
plotQQunif(simulationOutput_waterPres) # Left plot
plotResiduals(simulationOutput_waterPres) # Right plot

#### Visualize standing water presence ----

## Holding water sample sizes

# Sample sizes and summart stats
standingwater %>%
  group_by(AreaCB_Rclass, Time) %>%
  summarise(
    total_obs = n(),
    obs_NOT_holdingwater = sum(Dry == 1),
    percent_NOT_holdingwater = round(obs_NOT_holdingwater/total_obs *100, 2)
  )

# Calculate the proportion of catch basins NOT holding water by construction period, epiweek, and type
standingwater_figuredata <- standingwater %>%
  group_by(Time, EpiweekYear, AreaCB_Rclass) %>%
  summarise(
    prop_NOT_holding_water = sum(Dry == 1, na.rm = TRUE)/sum(Dry == 1 | Dry == 0, na.rm = TRUE),
    groupSize = n(),
    .groups = "drop"
  )

# Plot
holdingwater_fig <- ggplot(standingwater_figuredata, aes(x = Time, y = prop_NOT_holding_water, fill = AreaCB_Rclass, color = AreaCB_Rclass)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Construction Period", y = "Proportion of Catch Basins Not Holding Water") +
  scale_color_manual(values = c("black", "black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  scale_x_discrete(labels = c("Pre","Post")) +
  annotate("text", x = 0.7, y = -0.04, label = expression(italic(n) == 257), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = -0.04, label = expression(italic(n) == 102), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.3, y = -0.04, label = expression(italic(n) == 158), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.7, y = -0.04, label = expression(italic(n) == 528), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -0.04, label = expression(italic(n) == 236), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.3, y = -0.04, label = expression(italic(n) == 392), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme 

holdingwater_fig                          


#### Data organizing for analyses involving juvenile mosquitoes ----

wet_CBs <- icb_clean %>% 
  subset(Dry == 0 & Sampled == 1) # Select catch basins that were holding water and were sampled for juvenile mosquitoes

# Verify that samples that weren't processed were field neg and should have 0s in the juvenile mosquito-related columns

# Subset wet CBs where Processed == 0 (samples WEREN'T processed)
processed0 <- subset(wet_CBs, Processed == 0)

# Subset CBs in processed0 where all juvenile mosquito-related cols are NULL
juv_mos_cols <- c("L1L3_adj", "L4_Cx_total", "P_Cx_Cx_adj", "Other_adj")

# Checks each row to see if all values are NA in the juv_mos_cols
processed0_all_null <- processed0[apply(processed0[juv_mos_cols], 1, function(x) all(is.na(x))), ]

# Verify that the lengths of these two dataframes are the same
# This means that all wet CBs with NA for all juvenile mosquito-related columns are NEGATIVE for juvenile mosquitoes
# Aka because they were field negative, they weren't ever processed
length(processed0) == length(processed0_all_null)

# Create new column for positive for juvenile mosquitoes (1 = positive for juvenile mosquitoes; 0 = negative for juvenile mosquitoes)
# All columns are 0: juvenile_pos == 1
# All columns are NA: juvenile_pos == 1
# Some columns are NA but the rest are 0: juvenile_pos == 1
# Some columns are NA but at least one non-NA column isn’t 0: juvenile_pos == 0
# At least one column isn’t equal to 0: juvenile_pos == 0

wet_CBs$juvenile_pos <- apply(wet_CBs[, juv_mos_cols], 1, function(row) {
  
  # Captures all NULL values
  if (all(is.na(row))) {
    return(1)
  }
  
  # Captures all 0 values & some NULL but all non-NULL are zero
  if (all(row[!is.na(row)] == 0)) {
    return(1)
  }
  
  # Captures at least one non-zero value & some NULL but at least one non-NULL !0
  return(0)
})

# Verify that there are no NULL values in newly created juvenile_pos column
any(is.na(wet_CBs$juvenile_pos))

#### Juvenile mosquito presence in catch basin analysis ----

# Data for juvenile mosquito presence analysis
holdingmosquitoes <- wet_CBs

## Selecting best temperature variable
tempVars <- icb_raw %>% # Create dataframe with selected temp cols
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1)
tempVars <- colnames(tempVars) # Create a vector of the temp by extracting the col names

# Declare vectors
Predictor <- vector() # Holds the variable being used for prediction
Predictor_Coef <- vector() # Holds the coefficient of the predictor variable
Model_AIC <- vector() # Holds the AIC score for the model
Model_Formula <- vector() # Holds the model formula

# Specify data
for (i in tempVars) { # Specify prcpVars or tempVars
  formula <- as.formula(paste("juvenile_pos ~", i)) # Change response variable as needed!!!!
  print(formula)
  model <- glm(formula, data = holdingmosquitoes, family = binomial(link = "logit"))
  Predictor[i] = i
  Predictor_Coef[i] = model$coefficients[2]
  Model_AIC[i] = model$aic
  Model_Formula[i] = as.character(model$formula)
}

# Create dataframe to compare models
compareModels_df <- data.frame(cbind(Predictor, Predictor_Coef, Model_AIC, Model_Formula))
nrow(compareModels_df) # Triple check all predictors accounted for

# Find the model with the smallest AIC score
smallest_AIC <- min(compareModels_df$Model_AIC) # Determine the smallest AIC score
compareModels_df[Model_AIC == smallest_AIC, ] # Select the row with the smallest AIC score


## Selecting best precipitation variable
prcpVars <- icb_raw %>%
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20)
prcpVars <- colnames(prcpVars)

# Declare vectors
Predictor <- vector() # Holds the variable being used for prediction
Predictor_Coef <- vector() # Holds the coefficient of the predictor variable
Model_AIC <- vector() # Holds the AIC score for the model
Model_Formula <- vector() # Holds the model formula

# Specify data
for (i in prcpVars) { # Specify prcpVars or tempVars
  formula <- as.formula(paste("juvenile_pos ~", i)) # Change response variable as needed!!!!
  print(formula)
  model <- glm(formula, data = holdingmosquitoes, family = binomial(link = "logit"))
  Predictor[i] = i
  Predictor_Coef[i] = model$coefficients[2]
  Model_AIC[i] = model$aic
  Model_Formula[i] = as.character(model$formula)
}

# Create dataframe to compare models
compareModels_df <- data.frame(cbind(Predictor, Predictor_Coef, Model_AIC, Model_Formula))
nrow(compareModels_df) # Triple check all predictors accounted for

# Find the model with the smallest AIC score
smallest_AIC <- min(compareModels_df$Model_AIC) # Determine the smallest AIC score
compareModels_df[Model_AIC == smallest_AIC, ] # Select the row with the smallest AIC score


## Mixed effects logistic regression model
# Response variable: juvenile_pos (1 = no mosquitoes present; 0 = mosquitoes present)

# Predictor variables:
# NSA_TMAX_C
# NS_PRCP_mm_sum_d0_d4 
# AreaCB_Rclass
# Time
# AreaCB_Rclass * Time

# Random effects
# CB_ID
# EpiweekYear

# Check data types
str(holdingmosquitoes$juvenile_pos)
str(holdingmosquitoes$NSA_TMAX_C)
str(holdingmosquitoes$NS_PRCP_mm_sum_d0_d4)
str(holdingmosquitoes$AreaCB_Rclass)
str(holdingmosquitoes$Time)
str(holdingmosquitoes$CB_ID)
str(holdingmosquitoes$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) # Comparisons to the grand mean, NOT to only one baseline (reference) level

## Models

# Full model
full_model_mos_pres <- glmmTMB(juvenile_pos ~ NSA_TMAX_C + NS_PRCP_mm_sum_d0_d4 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                               data = holdingmosquitoes,
                               family = binomial(link= "logit"))

# Reduced model 1: dropped temp var
red_model1_mos_pres <- glmmTMB(juvenile_pos ~ NS_PRCP_mm_sum_d0_d4 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                               data = holdingmosquitoes,
                               family = binomial(link= "logit"))

# Reduced model 2: dropped prcp var
red_model2_mos_pres <- glmmTMB(juvenile_pos ~ NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                               data = holdingmosquitoes,
                               family = binomial(link= "logit"))

# Reduced model 3: dropped temp and prcp vars
red_model3_mos_pres <- glmmTMB(juvenile_pos ~ AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                               data = holdingmosquitoes,
                               family = binomial(link= "logit"))

# Compare models
model.sel(full_model_mos_pres, red_model1_mos_pres, red_model2_mos_pres, red_model3_mos_pres,
          rank = AICc)

# Selected model: full model

# Summary
summary_holdingmosquitoes <- summary(full_model_mos_pres)
summary_holdingmosquitoes
exp(summary_holdingmosquitoes$coefficients$cond)
confint(full_model_mos_pres, level = 0.95)
exp(confint(full_model_mos_pres, level = 0.95))

# Anova
Anova(full_model_mos_pres, type = 3)

# Emmeans
emmeans(full_model_mos_pres, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")

# By Time by AreaCB_Rclass
emobject_mos_pres <- emmeans(full_model_mos_pres, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_mos_pres, "revpairwise", by = "Time", type = "response", adjust = "sidak")

# By AreaCB_Rclass by Time
emobject_mos_pres2 <- emmeans(full_model_mos_pres, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_mos_pres2, "revpairwise", by = "AreaCB_Rclass", type = "response", adjust = "sidak")


## Model diagnostics

# Calculate residuals
simulationOutput_mosquitoPres <- simulateResiduals(fittedModel = full_model_mos_pres, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_mosquitoPres)
plotQQunif(simulationOutput_mosquitoPres) # Left plot
plotResiduals(simulationOutput_mosquitoPres) # Right plot

#### Visualize juvenile mosquito presence ----

# Sample sizes and summary stats
holdingmosquitoes %>%
  group_by(AreaCB_Rclass, Time) %>%
  summarise(
    total_obs = n(),
    obs_NOT_holdingmosquitoes = sum(juvenile_pos == 1),
    percent_NOT_holdingmosquitoes = round(obs_NOT_holdingmosquitoes/total_obs *100, 2)
  )


# Calculate the proportion of catch basins NOT holding juvenile mosquitoes by construction period, epiweek, and type
mosquitoesfigure_data <- holdingmosquitoes %>%
  group_by(Time, EpiweekYear, AreaCB_Rclass) %>%
  summarise(
    prop_NOT_holding_mosquitoes = sum(juvenile_pos == 1, na.rm = TRUE)/sum(juvenile_pos == 1 | juvenile_pos == 0, na.rm = TRUE),
    groupSize = n(),
    .groups = "drop"
  )

mosquitoes_fig <- ggplot(mosquitoesfigure_data, aes(x = Time, y = prop_NOT_holding_mosquitoes, fill = AreaCB_Rclass, color = AreaCB_Rclass)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Construction Period", y = "Proportion of Catch Basins Not Holding\nJuvenile Mosquitoes") +
  scale_color_manual(values = c("black", "black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  scale_x_discrete(labels = c("Pre","Post")) +
  scale_y_continuous(limits = c(-0.05,1), breaks = seq(0,1,0.2)) +
  annotate("text", x = 0.7, y = -0.04, label = expression(italic(n) == 211), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = -0.04, label = expression(italic(n) == 84), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.3, y = -0.04, label = expression(italic(n) == 138), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.7, y = -0.04, label = expression(italic(n) == 397), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -0.04, label = expression(italic(n) == 196), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.3, y = -0.04, label = expression(italic(n) == 152), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme

mosquitoes_fig


#### Standing water depth analysis ----

## Data for standing water depth analysis
waterdepth <- icb_clean %>% 
  subset(Dry == 0 & Inspected == 1) %>% # Select obs where the catch basin was inspected for water presence and holding water
  drop_na(Water_Depth_cm) # Drop obs with NA water depth data (this variable wasn't recorded in early 2013)

## Selecting best temperature variable
waterDepth_tempVars <- waterdepth %>% # Create dataframe with selected temp cols and Water_Depth_cm
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1, Water_Depth_cm)

# Correlation btwn tempVars and water depth
corrs_waterDepth_tempVars <- cor(x = waterDepth_tempVars[, colnames(waterDepth_tempVars) != "Water_Depth_cm"], 
                                 y = waterDepth_tempVars$Water_Depth_cm, 
                                 method = "spearman")

corrs_waterDepth_tempVars <- as.data.frame(corrs_waterDepth_tempVars)
colnames(corrs_waterDepth_tempVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_waterDepth_tempVars$absval_r <- abs(corrs_waterDepth_tempVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_waterDepth_tempVars$absval_r, na.rm = TRUE)

# Extract the temperature variable corresponding to the greatest correlation 
corrs_waterDepth_tempVars %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable
waterDepth_prcpVars <- waterdepth %>% # Create dataframe with selected prcp cols and Water_Depth_cm
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20, Water_Depth_cm)

# Correlation btwn prcpVars and water depth
corrs_waterDepth_prcpVars <- cor(x = waterDepth_prcpVars[, colnames(waterDepth_prcpVars) != "Water_Depth_cm"], 
                                 y = waterDepth_prcpVars$Water_Depth_cm,
                                 method = "pearson")

corrs_waterDepth_prcpVars <- as.data.frame(corrs_waterDepth_prcpVars)
colnames(corrs_waterDepth_prcpVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_waterDepth_prcpVars$absval_r <- abs(corrs_waterDepth_prcpVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_waterDepth_prcpVars$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_waterDepth_prcpVars %>%
  filter(absval_r == maxAbsval_r)

## Mixed effects linear regression model
# Response variable: Water_Depth_cm

# Predictor variables:
# NSA_TMAX_C
# NS_PRCP_mm_sum_Lag10  
# AreaCB_Rclass
# Time
# AreaCB_Rclass * Time

# Random effects
# CB_ID
# EpiweekYear

# Check data types
str(waterdepth$Water_Depth_cm)
str(waterdepth$NSA_TMAX_C)
str(waterdepth$NS_PRCP_mm_sum_Lag10)
str(waterdepth$AreaCB_Rclass)
str(waterdepth$Time)
str(waterdepth$CB_ID)
str(waterdepth$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) # Comparisons to the grand mean, NOT to only one baseline (reference) level

## Models

# Full model
full_model_water_depth <- glmmTMB(Water_Depth_cm ~ NSA_TMAX_C + NS_PRCP_mm_sum_Lag10 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                  data = waterdepth,
                                  family = gaussian(link = "identity"))

# Reduced model 1: dropped temp var
red_model1_water_depth <- glmmTMB(Water_Depth_cm ~ NS_PRCP_mm_sum_Lag10 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                  data = waterdepth,
                                  family = gaussian(link = "identity"))

# Reduced model 2: dropped prcp var
red_model2_water_depth <- glmmTMB(Water_Depth_cm ~ NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                  data = waterdepth,
                                  gaussian(link = "identity"))

# Reduced model 3: dropped temp and prcp vars
red_model3_water_depth <- glmmTMB(Water_Depth_cm ~ AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                                  data = waterdepth,
                                  gaussian(link = "identity"))

# Compare models
model.sel(full_model_water_depth, red_model1_water_depth, red_model2_water_depth, red_model3_water_depth,
          rank = AICc)

# Selected model: reduced model #3

# Summary
summary(red_model3_water_depth)
confint(red_model3_water_depth, level = 0.95)

# Anova
Anova(red_model3_water_depth, type = 3)

# Emmeans
emmeans(red_model3_water_depth, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")

# By Time by AreaCB_R  class
emobject_water_depth <- emmeans(red_model3_water_depth, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_water_depth, "revpairwise", by = "Time", type = "response", adjust = "sidak")

# By AreaCB_Rclass by Time
emobject_water_depth2 <- emmeans(red_model3_water_depth, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_water_depth2, "revpairwise", by = "AreaCB_Rclass", type = "response", adjust = "sidak")

## Model diagnostics

# Calculate residuals
simulationOutput_waterDepth <- simulateResiduals(fittedModel = red_model3_water_depth, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_waterDepth)
plotQQunif(simulationOutput_waterDepth) # Left plot
plotResiduals(simulationOutput_waterDepth) # Right plot

#### Visualize water depth ----

# Sample sizes and summary stats
waterdepth %>%
  group_by(AreaCB_Rclass, Time) %>%
  summarise(
    total_obs = n(),
    avg_waterdepth = round(mean(Water_Depth_cm), 2),
    se_waterdepth = sd(Water_Depth_cm) / sqrt(total_obs)
  )

waterdepthfigure_data <- waterdepth %>%
  group_by(Time, EpiweekYear, AreaCB_Rclass) %>%
  summarise(
    avg_waterdepth = mean(Water_Depth_cm),
    groupSize = n(),
    .groups = "drop"
  )

waterdepth_fig <- ggplot(waterdepthfigure_data, aes(x = Time, y = avg_waterdepth, fill = AreaCB_Rclass, color = AreaCB_Rclass)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Construction Period", y = "Average Water Depth (cm)") +
  scale_color_manual(values = c("black", "black", "black"), labels = legend_labs, name = "Cateh Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Cateh Basin Type") +
  scale_x_discrete(labels = c("Pre","Post")) +
  annotate("text", x = 0.7, y = -0.04, label = expression(italic(n) == 173), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = -0.04, label = expression(italic(n) == 73), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.3, y = -0.04, label = expression(italic(n) == 111), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.7, y = -0.04, label = expression(italic(n) == 401), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -0.04, label = expression(italic(n) == 196), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.3, y = -0.04, label = expression(italic(n) == 152), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme + theme(legend.title.position = "top", legend.title = element_text(hjust = 0.5))

waterdepth_fig

ggsave(
  plot = waterdepth_fig,
  filename = "figures/manuscript/Exp1_waterdepth.jpeg",
  device ="jpeg",
  units = "mm",
  height = 100, width = 88, dpi = 300
)

#### Data organizing for abundance models ----

# Data for abundance models 
abundance <- wet_CBs %>%
  mutate(L1L3_recode = ifelse(is.na(L1L3_adj), 0, L1L3_adj)) %>% # Convert NAs to 0s
  mutate(L4_recode = ifelse(is.na(L4_Cx_total ), 0, L4_Cx_total)) %>% # Convert NAs to 0s
  mutate(P_recode = ifelse(is.na(P_Cx_Cx_adj), 0, P_Cx_Cx_adj)) %>% # Convert NAs to 0s
  mutate(combined_abundance = L1L3_recode + L4_recode + P_recode) # Calculate combined juvenile mosquito abundance 

#### L1L3 abundance model ----

# Data for L1L3 model
L1L3 <- abundance

## Selecting best temperature variable
L3L3_tempVars <- L1L3 %>% # Create dataframe with selected temp cols and L1L3_recode
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1, L1L3_recode)

# Correlation btwn tempVars and L1L3_recode
corrs_L1L3_tempVars <- cor(x = L3L3_tempVars[, colnames(L3L3_tempVars) != "L1L3_recode"], 
                           y = L3L3_tempVars$L1L3_recode, 
                           method = "spearman")

corrs_L1L3_tempVars <- as.data.frame(corrs_L1L3_tempVars)
colnames(corrs_L1L3_tempVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_L1L3_tempVars$absval_r <- abs(corrs_L1L3_tempVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_L1L3_tempVars$absval_r, na.rm = TRUE)

# Extract the temperature variable corresponding to the greatest correlation 
corrs_L1L3_tempVars %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable
L1L3_prcpVars <- L1L3 %>% # Create dataframe with selected prcp cols and L1L3_recode
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20, L1L3_recode)

# Correlation btwn prcpVars and L1L3_recode
corrs_L1L3_prcpVars <- cor(x = L1L3_prcpVars[, colnames(L1L3_prcpVars) != "L1L3_recode"], 
                           y = L1L3_prcpVars$L1L3_recode,
                           method = "spearman")

corrs_L1L3_prcpVars <- as.data.frame(corrs_L1L3_prcpVars)
colnames(corrs_L1L3_prcpVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_L1L3_prcpVars$absval_r <- abs(corrs_L1L3_prcpVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_L1L3_prcpVars$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_L1L3_prcpVars %>%
  filter(absval_r == maxAbsval_r)


## Examining other predictor variables

# Salinity
sum(is.na(L1L3$Salinity_ppm)) # 38 NA obs

NAsalinity <- subset(L1L3, is.na(Salinity_ppm))

# Relative proportion of water covered by floating OM

# Possible values
unique(L1L3$Surf_Detrit)
sum(is.na(L1L3$Surf_Detrit)) # 8 missing obs

# Boxplot
ggplot(L1L3, aes(x = Surf_Detrit, y = L1L3_recode, fill = Surf_Detrit)) +
  geom_boxplot()

L1L3 %>%
  group_by(Surf_Detrit) %>%
  summarise(
    avg_L1L3 = mean(L1L3_recode),
    var_L1L3 = var(L1L3_recode),
    groupSize = n())

# ANOVA -- response: L1L3 abundance, groups: H, M, L %OM

# Kruskal-Wallis 
kruskal.test(L1L3_recode ~ Surf_Detrit, data = L1L3)

## Negative binomial glmm

# Filter out any missing data for predictors
L1L3 <- L1L3 %>%
  filter(!is.na(Surf_Detrit)) %>%
  filter(!is.na(Salinity_ppm))

# Response: L1L3_recode

# Predictor variables:
# NS_PRCP_mm_sum_d0_d4
# NSA_TMAX_C
# Salinity_ppm
# Surf_Detrit
# AreaCB_Rclass
# Time
# AreaCB_Rclass * Time

# Random effects:
# CB_ID
# EpiweekYear

# Check data types
str(L1L3$L1L3_recode)
str(L1L3$NS_PRCP_mm_sum_d0_d4)
str(L1L3$NSA_TMAX_C)
str(L1L3$Salinity_ppm)
str(L1L3$Surf_Detrit)
str(L1L3$AreaCB_Rclass)
str(L1L3$Time)
str(L1L3$CB_ID)
str(L1L3$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) # Comparisons to the grand mean, NOT to only one baseline (reference) level

## Models

full_model_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 1: removed temp var
red_model1_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 2: removed prcp var
red_model2_L1L3 <- glmmTMB(L1L3_recode ~  NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 3: removed salinity var
red_model3_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 4: removed surface detritus var
red_model4_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 5: removed temp and prcp vars
red_model5_L1L3 <- glmmTMB(L1L3_recode ~ Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 6: removed temp and salinity vars
red_model6_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 7: removed temp and surface detritus variables
red_model7_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 8: removed prcp and salinity vars
red_model8_L1L3 <- glmmTMB(L1L3_recode ~  NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 9: removed prcp and surface detritus vars
red_model9_L1L3 <- glmmTMB(L1L3_recode ~  NSA_TMAX_C + Salinity_ppm  + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 10: removed salinity and surface detritus variables
red_model10_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                            data = L1L3,
                            family = nbinom2(link = "log"))

# Reduced model 11: removed prcp, salinity, and surface detritus vars
red_model11_L1L3 <- glmmTMB(L1L3_recode ~  NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                            data = L1L3,
                            family = nbinom2(link = "log"))

# Reduced model 12: removed temp, salinity, and surface detritus vars
red_model12_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d4 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                            data = L1L3,
                            family = nbinom2(link = "log"))

# Reduced model 13: removed temp, prcp, and surface detritus vars
red_model13_L1L3 <- glmmTMB(L1L3_recode ~ Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                            data = L1L3,
                            family = nbinom2(link = "log"))

# Reduced model 14: removed temp, prcp, and salinity vars
red_model14_L1L3 <- glmmTMB(L1L3_recode ~ Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                            data = L1L3,
                            family = nbinom2(link = "log"))

# Reduced model 15: removed, temp, prcp, salinity, and surface detritus variables
red_model15_L1L3 <- glmmTMB(L1L3_recode ~ AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                            data = L1L3,
                            family = nbinom2(link = "log"))

# Compare models
model.sel(full_model_L1L3, red_model1_L1L3, red_model2_L1L3, red_model3_L1L3, red_model4_L1L3, red_model5_L1L3, red_model6_L1L3, red_model7_L1L3, 
          red_model8_L1L3, red_model9_L1L3, red_model10_L1L3, red_model12_L1L3, red_model13_L1L3, red_model14_L1L3, red_model15_L1L3,
          rank = AICc)


# Selected model: reduced model 3

# Summary
summary(red_model3_L1L3)

confint(red_model3_L1L3, level = 0.95)

# Anova
Anova(red_model3_L1L3, type = 3)

# Emmeans
emmeans(red_model3_L1L3, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")

# By Time by AreaCB_R  class
emobject_L1L3 <- emmeans(red_model3_L1L3, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_L1L3, "revpairwise", by = "Time", type = "response", adjust = "sidak")

# By AreaCB_Rclass by Time
emobject_L1L3_2 <- emmeans(red_model3_L1L3, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_L1L3_2, "revpairwise", by = "AreaCB_Rclass", type = "response", adjust = "sidak")

## Model diagnostics

# Calculate residuals
simulationOutput_L1L3 <- simulateResiduals(fittedModel = red_model3_L1L3, plot = FALSE)

# Test zero inflation
testZeroInflation(simulationOutput_L1L3)

# Plotting the scaled residuals
plot(simulationOutput_L1L3)
plotQQunif(simulationOutput_L1L3) # Left plot
plotResiduals(simulationOutput_L1L3) # Right plot

#### L4 abundance model ----

# Data for L4 model
L4 <- abundance

## Selecting best temperature variable
L4_tempVars <- L4 %>% # Create dataframe with selected temp cols and L4_recode
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1, L4_recode)

# Correlation btwn tempVars and L4_recode
corrs_L4_tempVars <- cor(x = L4_tempVars[, colnames(L4_tempVars) != "L4_recode"], 
                         y = L4_tempVars$L4_recode, 
                         method = "spearman")

corrs_L4_tempVars <- as.data.frame(corrs_L4_tempVars)
colnames(corrs_L4_tempVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_L4_tempVars$absval_r <- abs(corrs_L4_tempVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_L4_tempVars$absval_r, na.rm = TRUE)

# Extract the temperature variable corresponding to the greatest correlation 
corrs_L4_tempVars %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable
L4_prcpVars <- L4 %>% # Create dataframe with selected prcp cols and L4_recode
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20, L4_recode)

# Correlation btwn prcpVars and L4_recode
corrs_L4_prcpVars <- cor(x = L4_prcpVars[, colnames(L4_prcpVars) != "L4_recode"], 
                         y = L4_prcpVars$L4_recode,
                         method = "spearman")

corrs_L4_prcpVars <- as.data.frame(corrs_L4_prcpVars)
colnames(corrs_L4_prcpVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_L4_prcpVars$absval_r <- abs(corrs_L4_prcpVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_L4_prcpVars$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_L4_prcpVars %>%
  filter(absval_r == maxAbsval_r)


## Examining other predictor variables

# Salinity
sum(is.na(L4$Salinity_ppm)) # 38 NA obs

NAsalinity <- subset(L4, is.na(Salinity_ppm))

# Relative proportion of water covered by floating OM

# Possible values
unique(L4$Surf_Detrit)
sum(is.na(L4$Surf_Detrit)) # 8 missing obs

# Boxplot
ggplot(L4, aes(x = Surf_Detrit, y = L4_recode, fill = Surf_Detrit)) +
  geom_boxplot()

L4 %>%
  group_by(Surf_Detrit) %>%
  summarise(
    avg_L4 = mean(L4_recode),
    var_L4 = var(L4_recode),
    groupSize = n())

# ANOVA -- response: L4 abundance, groups: H, M, L %OM

# Kruskal-Wallis 
kruskal.test(L4_recode ~ Surf_Detrit, data = L4)


## Zero-inflated model

# Filter out any missing data for predictors
L4 <- L4 %>%
  filter(!is.na(Surf_Detrit)) %>%
  filter(!is.na(Salinity_ppm))

# Response: L4_recode

# Predictor variables:
# NS_PRCP_mm_sum_d0_d7 
# NSA_TMAX_C
# Salinity_ppm
# Surf_Detrit
# AreaCB_Rclass
# Time
# AreaCB_Rclass * Time

# Random effects:
# CB_ID
# EpiweekYear


# Check data types
str(L4$L4_recode)
str(L4$NS_PRCP_mm_sum_d0_d7)
str(L4$NSA_TMAX_C)
str(L4$Salinity_ppm)
str(L4$Surf_Detrit)
str(L4$AreaCB_Rclass)
str(L4$Time)
str(L4$CB_ID)
str(L4$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) # Comparisons to the grand mean, NOT to only one baseline (reference) level

## Models

full_model_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 1: removed temp var
red_model1_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 2: removed prcp var
red_model2_L4 <- glmmTMB(L4_recode ~  NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 3: removed salinity var
red_model3_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 4: removed surface detritus var
red_model4_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 5: removed temp and prcp vars
red_model5_L4 <- glmmTMB(L4_recode ~ Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 6: removed temp and salinity vars
red_model6_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 7: removed temp and surface detritus variables
red_model7_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 8: removed prcp and salinity vars
red_model8_L4 <- glmmTMB(L4_recode ~  NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 9: removed prcp and surface detritus vars
red_model9_L4 <- glmmTMB(L4_recode ~  NSA_TMAX_C + Salinity_ppm  + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 10: removed salinity and surface detritus variables
red_model10_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = L4,
                          family = nbinom2(link = "log"))

# Reduced model 11: removed prcp, salinity, and surface detritus vars
red_model11_L4 <- glmmTMB(L4_recode ~  NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = L4,
                          family = nbinom2(link = "log"))

# Reduced model 12: removed temp, salinity, and surface detritus vars
red_model12_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d0_d7 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = L4,
                          family = nbinom2(link = "log"))

# Reduced model 13: removed temp, prcp, and surface detritus vars
red_model13_L4 <- glmmTMB(L4_recode ~ Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = L4,
                          family = nbinom2(link = "log"))

# Reduced model 14: removed temp, prcp, and salinity vars
red_model14_L4 <- glmmTMB(L4_recode ~ Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = L4,
                          family = nbinom2(link = "log"))

# Reduced model 15: removed, temp, prcp, salinity, and surface detritus variables
red_model15_L4 <- glmmTMB(L4_recode ~ AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = L4,
                          family = nbinom2(link = "log"))

# Compare models
model.sel(full_model_L4, red_model1_L4, red_model2_L4, red_model3_L4, red_model4_L4, red_model5_L4, red_model6_L4, red_model7_L4, 
          red_model8_L4, red_model9_L4, red_model10_L4, red_model12_L4, red_model13_L4, red_model14_L4, red_model15_L4,
          rank = AICc)


# Selected model: reduced model 1

# Summary
summary(red_model1_L4)
confint(red_model1_L4, level = 0.95)

# Anova
Anova(red_model1_L4, type = 3)

# Emmeans
emmeans(red_model1_L4, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")

# By Time by AreaCB_R  class
emobject_L4 <- emmeans(red_model1_L4, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_L4, "revpairwise", by = "Time", type = "response", adjust = "sidak")

# By AreaCB_Rclass by Time
emobject_L4_2 <- emmeans(red_model1_L4, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_L4_2, "revpairwise", by = "AreaCB_Rclass", type = "response", adjust = "sidak")

## Model diagnostics

# Calculate residuals
simulationOutput_L4 <- simulateResiduals(fittedModel = red_model1_L4, plot = FALSE)

# Test zero inflation
testZeroInflation(simulationOutput_L4)

# Plotting the scaled residuals
plot(simulationOutput_L4)
plotQQunif(simulationOutput_L4) # Left plot
plotResiduals(simulationOutput_L4) # Right plot

#### Pupal abundance model ----

# Data for pupal model
pupae <- abundance

## Selecting best temperature variable
P_tempVars <- pupae %>% # Create dataframe with selected temp cols and P_recode
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1, P_recode)

# Correlation btwn tempVars and P_recode
corrs_P_tempVars <- cor(x = P_tempVars[, colnames(P_tempVars) != "P_recode"], 
                        y = P_tempVars$P_recode, 
                        method = "spearman")

corrs_P_tempVars <- as.data.frame(corrs_P_tempVars)
colnames(corrs_P_tempVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_P_tempVars$absval_r <- abs(corrs_P_tempVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_P_tempVars$absval_r, na.rm = TRUE)

# Extract the temperature variable corresponding to the greatest correlation 
corrs_P_tempVars %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable
P_prcpVars <- pupae %>% # Create dataframe with selected prcp cols and P_recode
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20, P_recode)

# Correlation btwn prcpVars and P_recode
corrs_P_prcpVars <- cor(x = P_prcpVars[, colnames(P_prcpVars) != "P_recode"], 
                        y = P_prcpVars$P_recode,
                        method = "spearman")

corrs_P_prcpVars <- as.data.frame(corrs_P_prcpVars)
colnames(corrs_P_prcpVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_P_prcpVars$absval_r <- abs(corrs_P_prcpVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_P_prcpVars$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_P_prcpVars %>%
  filter(absval_r == maxAbsval_r)


## Examining other predictor variables

# Salinity
sum(is.na(pupae$Salinity_ppm)) # 38 NA obs

NAsalinity <- subset(pupae, is.na(Salinity_ppm))

# Relative proportion of water covered by floating OM

# Possible values
unique(pupae$Surf_Detrit)
sum(is.na(pupae$Surf_Detrit)) # 8 missing obs

# Boxplot
ggplot(pupae, aes(x = Surf_Detrit, y = P_recode, fill = Surf_Detrit)) +
  geom_boxplot()

pupae %>%
  group_by(Surf_Detrit) %>%
  summarise(
    avg_P = mean(P_recode),
    var_P = var(P_recode),
    groupSize = n())

# ANOVA -- response: Pupal abundance, groups: H, M, L %OM

# Kruskal-Wallis 
kruskal.test(P_recode ~ Surf_Detrit, data = pupae)


## Zero-inflated model

# Filter out any missing data for predictors
pupae <- pupae %>%
  filter(!is.na(Surf_Detrit)) %>%
  filter(!is.na(Salinity_ppm))

# Response: P_recode

# Predictor variables:
# NS_PRCP_mm_sum_d0_d7 
# NSA_TMAX_C
# Salinity_ppm
# Surf_Detrit
# AreaCB_Rclass
# Time
# AreaCB_Rclass * Time

# Random effects:
# CB_ID
# EpiweekYear

# Check data types
str(pupae$P_recode)
str(pupae$NS_PRCP_mm_sum_d0_d7)
str(pupae$NSA_TMAX_C)
str(pupae$Salinity_ppm)
str(pupae$Surf_Detrit)
str(pupae$AreaCB_Rclass)
str(pupae$Time)
str(pupae$CB_ID)
str(pupae$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) # Comparisons to the grand mean, NOT to only one baseline (reference) level

## Models

full_model_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 1: removed temp var
red_model1_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 2: removed prcp var
red_model2_P <- glmmTMB(P_recode ~  NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 3: removed salinity var
red_model3_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 4: removed surface detritus var
red_model4_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 5: removed temp and prcp vars
red_model5_P <- glmmTMB(P_recode ~ Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 6: removed temp and salinity vars
red_model6_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 7: removed temp and surface detritus variables
red_model7_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 8: removed prcp and salinity vars
red_model8_P <- glmmTMB(P_recode ~  NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 9: removed prcp and surface detritus vars
red_model9_P <- glmmTMB(P_recode ~  NSA_TMAX_C + Salinity_ppm  + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 10: removed salinity and surface detritus variables
red_model10_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = pupae,
                         family = nbinom2(link = "log"))

# Reduced model 11: removed prcp, salinity, and surface detritus vars
red_model11_P <- glmmTMB(P_recode ~  NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = pupae,
                         family = nbinom2(link = "log"))

# Reduced model 12: removed temp, salinity, and surface detritus vars
red_model12_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d0_d7 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = pupae,
                         family = nbinom2(link = "log"))

# Reduced model 13: removed temp, prcp, and surface detritus vars
red_model13_P <- glmmTMB(P_recode ~ Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = pupae,
                         family = nbinom2(link = "log"))

# Reduced model 14: removed temp, prcp, and salinity vars
red_model14_P <- glmmTMB(P_recode ~ Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = pupae,
                         family = nbinom2(link = "log"))

# Reduced model 15: removed, temp, prcp, salinity, and surface detritus variables
red_model15_P <- glmmTMB(P_recode ~ AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                         data = pupae,
                         family = nbinom2(link = "log"))

# Compare models
model.sel(full_model_P, red_model1_P, red_model2_P, red_model3_P, red_model4_P, red_model5_P, red_model6_P, red_model7_P, 
          red_model8_P, red_model9_P, red_model10_P, red_model12_P, red_model13_P, red_model14_P, red_model15_P,
          rank = AICc)


# Selected model: reduced model 6

# Summary
summary(red_model6_P)
confint(red_model6_P, level = 0.95)

# Anova
Anova(red_model6_P, type = 3)

# Emmeans
emmeans(red_model6_P, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")

# By Time by AreaCB_R  class
emobject_P <- emmeans(red_model6_P, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_P, "revpairwise", by = "Time", type = "response", adjust = "sidak")

# By AreaCB_Rclass by Time
emobject_P_2 <- emmeans(red_model6_P, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_P_2, "revpairwise", by = "AreaCB_Rclass", type = "response", adjust = "sidak")

## Model diagnostics

# Calculate residuals
simulationOutput_P <- simulateResiduals(fittedModel = red_model6_P, plot = FALSE)

# Test zero inflation
testZeroInflation(simulationOutput_P)

# Plotting the scaled residuals
plot(simulationOutput_P)
plotQQunif(simulationOutput_P) # Left plot
plotResiduals(simulationOutput_P) # Right plot

#### Combined juvenile mosquito abundance model ----

# Data for combined model
combined <- abundance 

## Selecting best temperature variable
all_tempVars <- combined %>% # Create dataframe with selected temp cols and combined_abundance
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1, combined_abundance)

# Correlation btwn tempVars and combined_abundance
corrs_all_tempVars <- cor(x = all_tempVars[, colnames(all_tempVars) != "combined_abundance"], 
                          y = all_tempVars$combined_abundance, 
                          method = "spearman")

corrs_all_tempVars <- as.data.frame(corrs_all_tempVars)
colnames(corrs_all_tempVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_all_tempVars$absval_r <- abs(corrs_all_tempVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_all_tempVars$absval_r, na.rm = TRUE)

# Extract the temperature variable corresponding to the greatest correlation 
corrs_all_tempVars %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable
all_prcpVars <- combined %>% # Create dataframe with selected prcp cols and combined_abundance
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20, combined_abundance)

# Correlation btwn prcpVars and combined_abundance
corrs_all_prcpVars <- cor(x = all_prcpVars[, colnames(all_prcpVars) != "combined_abundance"], 
                          y = all_prcpVars$combined_abundance,
                          method = "spearman")

corrs_all_prcpVars <- as.data.frame(corrs_all_prcpVars)
colnames(corrs_all_prcpVars)[1] <- "r"

# Make new column with absolute value of r values
corrs_all_prcpVars$absval_r <- abs(corrs_all_prcpVars$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_all_prcpVars$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_all_prcpVars %>%
  filter(absval_r == maxAbsval_r)


## Examining other predictor variables

# Salinity
sum(is.na(combined$Salinity_ppm)) # 38 NA obs

NAsalinity <- subset(combined, is.na(Salinity_ppm))

# Relative proportion of water covered by floating OM

# Possible values
unique(combined$Surf_Detrit)
sum(is.na(combined$Surf_Detrit)) # 8 missing obs

# Boxplot
ggplot(combined, aes(x = Surf_Detrit, y = combined_abundance, fill = Surf_Detrit)) +
  geom_boxplot()

# ANOVA -- response: Pupal abundance, groups: H, M, L %OM

# Kruskal-Wallis 
kruskal.test(combined_abundance ~ Surf_Detrit, data = combined)


## Negative binomial count model

# Filter out any missing data for predictors
combined <- combined %>%
  filter(!is.na(Surf_Detrit)) %>%
  filter(!is.na(Salinity_ppm))

# Response: combined_abundance

# Predictor variables:
# NS_PRCP_mm_sum_d0_d4  
# NSA_TMAX_C
# Salinity_ppm
# Surf_Detrit
# AreaCB_Rclass
# Time
# AreaCB_Rclass * Time

# Random effects:
# CB_ID
# EpiweekYear

# Check data types
str(combined$combined_abundance)
str(combined$NS_PRCP_mm_sum_d0_d4)
str(combined$NSA_TMAX_C)
str(combined$Salinity_ppm)
str(combined$Surf_Detrit)
str(combined$AreaCB_Rclass)
str(combined$Time)
str(combined$CB_ID)
str(combined$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) # Comparisons to the grand mean, NOT to only one baseline (reference) level

## Models

full_model_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 1: removed temp var
red_model1_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 2: removed prcp var
red_model2_all <- glmmTMB(combined_abundance ~  NSA_TMAX_C + Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 3: removed salinity var
red_model3_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 4: removed surface detritus var
red_model4_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 5: removed temp and prcp vars
red_model5_all <- glmmTMB(combined_abundance ~ Salinity_ppm + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 6: removed temp and salinity vars
red_model6_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 7: removed temp and surface detritus variables
red_model7_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 8: removed prcp and salinity vars
red_model8_all <- glmmTMB(combined_abundance ~  NSA_TMAX_C + Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 9: removed prcp and surface detritus vars
red_model9_all <- glmmTMB(combined_abundance ~  NSA_TMAX_C + Salinity_ppm  + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                          data = combined,
                          family = nbinom2(link = "log"))

# Reduced model 10: removed salinity and surface detritus variables
red_model10_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = combined,
                           family = nbinom2(link = "log"))

# Reduced model 11: removed prcp, salinity, and surface detritus vars
red_model11_all <- glmmTMB(combined_abundance ~  NSA_TMAX_C + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = combined,
                           family = nbinom2(link = "log"))

# Reduced model 12: removed temp, salinity, and surface detritus vars
red_model12_all <- glmmTMB(combined_abundance ~ NS_PRCP_mm_sum_d0_d4 + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = combined,
                           family = nbinom2(link = "log"))

# Reduced model 13: removed temp, prcp, and surface detritus vars
red_model13_all <- glmmTMB(combined_abundance ~ Salinity_ppm + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = combined,
                           family = nbinom2(link = "log"))

# Reduced model 14: removed temp, prcp, and salinity vars
red_model14_all <- glmmTMB(combined_abundance ~ Surf_Detrit + AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = combined,
                           family = nbinom2(link = "log"))

# Reduced model 15: removed, temp, prcp, salinity, and surface detritus variables
red_model15_all <- glmmTMB(combined_abundance ~ AreaCB_Rclass + Time + AreaCB_Rclass*Time + (1|CB_ID) + (1|EpiweekYear),
                           data = combined,
                           family = nbinom2(link = "log"))

# Compare models
model.sel(full_model_all, red_model1_all, red_model2_all, red_model3_all, red_model4_all, red_model5_all, red_model6_all, red_model7_all, 
          red_model8_all, red_model9_all, red_model10_all, red_model12_all, red_model13_all, red_model14_all, red_model15_all,
          rank = AICc)


# Selected model: full model

# Summary
summary(full_model_all)
confint(full_model_all, level = 0.95)

# Anova
Anova(full_model_all, type = 3)

# Emmeans
emmeans(full_model_all, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")

# By Time by AreaCB_R  class
emobject_all <- emmeans(full_model_all, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_all, "revpairwise", by = "Time", type = "response", adjust = "sidak")

# By AreaCB_Rclass by Time
emobject_all_2 <- emmeans(full_model_all, ~ AreaCB_Rclass*Time, type = "response", adjust = "sidak")
contrast(emobject_all_2, "revpairwise", by = "AreaCB_Rclass", type = "response", adjust = "sidak")

## Model diagnostics

# Calculate residuals
simulationOutput_all <- simulateResiduals(fittedModel = full_model_all, plot = FALSE)

# Check zero inflation
testZeroInflation(simulationOutput_all)

# Plotting the scaled residuals
plot(simulationOutput_all)
plotQQunif(simulationOutput_all) # Left plot
plotResiduals(simulationOutput_all) # Right plot


#### Visualize abundance results ----

# Sample sizes and summary stats
combined %>%
  group_by(AreaCB_Rclass, Time) %>%
  summarise(
    total_obs = n(),
    avg_combined = round(mean(combined_abundance), 2),
    se_combined = sd(combined_abundance) / sqrt(total_obs)
  )

L1L3 %>%
  group_by(AreaCB_Rclass, Time) %>%
  summarise(
    total_obs = n(),
    avg_L1L3 = round(mean(L1L3_recode), 2),
    se_L1L3 = sd(L1L3_recode) / sqrt(total_obs)
  )

L4 %>%
  group_by(AreaCB_Rclass, Time) %>%
  summarise(
    total_obs = n(),
    avg_L4 = round(mean(L4_recode), 2),
    se_L4 = sd(L4_recode) / sqrt(total_obs)
  )

pupae %>%
  group_by(AreaCB_Rclass, Time) %>%
  summarise(
    total_obs = n(),
    avg_P = round(mean(P_recode), 2),
    se_P = sd(P_recode) / sqrt(total_obs)
  )


#### Average Combined Juvenile Mosquito Abundance by Epiweek ----

combinedfigure_data <- combined %>%
  group_by(Time, EpiweekYear, AreaCB_Rclass) %>%
  summarise(
    avg_combined = mean(combined_abundance),
    groupSize = n(),
    .groups = "drop"
  )

combined_fig <- ggplot(combinedfigure_data, aes(x = Time, y = avg_combined, fill = AreaCB_Rclass, color = AreaCB_Rclass)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Construction Period", y = "Average Combined Juvenile Mosquito Abundance") +
  scale_color_manual(values = c("black", "black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  scale_x_discrete(labels = c("Pre","Post")) +
  scale_y_continuous(breaks = seq(0, 2000, 250), limits = c(-50, 1900)) +
  annotate("text", x = 0.7, y = -50, label = expression(italic(n) == 194), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = -50, label = expression(italic(n) == 82), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.3, y = -50, label = expression(italic(n) == 138), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.7, y = -50, label = expression(italic(n) == 386), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -50, label = expression(italic(n) == 193), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.3, y = -50, label = expression(italic(n) == 147), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme

combined_fig


#### Average L1L3 Abundance by Epiweek ----

L1L3figure_data <- L1L3 %>%
  group_by(Time, EpiweekYear, AreaCB_Rclass) %>%
  summarise(
    avg_L1L3 = mean(L1L3_recode),
    groupSize = n(),
    .groups = "drop"
  )


L1L3_fig <- ggplot(L1L3figure_data, aes(x = Time, y = avg_L1L3, fill = AreaCB_Rclass, color = AreaCB_Rclass)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Construction Period", y = "Average 1st to 3rd Instar Larval Abundance") +
  scale_color_manual(values = c("black", "black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  scale_x_discrete(labels = c("Pre","Post")) +
  scale_y_continuous(breaks = seq(0, 2000, 250)) +
  annotate("text", x = 0.7, y = -60, label = expression(italic(n) == 194), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = -60, label = expression(italic(n) == 82), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.3, y = -60, label = expression(italic(n) == 138), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.7, y = -60, label = expression(italic(n) == 386), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -60, label = expression(italic(n) == 193), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.3, y = -60, label = expression(italic(n) == 147), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme + theme(legend.title.position = "top", legend.title = element_text(hjust = 0.5))

L1L3_fig

ggsave(
  plot = L1L3_fig,
  filename = "figures/manuscript/Exp1_L1L3abundance.jpeg",
  device ="jpeg",
  units = "mm",
  height = 100, width = 88, dpi = 300
)



#### Average L4 Abundance by Epiweek ----

L4figure_data <- L4 %>%
  group_by(Time, EpiweekYear, AreaCB_Rclass) %>%
  summarise(
    avg_L4 = mean(L4_recode),
    groupSize = n(),
    .groups = "drop"
  )


L4_fig <- ggplot(L4figure_data, aes(x = Time, y = avg_L4, fill = AreaCB_Rclass, color = AreaCB_Rclass)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Construction Period", y = "Average 4th Instar Larval Abundance") +
  scale_color_manual(values = c("black", "black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, "Catch Basin Type") +
  scale_x_discrete(labels = c("Pre","Post")) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  annotate("text", x = 0.7, y = -10, label = expression(italic(n) == 194), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = -10, label = expression(italic(n) == 82), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.3, y = -10, label = expression(italic(n) == 138), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.7, y = -10, label = expression(italic(n) == 386), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -10, label = expression(italic(n) == 193), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.3, y = -10, label = expression(italic(n) == 147), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme

L4_fig


#### Average Pupae Abundance by Epiweek ----

pupaefigure_data <- pupae %>%
  group_by(Time, EpiweekYear, AreaCB_Rclass) %>%
  summarise(
    avg_pupae = mean(P_recode),
    groupSize = n(),
    .groups = "drop"
  )

pupaefigure_outlier_removed <- pupaefigure_data %>%
  filter(avg_pupae < 100)


pupae_fig_outlierremoved <- ggplot(pupaefigure_outlier_removed, aes(x = Time, y = avg_pupae, fill = AreaCB_Rclass, color = AreaCB_Rclass)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Construction Period", y = "Average Pupal Abundance") +
  scale_color_manual(values = c("black", "black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  scale_x_discrete(labels = c("Pre","Post")) +
  #scale_y_continuous(breaks = seq(0, 50, 10)) +
  annotate("text", x = 0.7, y = -3, label = expression(italic(n) == 194), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = -3, label = expression(italic(n) == 82), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.3, y = -3, label = expression(italic(n) == 138), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1.7, y = -3, label = expression(italic(n) == 386), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -3, label = expression(italic(n) == 193), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.3, y = -3, label = expression(italic(n) == 147), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme + theme(legend.title.position = "top", legend.title = element_text(hjust = 0.5))

pupae_fig_outlierremoved

ggsave(
  plot = pupae_fig_outlierremoved,
  filename = "figures/manuscript/Exp1_pupaeabundance_no_outlier.jpeg",
  device ="jpeg",
  units = "mm",
  height = 100, width = 88, dpi = 300
)


#### Combination figures ----

Fig3 <- holdingwater_fig + mosquitoes_fig + combined_fig + L4_fig +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &#, caption = "Treatment ICBs refers to catch basins that were CCBs in the pre-construction period and ICBs in the post-construction period.") &
  theme(legend.position = "bottom", 
        plot.tag.position = c(0.02, 1.01), 
        plot.tag = element_text(size = 9, face = "bold", family = "HelveticaNeueforSAS"),
        plot.caption = element_text(size = 8, face = "italic", family = "HelveticaNeueforSAS", hjust = 0),
        plot.caption.position = "plot")
Fig3

#Fig3 <- holdingwater_fig + mosquitoes_fig + combined_fig + L4_fig +
#  plot_layout(guides = "collect", axis_titles = "collect") +
#  plot_annotation(tag_levels = "a") &
#  theme(legend.position = "bottom", 
#       plot.tag.position = c(0.02, 1.01), 
#       plot.tag = element_text(size = 9, face = "bold", family = "HelveticaNeueforSAS"),
#       plot.caption.position = "plot")
#Fig3

ggsave(
  plot = Fig3,
  filename = "figures/manuscript/Figure3.jpeg",
  device ="jpeg",
  units = "mm",
  height = 185, width = 180, dpi = 300
)


#### Compare combined juvenile abundance among treatment CCBs in the pre- and post-construction periods ----

# Subset CCBs in the treatment area
trt_CCBs <- abundance %>%
  subset(AreaCB_Rclass == "Treatment-C-C")

trt_CCBs$Time <- factor(trt_CCBs$Time, levels = c("Pre ICB Construction", "Post ICB Construction"))

n_distinct(trt_CCBs$CB_ID) # 17 unique CCBs in the treatment area

# Average combined abundance values by CB ID and Time (CB 527 wasn't holding water in the pre-construction period - excluded from analysis)
trt_CBs_abundance_pre_post <- trt_CCBs %>%
  group_by(CB_ID, Time) %>%
  summarise(avg_combined_abundance = mean(combined_abundance)) %>%
  pivot_wider(
    names_from = Time,
    values_from = avg_combined_abundance
  )

# Adjsut col names
colnames(trt_CBs_abundance_pre_post) <- c("CB_ID", "Pre", "Post")

# Checking normality - pre-construction 
shapiro.test(trt_CBs_abundance_pre_post$Pre)
qqnorm(trt_CBs_abundance_pre_post$Pre)
qqline(trt_CBs_abundance_pre_post$Pre)

# Checking normality - post-construction 
shapiro.test(trt_CBs_abundance_pre_post$Post)
qqnorm(trt_CBs_abundance_pre_post$Post)
qqline(trt_CBs_abundance_pre_post$Post)

# Paired Wilcoxon test
wilcox.test(x = trt_CBs_abundance_pre_post$Pre, y = trt_CBs_abundance_pre_post$Post,
            alternative = "two.sided", paired = TRUE)

median(trt_CBs_abundance_pre_post$Pre, na.rm = TRUE)
median(trt_CBs_abundance_pre_post$Post, na.rm = TRUE)

# Convert data into long form for plotting
trt_CBs_abundance_pre_post_PLOTTING <- trt_CCBs %>%
  group_by(CB_ID, Time) %>%
  summarise(avg_combined_abundance = mean(combined_abundance))

# Change data types
trt_CBs_abundance_pre_post_PLOTTING$CB_ID <- factor(trt_CBs_abundance_pre_post_PLOTTING$CB_ID)
trt_CBs_abundance_pre_post_PLOTTING$Time <- factor(trt_CBs_abundance_pre_post_PLOTTING$Time,
                                                   levels = c("Pre ICB Construction", "Post ICB Construction"))

# Histograms
ggplot(trt_CBs_abundance_pre_post_PLOTTING, aes(x = avg_combined_abundance)) +
  geom_histogram() +
  facet_wrap(vars(Time), nrow = 2) +
  scale_y_continuous(expand = c(0,0)) +
  figtheme


# Plot
paired_trt_CCBs_abundance <- ggplot(trt_CBs_abundance_pre_post_PLOTTING, aes(x = Time, y = avg_combined_abundance, group = CB_ID)) +
  geom_line(linewidth = 0.25) +
  geom_point(size = 0.75) +
  scale_x_discrete(labels = c("Pre", "Post")) +
  scale_y_continuous(limits = c(0,1325), breaks = seq(0, 1500, 250)) +
  labs(x = "Construction Period", y = "Average Combined Juvenile Mosquito\nAbundance by Catch Basin ID") +
  annotate("text", x = 2.3, y = 0, label = expression(italic(n) == 16), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme

paired_trt_CCBs_abundance

ggsave(
  plot = paired_trt_CCBs_abundance,
  filename = "figures/manuscript/Exp1_paired_trt_CCBs_combinedabundance.jpeg",
  device ="jpeg",
  units = "mm",
  height = 100, width = 88, dpi = 300
)


