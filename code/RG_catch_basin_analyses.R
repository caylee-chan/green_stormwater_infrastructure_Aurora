# Experiment #2: RG analyses (catch basins)
# Caylee Chan
# Created: 2 Jan 2025
# Updated: 22 Jan 2026
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
legend_colors <- c("grey27", "grey60", "green4")

# Set N size size
N_size_size <- 2

#### Clean raw data ----

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

#### Standing water presence in all catch basins analysis ----

## Data organizing
standingwater_all <- rg_clean %>%
  subset(Inspected == 1) %>% # Subset observations where catch basin was inspected for the presence of standing water
  drop_na(Dry) # Drop observations where Dry is NA

## Selecting best temperature variable
tempVars <- rg_raw %>% # Create dataframe with selected temp cols
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
  model <- glm(formula, data = standingwater_all, family = binomial(link = "logit"))
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
prcpVars <- rg_raw %>%
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
  model <- glm(formula, data = standingwater_all, family = binomial(link = "logit"))
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
# NS_PRCP_mm_sum_Lag7  
# interTypecbType

# Random effects
# (Inter_ID/CB_ID_RG)
# EpiweekYear

# Check data types
str(standingwater_all$Dry)
str(standingwater_all$NSA_TMIN_C)
str(standingwater_all$NS_PRCP_mm_sum_Lag7)
str(standingwater_all$interTypecbType)
unique(standingwater_all$interTypecbType)
str(standingwater_all$Inter_ID)
str(standingwater_all$CB_ID_RG)
str(standingwater_all$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models

# Full model
full_model_water_pres <- glmmTMB(Dry ~ NSA_TMIN_C + NS_PRCP_mm_sum_Lag7 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                                 data = standingwater_all,
                                 family = binomial(link= "logit"))

# Reduced model #1: dropped temp variable
red_model1_water_pres <- glmmTMB(Dry ~ NS_PRCP_mm_sum_Lag7 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                                 data = standingwater_all,
                                 family = binomial(link= "logit"))

# Reduced model #2: dropped prcp variable
red_model2_water_pres <- glmmTMB(Dry ~ NSA_TMIN_C + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                                 data = standingwater_all,
                                 family = binomial(link= "logit"))

# Reduced model #3: dropped temp and prcp variables
red_model3_water_pres <- glmmTMB(Dry ~ interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                                 data = standingwater_all,
                                 family = binomial(link= "logit"))

# Model selection
model.sel(full_model_water_pres, red_model1_water_pres, red_model2_water_pres, red_model3_water_pres,
          rank = AICc)


# Selected model: reduced model #3

# Summary
summary_all_CBs_holdingwater <- summary(red_model3_water_pres)
summary_all_CBs_holdingwater

exp(summary_all_CBs_holdingwater$coefficients$cond)

exp(confint(red_model3_water_pres, level = 0.95))

# Anova
Anova(red_model3_water_pres, type = 3)

# Emmeans

emobject_holdingwater <- emmeans(red_model3_water_pres, specs = pairwise ~ interTypecbType, type = "response", adjust = "tukey")
contrast(emobject_holdingwater, "revpairwise", type = "response", adjust = "tukey")


## Model diagnostics

# Calculate residuals
simulationOutput_waterPresAll <- simulateResiduals(fittedModel = red_model3_water_pres, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_waterPresAll)
plotQQunif(simulationOutput_waterPresAll) # Left plot
plotResiduals(simulationOutput_waterPresAll) # Right plot


#### Visualize standing water presence (all catch basins) ----

standingwater_all$Dry_factor <- factor(standingwater_all$Dry)

# Sample sizes and summary stats
standingwater_all %>%
  group_by(interTypecbType) %>%
  summarise(
    total_obs = n(),
    total_notholdingwater = sum(Dry == 1),
    percent_notholdingwater = (total_notholdingwater/total_obs) *100
  )


standing_water_interTypecbType <- ggplot(standingwater_all, aes(x = interTypecbType, fill = Dry_factor)) +
  geom_bar(position = "fill", color = "black") +
  labs(x = "Catch Basin Type", y = "Relative Frequency") +
  scale_fill_manual(values = c("dodgerblue1","wheat3"), labels = c("Holding\nWater", "Not Holding\nWater")) +
  scale_x_discrete(labels = c("Control\nConventional", 
                              "Rain Garden\nConventional",
                              "Rain Garden\nInfiltration",
                              "Rain Garden\nOverflow")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = seq(0,1,0.25)) +
  annotate("text", x = 1, y = 1.02, label = expression(italic(n) == 421), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = 1.02, label = expression(italic(n) == 128), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 3, y = 1.02, label = expression(italic(n) == 985), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 4, y = 1.02, label = expression(italic(n) == 221), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = 1.07, label = "Control Intersection", size = 2, family = "HelveticaNeueforSAS") +
  annotate("text", x = 3, y = 1.07, label = "Rain Garden Intersection", size = 2, family = "HelveticaNeueforSAS") +
  geom_bracket(data = standingwater_all, xmin = 1.55, xmax = 4.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # RG Intersection grouping
  geom_bracket(data = standingwater_all, xmin = 0.55, xmax = 1.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # Control Intersection grouping
  figtheme + theme(legend.title = element_blank(), legend.key.size = unit(4, "mm"))


standing_water_interTypecbType

#### Data organizing for analyses involving juvenile mosquitoes ----

wet_CBs <- rg_clean %>% 
  subset(Dry == 0 & Sampled == 1) %>% # Select catch basins that were holding water and were sampled for juvenile mosquitoes
  mutate(interTypecbType = as.character(interTypecbType)) %>% # Change data type for subsetting
  subset(interTypecbType != "Rain Garden-RGI") %>% # Drop single obs where a RGI was holding water and sampled for mosquitoes
  mutate(interTypecbType = factor(interTypecbType)) # Convert back to factor

# Verify that samples that weren't processed were field neg and should have 0s in the juvenile mosquito-related columns

# Subset wet CBs where Processed == 0 (samples WEREN'T processed)
processed0 <- subset(wet_CBs, Processed == 0)

# Subset CBs in processed0 where all juvenile mosquito-related cols are NULL
juv_mos_cols <- c("L1L3_adj", "L4_Cx_total", "P_Cx_Cx_adj", "Other_adj")
# Checks each row to see if all values are NA in the juv_mos_cols
processed0_all_null <- processed0[apply(processed0[juv_mos_cols], 1, function(x) all(is.na(x))), ]

# Verify that the lengths of these two dataframes are the same
# This means that all wet CBs with NA for all juvenile mosquito-related columns are NEGATIVE for juvenile mosquitoes
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

## Selecting best temperature variable
tempVars <- rg_raw %>% # Create dataframe with selected temp cols
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
  model <- glm(formula, data = mosquitoes_all, family = binomial(link = "logit"))
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
prcpVars <- rg_raw %>%
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
  model <- glm(formula, data = mosquitoes_all, family = binomial(link = "logit"))
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
# Response variable: juvenile_pos (1 = negative for juveniles; 0 = positive for juveniles)

# Predictor variables:
# NSA_avgT_d0_d4
# NS_PRCP_mm_sum_d7_d13    
# interTypecbType

# Random effects
# (Inter_ID/CB_ID_RG)
# EpiweekYear

# Check data types
str(mosquitoes_all$juvenile_pos)
str(mosquitoes_all$NSA_avgT_d0_d4)
str(mosquitoes_all$NS_PRCP_mm_sum_d7_d13)
str(mosquitoes_all$interTypecbType)
str(mosquitoes_all$Inter_ID)
str(mosquitoes_all$CB_ID_RG)
str(mosquitoes_all$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models

# Full model
full_model_mos_pres <- glmmTMB(juvenile_pos ~ NSA_avgT_d0_d4 + NS_PRCP_mm_sum_d7_d13 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                               data = mosquitoes_all,
                               family = binomial(link= "logit"))

# Reduced model #1: dropped temp variable
red_model1_mos_pres <- glmmTMB(juvenile_pos ~ NS_PRCP_mm_sum_d7_d13 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                               data = mosquitoes_all,
                               family = binomial(link= "logit"))

# Reduced model #2: dropped prcp variable
red_model2_mos_pres <- glmmTMB(juvenile_pos ~ NSA_avgT_d0_d4 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                               data = mosquitoes_all,
                               family = binomial(link= "logit"))

# Reduced model #3: dropped temp and prcp variables
red_model3_mos_pres <- glmmTMB(juvenile_pos ~ interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                               data = mosquitoes_all,
                               family = binomial(link= "logit"))

# Model selection
model.sel(full_model_mos_pres, red_model1_mos_pres, red_model2_mos_pres, red_model3_mos_pres,
          rank = AICc)


# Selected model: reduced model #1

# Summary
summary_allCBs_holdingmosquitoes <- summary(red_model1_mos_pres)
summary_allCBs_holdingmosquitoes

exp(summary_allCBs_holdingmosquitoes$coefficients$cond)
exp(confint(red_model1_mos_pres))

# Anova
Anova(red_model1_mos_pres, type = 3)

# Emmeans
emobject_holdingmosquitoes <- emmeans(red_model1_mos_pres, specs = pairwise ~ interTypecbType, adjust = "tukey", type = "response")
contrast(emobject_holdingmosquitoes, "revpairwise", type = "response", adjust = "tukey")

# Model diagnostics

# Calculate residuals
simulationOutput_mosquitoPresAll <- simulateResiduals(fittedModel = red_model1_mos_pres, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_mosquitoPresAll)
plotQQunif(simulationOutput_mosquitoPresAll) # Left plot
plotResiduals(simulationOutput_mosquitoPresAll) # Right plot

#### Visualize juvenile mosquito presence ----
mosquitoes_all$juvenile_pos_factor <- factor(mosquitoes_all$juvenile_pos)

# Samples sizes and summary stats
mosquitoes_all %>%
  group_by(interTypecbType) %>%
  summarise(
    total_obs = n(),
    total_notholdingmosquitoes = sum(juvenile_pos == 1),
    percent_notholdingmosquitoes = (total_notholdingmosquitoes/total_obs) *100
  )


juvenile_pos_interTypecbType <- ggplot(mosquitoes_all, aes(x = interTypecbType, fill = juvenile_pos_factor)) +
  geom_bar(position = "fill", color = "black") +
  labs(x = "Catch Basin Type", y = "Relative Frequency") +
  scale_fill_manual(values = c("firebrick4","wheat3"), labels = c("Holding\nJuvenile Mosquitoes", "Not Holding\nJuvenile Mosquitoes")) +
  scale_x_discrete(labels = c("Control\nConventional", 
                              "Rain Garden\nConventional",
                              "Rain Garden Overflow")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1), breaks = seq(0,1,0.25)) +
  annotate("text", x = 1, y = 1.02, label = expression(italic(n) == 109), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = 1.02, label = expression(italic(n) == 22), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 3, y = 1.02, label = expression(italic(n) == 83), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 1, y = 1.07, label = "Control Intersection", size = 2, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2.5, y = 1.07, label = "Rain Garden Intersection", size = 2, family = "HelveticaNeueforSAS") +
  geom_bracket(data = mosquitoes_all, xmin = 1.55, xmax = 3.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # RG Intersection grouping
  geom_bracket(data = mosquitoes_all, xmin = 0.55, xmax = 1.45, y.position = 1.04, label = "", tip.length = c(0.5, 0.5)) + # Control Intersection grouping
  figtheme + theme(legend.title = element_blank(), legend.key.size = unit(4, "mm"))

juvenile_pos_interTypecbType


#### Combined figure (standing water and juvenile mosquito presence) ----

Fig5 <- standing_water_interTypecbType + juvenile_pos_interTypecbType +
  plot_layout(axis_titles = "collect", axes = "collect_y") +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag.position = c(0.02, 0.97), plot.tag = element_text(size = 9, face = "bold", family = "HelveticaNeueforSAS"))

Fig5

ggsave(
  plot = Fig5,
  filename = "figures/manuscript/Figure5.jpeg",
  device ="jpeg",
  units = "mm",
  height = 120, width = 180, dpi = 300
)

#### Data organizing for abundance models ----

# Data for abundance models 
abundance <- wet_CBs %>%
  mutate(L1L3_recode = ifelse(is.na(L1L3_adj), 0, L1L3_adj)) %>% # Convert NAs to 0s
  mutate(L4_recode = ifelse(is.na(L4_Cx_total ), 0, L4_Cx_total)) %>% # Convert NAs to 0s
  mutate(P_recode = ifelse(is.na(P_Cx_Cx_adj), 0, P_Cx_Cx_adj))  # Convert NAs to 0s

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
sum(is.na(L1L3$Salinity_ppm)) # 2 NA obs

# Relative proportion of water covered by floating OM

# Possible values
unique(L1L3$Surf_Detrit)
sum(is.na(L1L3$Surf_Detrit)) # 2 missing obs

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
  filter(!is.na(Salinity_ppm))

# Response: L1L3_recode

# Predictor variables:
# NS_PRCP_mm_sum_d0_d1 
# NSA_TMIN_C 
# Salinity_ppm
# interTypecbType

# Random effects:
# Inter_ID/CB_ID_RG
# EpiweekYear

# Check data types
str(L1L3$L1L3_recode)
str(L1L3$NS_PRCP_mm_sum_d0_d1)
str(L1L3$NSA_TMIN_C)
str(L1L3$Salinity_ppm)
str(L1L3$interTypecbType)
str(L1L3$Inter_ID)
str(L1L3$CB_ID_RG)
str(L1L3$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly"))

## Models

# Full model
full_model_L1_L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d1 + NSA_TMIN_C + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                            data = L1L3,
                            family = nbinom2(link = "log"))

# Reduced model 1: dropped temp var
red_model1_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d1 + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 2: dropped prcp var
red_model2_L1L3 <- glmmTMB(L1L3_recode ~ NSA_TMIN_C + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 3: dropped salinity var
red_model3_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d1 + NSA_TMIN_C + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 4: dropped temp and prcp vars
red_model4_L1L3 <- glmmTMB(L1L3_recode ~ Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 5: dropped temp and salinity vars
red_model5_L1L3 <- glmmTMB(L1L3_recode ~ NS_PRCP_mm_sum_d0_d1 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 6: dropped prcp and salinity vars
red_model6_L1L3 <- glmmTMB(L1L3_recode ~  NSA_TMIN_C + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Reduced model 7: dropped temp, prcp, and salinity vars
red_model7_L1L3 <- glmmTMB(L1L3_recode ~ interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                           data = L1L3,
                           family = nbinom2(link = "log"))

# Compare models
model.sel(full_model_L1_L3, red_model1_L1L3, red_model2_L1L3, red_model3_L1L3, red_model4_L1L3, red_model5_L1L3, red_model6_L1L3, red_model7_L1L3,
          rank = AICc)

# Selected model: reduced model 5

# Summary
summary(red_model5_L1L3)

confint(red_model5_L1L3)

# Anova
Anova(red_model5_L1L3, type = 3)

# Emmeans
emobject_L1L3 <- emmeans(red_model5_L1L3, specs = pairwise ~ interTypecbType, adjust = "tukey", type = "response")
contrast(emobject_L1L3, "revpairwise", type = "response", adjust = "tukey")

## Model diagnostics

# Calculate residuals
simulationOutput_L1L3 <- simulateResiduals(fittedModel = red_model5_L1L3, plot = FALSE)

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
sum(is.na(L4$Salinity_ppm)) # 2 NA obs

# Relative proportion of water covered by floating OM

# Possible values
unique(L4$Surf_Detrit)
sum(is.na(L4$Surf_Detrit)) # 2 missing obs

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

## negative binomial glmm

# Filter out any missing data for predictors
L4 <- L4 %>%
  filter(!is.na(Salinity_ppm))

# Response: L4_recode

# Predictor variables:
# NS_PRCP_mm_sum_d7_d13  
# NSA_avgT_d0_d2 
# Salinity_ppm
# interTypecbType

# Random effects:
# Inter_ID/CB_ID_RG
# EpiweekYear

# Check data types
str(L4$L4_recode)
str(L4$NS_PRCP_mm_sum_d7_d13)
str(L4$NSA_avgT_d0_d2)
str(L4$Salinity_ppm)
str(L4$interTypecbType)
str(L4$CB_ID_RG)
str(L4$Inter_ID)
str(L4$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) # Compare to ref level (control intersection conventional catch basins)

## Models

# Full model
full_model_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d7_d13 + NSA_avgT_d0_d2 + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 1: dropped temp var
red_model1_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d7_d13 + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 2: dropped prcp var
red_model2_L4 <- glmmTMB(L4_recode ~ NSA_avgT_d0_d2 + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 3: dropped salinity var
red_model3_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d7_d13 + NSA_avgT_d0_d2 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 4: dropped temp and prcp vars
red_model4_L4 <- glmmTMB(L4_recode ~ Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 5: dropped temp and salinity vars
red_model5_L4 <- glmmTMB(L4_recode ~ NS_PRCP_mm_sum_d7_d13 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 6: dropped prcp and salinity vars
red_model6_L4 <- glmmTMB(L4_recode ~  NSA_avgT_d0_d2 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Reduced model 7: dropped temp, prcp, and salinity vars
red_model7_L4 <- glmmTMB(L4_recode ~ interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                         data = L4,
                         family = nbinom2(link = "log"))

# Compare models
model.sel(full_model_L4, red_model1_L4, red_model2_L4, red_model3_L4, red_model4_L4, red_model5_L4, red_model6_L4, red_model7_L4,
          rank = AICc)

# Selected model: reduced model 3

# Summary
summary(red_model3_L4)

confint(red_model3_L4)

# Anova
Anova(red_model3_L4, type = 3)

# Emmeans
emobject_L14 <- emmeans(red_model3_L4, specs = pairwise ~ interTypecbType, adjust = "tukey", type = "response")
contrast(emobject_L14, "revpairwise", type = "response", adjust = "tukey")

## Model diagnostics

# Calculate residuals
simulationOutput_L4 <- simulateResiduals(fittedModel = red_model3_L4, plot = FALSE)

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
sum(is.na(pupae$Salinity_ppm)) # 2 NA obs

# Relative proportion of water covered by floating OM

# Possible values
unique(pupae$Surf_Detrit)
sum(is.na(pupae$Surf_Detrit)) # 2 missing obs

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

## Negative binomial glmm

# Filter out any missing data for predictors
pupae <- pupae %>%
  filter(!is.na(Salinity_ppm))

# Response: P_recode

# Predictor variables:
# NS_PRCP_mm_sum_d7_d13  
# NSA_TMIN_C 
# Salinity_ppm
# interTypecbType

# Random effects:
# Inter_ID/CB_ID_RG
# EpiweekYear

# Check data types
str(pupae$P_recode)
str(pupae$NS_PRCP_mm_sum_d7_d13)
str(pupae$NSA_TMIN_C)
str(pupae$Salinity_ppm)
str(pupae$interTypecbType)
str(pupae$Inter_ID)
str(pupae$CB_ID_RG)
str(pupae$EpiweekYear)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models

# Full model
full_model_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d7_d13 + NSA_TMIN_C + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 1: dropped temp var
red_model1_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d7_d13 + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 2: dropped prcp var
red_model2_P <- glmmTMB(P_recode ~ NSA_TMIN_C + Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 3: dropped salinity var
red_model3_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d7_d13 + NSA_TMIN_C + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 4: dropped temp and prcp vars
red_model4_P <- glmmTMB(P_recode ~ Salinity_ppm + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 5: dropped temp and salinity vars
red_model5_P <- glmmTMB(P_recode ~ NS_PRCP_mm_sum_d7_d13 + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 6: dropped prcp and salinity vars
red_model6_P <- glmmTMB(P_recode ~  NSA_TMIN_C + interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Reduced model 7: dropped temp, prcp, and salinity vars
red_model7_P <- glmmTMB(P_recode ~ interTypecbType + (1|Inter_ID/CB_ID_RG) + (1|EpiweekYear),
                        data = pupae,
                        family = nbinom2(link = "log"))

# Compare models
model.sel(full_model_P, red_model1_P, red_model2_P, red_model3_P, red_model4_P, red_model5_P, red_model6_P, red_model7_P,
          rank = AICc)

# Selected model: reduced model 3

# Summary
summary(red_model3_P)

confint(red_model3_P, level = 0.95)

# Anova
Anova(red_model3_P, type = 3)

# Emmeans
emobject_P <- emmeans(red_model3_P, specs = pairwise ~ interTypecbType, adjust = "tukey", type = "response")
contrast(emobject_P, "revpairwise", type = "response", adjust = "tukey")


# Model diagnostics

# Calculate residuals
simulationOutput_P <- simulateResiduals(fittedModel = red_model3_P, plot = FALSE)

# Test zero inflation
testZeroInflation(simulationOutput_P)

# Plotting the scaled residuals
plot(simulationOutput_P)
plotQQunif(simulationOutput_P) # Left plot
plotResiduals(simulationOutput_P) # Right plot

#### Abundance figures ----

L1L3 %>%
  group_by(interTypecbType) %>%
  summarise(
    total_obs = n(),
    avg_L1L3 = mean(L1L3_recode),
    se_L1L3 = sd(L1L3_recode) / sqrt(total_obs)
  )

L4 %>%
  group_by(interTypecbType) %>%
  summarise(
    total_obs = n(),
    avg_L4 = mean(L4_recode),
    se_L4 = sd(L4_recode) / sqrt(total_obs)
  )

pupae %>%
  group_by(interTypecbType) %>%
  summarise(
    total_obs = n(),
    avg_P = mean(P_recode),
    se_P = sd(P_recode) / sqrt(total_obs)
  )


# Catch basin labels 
CB_labs = c("Control\nConventional", "Rain Garden\nConventional", "Rain Garden\nOverflow")

# L1L3
L1L3_rg <- ggplot(L1L3, aes(x = interTypecbType, y = L1L3_recode, fill = interTypecbType, color = interTypecbType)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +
  labs(x = "Catch Basin Type", y = "Average 1st to 3rd Instar\nLarval Abundance") +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("grey27", "green4", "mediumblue")) +
  scale_x_discrete(labels = CB_labs) +
  scale_y_continuous(breaks = seq(0,2500,500)) +
  annotate("text", x = 1, y = -115, label = expression(italic(n) == 108), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -115, label = expression(italic(n) == 21), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 3, y = -115, label = expression(italic(n) == 83), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme + theme(legend.position = "none")

L1L3_rg


# L4
L4_rg <- ggplot(L4, aes(x = interTypecbType, y = L4_recode, fill = interTypecbType, color = interTypecbType)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +
  labs(x = "Catch Basin Type", y = "Average 4th Instar\nLarval Abundance") +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("grey27", "green4", "mediumblue")) +
  scale_x_discrete(labels = CB_labs) +
  annotate("text", x = 1, y = -15, label = expression(italic(n) == 108), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -15, label = expression(italic(n) == 21), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 3, y = -15, label = expression(italic(n) == 83), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme + theme(legend.position = "none")

L4_rg

# Pupae
pupae_rg <- ggplot(pupae, aes(x = interTypecbType, y = P_recode, fill = interTypecbType, color = interTypecbType)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +
  labs(x = "Catch Basin Type", y = "Average Pupal Abundance") +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("grey27", "green4", "mediumblue")) +
  scale_x_discrete(labels = CB_labs) +
  annotate("text", x = 1, y = -3, label = expression(italic(n) == 108), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 2, y = -3, label = expression(italic(n) == 21), size = N_size_size, family = "HelveticaNeueforSAS") +
  annotate("text", x = 3, y = -3, label = expression(italic(n) == 83), size = N_size_size, family = "HelveticaNeueforSAS") +
  figtheme + theme(legend.position = "none")

pupae_rg

abundance_combo <- L1L3_rg + L4_rg + pupae_rg +
  plot_layout(guides = "collect", axis_titles = "collect", nrow = 3) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag.position = c(0.02, 0.97), plot.tag = element_text(size = 9, face = "bold", family = "HelveticaNeueforSAS"))
abundance_combo


# Save to Aurora Box folder
ggsave(
  plot = abundance_combo,
  filename = "figures/manuscript/Exp2_abundance_figs.jpeg",
  device ="jpeg",
  units = "mm",
  height = 180, width = 88, dpi = 300
)

