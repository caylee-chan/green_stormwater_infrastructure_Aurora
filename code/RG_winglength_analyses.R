# Experiment #2: Wing length analysis
# Caylee Chan
# Created: 14 Jan 2025
# Updated:23 Jan 2026
# Notes: 


# Libraries:
library(dplyr)
library(ggplot2)
library(tidyr)
library(car)
library(emmeans)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(extrafont)

# Set wd
setwd("C:/Users/cayle.LAPTOP-QMLQMN4A/OneDrive - University of Illinois - Urbana/green_stormwater_infrastructure_Aurora/green_stormwater_infrastructure_Aurora")


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

# Set legend labels and colors
legend_labs <- c("Conventional\nCatch Basin", "Rain Garden\nOverflow\nCatch Basin")
legend_colors <- c("grey27", "green4")

# Set N size size
N_size_size <- 2


#### Data organizing ----

# Raw winglength data
raw_wing <- read.csv("data/Rain gardens/pupal sampling/starvation_resistance/Aurora_2015_starvation_resistance_20241007.csv")

# Subset Cx. pipiens and remove NA wing lengths
wing_clean_all <- raw_wing %>%
  filter(!is.na(Winglength_mm)) %>%
  subset(Species == "Cx. pipiens") %>%
  mutate(CB_Class = factor(CB_Class)) %>% # Change data type
  mutate(CBID = factor(CBID)) %>% # Change data type
  mutate(Epiweek = factor(Epiweek)) # Change data type

# Raw weather data
weather_raw <- read.csv("data/weather/Weather_daily_lags_Apr2013toOct2015.csv")

weatherClean <- weather_raw %>%
  select(DATE:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NSA_minT_d0_d1)

# Merge wing length and weather data
wingWeather_all <- merge(x = wing_clean_all,
                         y = weatherClean,
                         by.x = "DateColl",
                         by.y = "DATE",
                         all.x = FALSE)

## Subset females and males
females_wing <- subset(wingWeather_all, Sex == "F")
males_wing <- subset(wingWeather_all, Sex == "M")

#### Summary statistics and sample sizes and figure ----

females_wing %>%
  group_by(CB_Class) %>%
  summarise(
    total_obs = n(),
    mean_wing = mean(Winglength_mm),
    se_wing = sd(Winglength_mm) / sqrt(total_obs)
  )

males_wing %>%
  group_by(CB_Class) %>%
  summarise(
    total_obs = n(),
    mean_wing = mean(Winglength_mm),
    se_wing = sd(Winglength_mm) / sqrt(total_obs)
  )

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

#### Wing length model (FEMALES) ----

## Selecting best temperature variable - females
wing_tempVarsF <- females_wing %>% # Create dataframe with selected temp cols and Winglength_mm
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1, Winglength_mm)

# Correlation btwn tempVars and Winglength_mm
corrs_wing_tempVarsF <- cor(x = wing_tempVarsF[, colnames(wing_tempVarsF) != "Winglength_mm"], 
                            y = wing_tempVarsF$Winglength_mm, 
                            method = "spearman")

corrs_wing_tempVarsF <- as.data.frame(corrs_wing_tempVarsF)
colnames(corrs_wing_tempVarsF)[1] <- "r"

# Make new column with absolute value of r values
corrs_wing_tempVarsF$absval_r <- abs(corrs_wing_tempVarsF$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_wing_tempVarsF$absval_r, na.rm = TRUE)

# Extract the temp variable corresponding to the greatest correlation 
corrs_wing_tempVarsF %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable - females
wing_prcpVarsF <- females_wing %>% # Create dataframe with selected prcp cols and Winglength_mm
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20, Winglength_mm)

# Correlation btwn prcpVars and Winglength_mm
corrs_wing_prcpVarsF <- cor(x = wing_prcpVarsF[, colnames(wing_prcpVarsF) != "Winglength_mm"], 
                            y = wing_prcpVarsF$Winglength_mm,
                            method = "spearman")

corrs_wing_prcpVarsF <- as.data.frame(corrs_wing_prcpVarsF)
colnames(corrs_wing_prcpVarsF)[1] <- "r"

# Make new column with absolute value of r values
corrs_wing_prcpVarsF$absval_r <- abs(corrs_wing_prcpVarsF$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_wing_prcpVarsF$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_wing_prcpVarsF %>%
  filter(absval_r == maxAbsval_r)

## GLMM - females
# Response variable: Winglength_mm

# Predictor variables:
# NSA_avgT_d0_d1  
# NS_PRCP_mm_sum_Lag9    
# CB_Class

# Random effects:
# CB_ID_G
# Epiweek

# Verify data types
str(females_wing$Winglength_mm)
str(females_wing$NSA_avgT_d0_d1)
str(females_wing$NS_PRCP_mm_sum_Lag9)
str(females_wing$CB_Class)
str(females_wing$CBID)
str(females_wing$Epiweek)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models - females

# Full model - females
full_model_F <- glmmTMB(Winglength_mm ~ NSA_avgT_d0_d1 + NS_PRCP_mm_sum_Lag9 + CB_Class + (1|CBID) + (1|Epiweek),
                        data = females_wing,
                        family = gaussian(link = "identity"))

# Reduced model #1: dropped temperature variable
red_model1_F <- glmmTMB(Winglength_mm ~ NS_PRCP_mm_sum_Lag9 + CB_Class + (1|CBID) + (1|Epiweek),
                        data = females_wing,
                        family = gaussian(link = "identity"))

# Reduced model #2: dropped precipitation variable
red_model2_F <- glmmTMB(Winglength_mm ~ NSA_avgT_d0_d1 + CB_Class + (1|CBID) + (1|Epiweek),
                        data = females_wing,
                        family = gaussian(link = "identity"))

# Reduced model #3: dropped temperature and precipitation variables
red_model3_F <- glmmTMB(Winglength_mm ~ CB_Class + (1|CBID) + (1|Epiweek),
                        data = females_wing,
                        family = gaussian(link = "identity"))

# Compare models
model.sel(full_model_F, red_model1_F, red_model2_F, red_model3_F,
          rank = AICc)


# Selected model: full model

# Summary
summary(full_model_F)

confint(full_model_F, level = 0.95)

# Anova
Anova(full_model_F, type = 3)

# Emmeans
emobject_females <- emmeans(full_model_F, specs = pairwise ~ CB_Class, type = "response", adjust = "tukey")
contrast(emobject_females, "revpairwise", type = "response", adjust = "tukey")

## Model diagnostics - females

# Calculate residuals
simulationOutput_wingFemales <- simulateResiduals(fittedModel = full_model_F, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_wingFemales)
plotQQunif(simulationOutput_wingFemales) # Left plot
plotResiduals(simulationOutput_wingFemales) # Right plot

#### Wing length model (MALES) ----

## Selecting best temperature variable- males
wing_tempVarsM <- males_wing %>% # Create dataframe with selected temp cols and Winglength_mm
  select(NSA_TMAX_C:NSA_TAVG_C, NSA_avgT_d0_d1:NSA_minT_d0_d1, Winglength_mm)

# Correlation btwn tempVars and Winglength_mm
corrs_wing_tempVarsM <- cor(x = wing_tempVarsM[, colnames(wing_tempVarsM) != "Winglength_mm"], 
                            y = wing_tempVarsM$Winglength_mm, 
                            method = "spearman")

corrs_wing_tempVarsM <- as.data.frame(corrs_wing_tempVarsM)
colnames(corrs_wing_tempVarsM)[1] <- "r"

# Make new column with absolute value of r values
corrs_wing_tempVarsM$absval_r <- abs(corrs_wing_tempVarsM$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_wing_tempVarsM$absval_r, na.rm = TRUE)

# Extract the temp variable corresponding to the greatest correlation 
corrs_wing_tempVarsM %>%
  filter(absval_r == maxAbsval_r)

## Selecting best precipitation variable - males
wing_prcpVarsM <- males_wing %>% # Create dataframe with selected prcp cols and Winglength_mm
  select(NS_AVG_PRCP_mm:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NS_PRCP_mm_sum_d14_d20, Winglength_mm)

# Correlation btwn prcpVars and Winglength_mm
corrs_wing_prcpVarsM <- cor(x = wing_prcpVarsM[, colnames(wing_prcpVarsM) != "Winglength_mm"], 
                            y = wing_prcpVarsM$Winglength_mm,
                            method = "spearman")

corrs_wing_prcpVarsM <- as.data.frame(corrs_wing_prcpVarsM)
colnames(corrs_wing_prcpVarsM)[1] <- "r"

# Make new column with absolute value of r values
corrs_wing_prcpVarsM$absval_r <- abs(corrs_wing_prcpVarsM$r)

# Find the max absolute value for all r values
maxAbsval_r <- max(corrs_wing_prcpVarsM$absval_r, na.rm = TRUE)

# Extract the prcp variable corresponding to the greatest correlation 
corrs_wing_prcpVarsM %>%
  filter(absval_r == maxAbsval_r)

## GLMM - males
# Response variable: Winglength_mm

# Predictor variables:
# NSA_minT_d0_d1    
# NS_PRCP_mm_sum_d5_d7     
# CB_Class

# Random effects:
# CB_ID_G
# Epiweek

# Verify data types
str(males_wing$Winglength_mm)
str(males_wing$NSA_minT_d0_d1)
str(males_wing$NS_PRCP_mm_sum_d5_d7)
str(males_wing$CB_Class)
str(males_wing$CBID)
str(males_wing$Epiweek)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models - males

# Full model - males
full_model_M <- glmmTMB(Winglength_mm ~ NSA_minT_d0_d1 + NS_PRCP_mm_sum_d5_d7 + CB_Class + (1|CBID) + (1|Epiweek),
                        data = males_wing,
                        family = gaussian(link = "identity"))

# Reduced model #1: dropped temperature variable
red_model1_M <- glmmTMB(Winglength_mm ~ NS_PRCP_mm_sum_d5_d7 + CB_Class + (1|CBID) + (1|Epiweek),
                        data = males_wing,
                        family = gaussian(link = "identity"))

# Reduced model #2: dropped precipitation variable
red_model2_M <- glmmTMB(Winglength_mm ~ NSA_minT_d0_d1 + CB_Class + (1|CBID) + (1|Epiweek),
                        data = males_wing,
                        family = gaussian(link = "identity"))

# Reduced model #3: dropped temperature and precipitation variables
red_model3_M <- glmmTMB(Winglength_mm ~ CB_Class + (1|CBID) + (1|Epiweek),
                        data = males_wing,
                        family = gaussian(link = "identity"))

# Compare models
model.sel(full_model_M, red_model1_M, red_model2_M, red_model3_M, 
          rank = AICc)

# Selected model: reduced model 1

# Summary
summary(red_model1_M)

confint(red_model1_M, level = 0.95)

# Anova
Anova(red_model1_M, type = 3)

# Emmeans
emobject_males <- emmeans(red_model1_M, specs = pairwise ~ CB_Class, type = "response", adjust = "tukey")
contrast(emobject_males, "revpairwise", type = "response", adjust = "tukey")

# Model diagnostics - males

# Calculate residuals
simulationOutput_wingMales <- simulateResiduals(fittedModel = red_model1_M, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_wingMales)
plotQQunif(simulationOutput_wingMales) # Left plot
plotResiduals(simulationOutput_wingMales) # Right plot

