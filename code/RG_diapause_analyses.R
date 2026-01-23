# Experiment #2: Diapause analysis
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
library(formula.tools)
library(glmmTMB)
library(MuMIn)
library(extrafont)
library(DHARMa)

# Set wd
setwd("C:/Users/cayle.LAPTOP-QMLQMN4A/OneDrive - University of Illinois - Urbana/green_stormwater_infrastructure_Aurora/green_stormwater_infrastructure_Aurora")


#### Determining date windows for water temp ----

# Read in raw diapause data
raw_diapause <- read.csv("data/Rain gardens/pupal sampling/diapause/Aurora_2015_diapause_data_20241007.csv")

unique(raw_diapause$Coll_Date)
unique(raw_diapause$Species)

# Experiment #1 - collection date = 2 Sept 2015
exp1 <- subset(raw_diapause, Coll_Date == "2-Sep-15")
unique(exp1$F_emerge_last)

# Experiment #2 (1) - collection date = 8 Sept 2015
exp2_1 <- subset(raw_diapause, Coll_Date == "8-Sep-15")
unique(exp2_1$F_emerge_last)

# Experiment #2 (2) - collection date = 9 Sept 2015
exp2_2 <- subset(raw_diapause, Coll_Date == "9-Sep-15")
unique(exp2_2$F_emerge_last)

# Experiment #3 - collection date = 16 Sept 2015
exp3 <- subset(raw_diapause, Coll_Date == "16-Sep-15")
unique(exp3$F_emerge_last)

#### Calculating water temperature -----------------------------------------------------------------------------------------------------------------------

# Read in raw temp data
temp <- read.csv("data/Rain gardens/pupal sampling/HOBO_2015_CCB_RGO_1hr_int.csv")

# Create datetime col
temp$Date <- as.Date(temp$Date, "%d-%b-%y")
temp$Time_mean <- strptime(temp$Time_mean, "%I:%M:%S %p")
temp$Time <- format(temp$Time_mean, "%H:%M:%S")

temp$datetime <- as.POSIXct(paste(temp$Date, temp$Time), format = "%Y-%m-%d %H:%M:%S")

## Experiment #1 - collection date: 2 Sept 2015

# Temp 1: 31 Aug - 2 Sept @12PM

# Select dates
exp1_temp1Dates <- subset(temp, datetime >= "2015-08-31 00:00:00" & datetime <= "2015-09-02 12:00:00") 

# RGO temp1
exp1_temp1_RGO <- mean(exp1_temp1Dates$RGO_avg)
exp1_temp1_RGO

# CCB temp1
exp1_temp1_CCB <- mean(exp1_temp1Dates$CCB_avg)
exp1_temp1_CCB

# Temp 2: 2 Sept @12PM - 5 Sept

# Select dates
exp1_temp2Dates <- subset(temp, datetime >= "2015-09-02 12:00:00" & datetime <= "2015-09-05 23:59:59")

# Experiment 1 temp2
exp1_temp2 <- mean(exp1_temp2Dates$D1_Temp_C)
exp1_temp2

# tempFinal Experiment #1
exp1RGO_tempFinal <- mean(c(exp1_temp1_RGO, exp1_temp2))
exp1RGO_tempFinal

exp1CCB_tempFinal <- mean(c(exp1_temp1_CCB, exp1_temp2))
exp1CCB_tempFinal


## Experiment 2a - collection date: 8 Sept 2015

# Temp 1: 6 Sept - 8 Sept @12PM

# Select dates
exp2a_temp1Dates <- subset(temp, datetime >= "2015-09-06 00:00:00" & datetime <= "2015-09-08 12:00:00") 

# RGO temp1
exp2a_temp1_RGO <- mean(exp2a_temp1Dates$RGO_avg)
exp2a_temp1_RGO

# CCB temp1
exp2a_temp1_CCB <- mean(exp2a_temp1Dates$CCB_avg)
exp2a_temp1_CCB

# Temp 2: 8 Sept @12PM - 12 Sept

# Select dates
exp2a_temp2Dates <- subset(temp, datetime >= "2015-09-08 12:00:00" & datetime <= "2015-09-12 23:59:59")

# Experiment 2a temp2
exp2a_temp2 <- mean(exp2a_temp2Dates$D2_Temp_C)
exp2a_temp2

# tempFinal Experiemnt 2a
exp2aRGO_tempFinal <- mean(c(exp2a_temp1_RGO, exp2a_temp2))
exp2aRGO_tempFinal

exp2aCCB_tempFinal <- mean(c(exp2a_temp1_CCB, exp2a_temp2))
exp2aCCB_tempFinal


## Experiment 2b - collection date: 9 Sept 2015

# Temp 1: 7 Sept - 9 Sept @12PM

# Select dates
exp2b_temp1Dates <- subset(temp, datetime >= "2015-09-07 00:00:00" & datetime <= "2015-09-09 12:00:00") 

# RGO temp1
exp2b_temp1_RGO <- mean(exp2b_temp1Dates$RGO_avg)
exp2b_temp1_RGO

# CCB temp1
exp2b_temp1_CCB <- mean(exp2b_temp1Dates$CCB_avg)
exp2b_temp1_CCB

# Temp 2: 9 Sept @12PM - 12 Sept

# Select dates
exp2b_temp2Dates <- subset(temp, datetime >= "2015-09-09 12:00:00" & datetime <= "2015-09-12 23:59:59")

# Experiment 2b temp2
exp2b_temp2 <- mean(exp2b_temp2Dates$D2_Temp_C)
exp2b_temp2

# tempFinal Experiment 2b
exp2bRGO_tempFinal <- mean(c(exp2b_temp1_RGO, exp2b_temp2))
exp2bRGO_tempFinal

exp2bCCB_tempFinal <- mean(c(exp2b_temp1_CCB, exp2b_temp2))
exp2bCCB_tempFinal



## Experiment 3 - collection date: 16 Sept 2015

# Temp 1: 14 Sept - 16 Sept @12PM

# Select dates
exp3_temp1Dates <- subset(temp, datetime >= "2015-09-14 00:00:00" & datetime <= "2015-09-16 12:00:00") 

# RGO temp1
exp3_temp1_RGO <- mean(exp3_temp1Dates$RGO_avg)
exp3_temp1_RGO

# CCB temp1
exp3_temp1_CCB <- mean(exp3_temp1Dates$CCB_avg, na.rm = TRUE)
exp3_temp1_CCB

# Temp 2: 16 Sept @12PM - 19 Sept

# Select dates - can't select using datetime because original time not recorded after Sept 16
exp3_temp2Dates <- subset(temp, (Date == "2015-09-16" & (Time >= 12)) | Date == "2015-09-17" | Date == "2015-09-18" | Date == "2015-09-19")

# Experiment 3 temp2
exp3_temp2 <- mean(exp3_temp2Dates$D3_Temp_C)
exp3_temp2

# tempFinal Experiment 3
exp3RGO_tempFinal <- mean(c(exp3_temp1_RGO, exp3_temp2))
exp3RGO_tempFinal

exp3CCB_tempFinal <- mean(c(exp3_temp1_CCB, exp3_temp2))
exp3CCB_tempFinal

#### Assigning water temperature data to diapause dataframe

unique(raw_diapause$Coll_Date)

raw_diapause$tempFinal <- ifelse(raw_diapause$Coll_Date == "2-Sep-15" & raw_diapause$CB_CLASS == "RGO", exp1RGO_tempFinal,
                                 ifelse(raw_diapause$Coll_Date == "2-Sep-15" & raw_diapause$CB_CLASS == "C", exp1CCB_tempFinal,
                                        ifelse(raw_diapause$Coll_Date == "8-Sep-15" & raw_diapause$CB_CLASS == "RGO", exp2aRGO_tempFinal,
                                               ifelse(raw_diapause$Coll_Date == "8-Sep-15" & raw_diapause$CB_CLASS == "C", exp2aCCB_tempFinal,
                                                      ifelse(raw_diapause$Coll_Date == "9-Sep-15" & raw_diapause$CB_CLASS == "RGO", exp2bRGO_tempFinal,
                                                             ifelse(raw_diapause$Coll_Date == "9-Sep-15" & raw_diapause$CB_CLASS == "C", exp2bCCB_tempFinal,
                                                                    ifelse(raw_diapause$Coll_Date == "16-Sep-15" & raw_diapause$CB_CLASS == "RGO", exp3RGO_tempFinal,
                                                                           ifelse(raw_diapause$Coll_Date == "16-Sep-15" & raw_diapause$CB_CLASS == "C", exp3CCB_tempFinal, NA))))))))
# Save as csv 
write.csv(raw_diapause,"data/Rain gardens/pupal sampling/diapause/Aurora_2015_diapause_data_20241203_w_watertemp.csv", row.names = FALSE)

#### Diapause logistic regression model ----

## Read in raw diapause data with calculated temperatures
diapause <- read.csv("data/Rain gardens/pupal sampling/diapause/Aurora_2015_diapause_data_20241203_w_watertemp.csv")

# Data organizing
diapause <- diapause %>%
  mutate(CB_CLASS = factor(CB_CLASS)) %>% # Change data type
  mutate(CB_ID = factor(CB_ID)) # Change data type

## Mixed effects logistic regression model
# Response variable: Diapause (1 = diapausing female; 0 = non-diapausing female)

# Predictor variables:
# tempFinal
# Coll_daylength
# CB_CLASS

# Random effects:
# CB_ID

# Check data types
str(diapause$Diapause)
str(diapause$tempFinal)
str(diapause$Coll_daylength)
str(diapause$CB_CLASS)
str(diapause$CB_ID)

# Set contrasts
options(contrasts = c("contr.sum","contr.poly")) 

## Models

# Full model
full_model <- glmmTMB(Diapause ~ tempFinal + Coll_daylength + CB_CLASS + (1|CB_ID),
                      data = diapause,
                      family = binomial(link= "logit"))

# Reduced model #1: dropped tempFinal
red_model1 <- glmmTMB(Diapause ~ Coll_daylength + CB_CLASS + (1|CB_ID),
                      data = diapause,
                      family = binomial(link= "logit"))

# Reduced model #2: dropped Coll_daylength
red_model2 <- glmmTMB(Diapause ~ tempFinal + CB_CLASS + (1|CB_ID),
                      data = diapause,
                      family = binomial(link= "logit"))

# Compare models
model.sel(full_model, red_model1, red_model2,
          rank = AICc)

# Selected model: full model

# Summary
summary(full_model)

summary_diapause <- summary(full_model)
exp(summary_diapause$coefficients$cond)

exp(confint(full_model, level = 0.95))

# Anova
Anova(full_model, type = 3)

# Emmeans
emobject_diapause <- emmeans(full_model, specs = pairwise ~ CB_CLASS, type = "response", adjust = "tukey")
contrast(emobject_diapause, "revpairwise", type = "response", adjust = "tukey")

## Model diagnostics

# Calculate residuals
simulationOutput_diapause <- simulateResiduals(fittedModel = full_model, plot = FALSE)

# Plotting the scaled residuals
plot(simulationOutput_diapause)
plotQQunif(simulationOutput_diapause) # Left plot
plotResiduals(simulationOutput_diapause) # Right plot

#### Plot ----
diapause %>%
  group_by(CB_CLASS) %>%
  summarise(
    total_obs = n(),
    total_indiapause = sum(Diapause == 1),
    percent_indiapause = (total_indiapause / total_obs) * 100
  )


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
  legend.title = element_blank())

# Set N size size
N_size_size <- 2

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

ggsave(
  plot = diapause_viz,
  filename = "figures/manuscript/Exp2_diapause.jpeg",
  device ="jpeg",
  units = "mm",
  height = 100, width = 88, dpi = 300
)



#### Summary stats -------------
diapause %>%
  filter(CB_CLASS == "C") %>%
  summarise(
    percent_indiapause = sum(Diapause == 1) / n(),
    total_indiapause = sum(Diapause == 1),
    total_obs = n()
  )

diapause %>%
  filter(CB_CLASS == "RGO") %>%
  summarise(
    percent_indiapause = sum(Diapause == 1) / n(),
    total_indiapause = sum(Diapause == 1),
    total_obs = n()
  )

