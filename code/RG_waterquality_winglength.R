# Experiment #2: Water quality & wing length analysis
# Caylee Chan
# Created: 8 July 2025
# Updated: 1 June 2026
# Notes:


# Libraries:
library(dplyr)
library(ggplot2)
library(tidyr)
library(car)
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

# Read in raw water chemistry data
raw_water_chem <- read.csv("data/Rain gardens/pupal sampling/water_chemistry/GICSOCP_IPCA_water_chem_assays_raw.csv")

# Remove RGIs, select TrueValues != 99999.00
waterqual_clean <- raw_water_chem %>%
  subset(CB_type != "RGO-I") %>%
  subset(TrueValue != 99999.00) %>%
  # Calculate the mean TrueValue by collection date and CB ID by analyte (i.e., results in singular value when multiple dilutions used)
  group_by(Coll_Date, Epiweek, CB_ID_G, Inter_type, CB_type, Analyte) %>%
  summarise(
    avg_conc = mean(TrueValue)
  )

# Pivot wider to create columns for the concentrations of each analyte
waterqual_clean_wide <- pivot_wider(waterqual_clean, names_from = Analyte, values_from = avg_conc)
waterqual_clean_wide$Epiweek <- factor(waterqual_clean_wide$Epiweek)

# Number of CCBs and RGOs sampled for water quality data
n_distinct(waterqual_clean_wide$CB_ID_G[waterqual_clean_wide$CB_type == "C"]) # 17 CCBs
n_distinct(waterqual_clean_wide$CB_ID_G[waterqual_clean_wide$CB_type == "RGO"]) # 8 RGOs

#### NH3 N plots and stats ----

# NH3 N concentrations by Epiweek and catch basin type
NH3_N_byepiweek <- ggplot(data = subset(waterqual_clean_wide, !is.na(NH3_N)), aes(x = Epiweek, y = NH3_N, color = CB_type, fill = CB_type)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Epiweek", y = "Ammonia Nitrogen (mg/L)") +
  scale_color_manual(values = c("black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  figtheme

NH3_N_byepiweek

# Plot average NH3_N concentration by Epiweek and catch basin type
NH3_plotting <- waterqual_clean_wide %>%
  filter(!is.na(NH3_N)) %>%
  group_by(Epiweek, CB_type) %>%
  summarise(
    mean_NH3_byepiweek = mean(NH3_N)
  )

NH3_N_avg_byepiweek <- ggplot(NH3_plotting, aes(x = CB_type, y = mean_NH3_byepiweek, fill = CB_type)) +
  labs(x = "Catch Basin Type", y = "Average NH3 N \nConcentration (mg/L)by Epiweek") +
  scale_fill_manual(values = c("grey27", "green4")) +
  geom_boxplot() +
  theme_classic()
NH3_N_avg_byepiweek

# Nonparametric t-test
NH3_t_test_data <- pivot_wider(NH3_plotting,
                               names_from = CB_type,
                               names_prefix = "NH3_3_",
                               values_from = mean_NH3_byepiweek)

# Checking normality - NH3 in CCBs
shapiro.test(NH3_t_test_data$NH3_3_C)
qqnorm(NH3_t_test_data$NH3_3_C)
qqline(NH3_t_test_data$NH3_3_C)

# Checking normality - NH3 in RGOs
shapiro.test(NH3_t_test_data$NH3_3_RGO)
qqnorm(NH3_t_test_data$NH3_3_RGO)
qqline(NH3_t_test_data$NH3_3_RGO)

# Paired Wilcoxon test
wilcox.test(x = NH3_t_test_data$NH3_3_C, y = NH3_t_test_data$NH3_3_RGO,
            alternative = "two.sided", paired = FALSE)

t.test(x = NH3_t_test_data$NH3_3_C, y = NH3_t_test_data$NH3_3_RGO,
       alternative = "two.sided", paired = FALSE)

# Averaged concentrations
NH3_plotting %>%
  group_by(CB_type) %>%
  summarise(
    mean_NH3_N = mean(mean_NH3_byepiweek),
    size_n = n(),
    se_NH3_N = sd(mean_NH3_byepiweek)/sqrt(size_n)
  )

#### PO4 ----

# PO4 concentrations by Epiweek and catch basin type
PO4_byepiweek <- ggplot(data = subset(waterqual_clean_wide, !is.na(PO4)), aes(x = Epiweek, y = PO4, color = CB_type, fill = CB_type)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Epiweek", y = "Reactive phosohprus (mg/L)") +
  scale_color_manual(values = c("black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  figtheme

PO4_byepiweek

# Plot average PO4 concentration by Epiweek and catch basin type
PO4_plotting <- waterqual_clean_wide %>%
  filter(!is.na(PO4)) %>%
  group_by(Epiweek, CB_type) %>%
  summarise(
    mean_PO4_byepiweek = mean(PO4)
  )

PO4_avg_byepiweek <- ggplot(PO4_plotting, aes(x = CB_type, y = mean_PO4_byepiweek, fill = CB_type)) +
  labs(x = "Catch Basin Type", y = "Average PO4\nConcentration (mg/L)by Epiweek") +
  scale_fill_manual(values = c("grey27", "green4")) +
  geom_boxplot() +
  theme_classic()
PO4_avg_byepiweek

# Nonparametric t-test
PO4_t_test_data <- pivot_wider(PO4_plotting,
                               names_from = CB_type,
                               names_prefix = "PO4_",
                               values_from = mean_PO4_byepiweek)

# Checking normality - NH3 in CCBs
shapiro.test(PO4_t_test_data$PO4_C)
qqnorm(PO4_t_test_data$PO4_C)
qqline(PO4_t_test_data$PO4_C)

# Checking normality - NH3 in RGOs
shapiro.test(PO4_t_test_data$PO4_RGO)
qqnorm(PO4_t_test_data$PO4_RGO)
qqline(PO4_t_test_data$PO4_RGO)

# Paired Wilcoxon test
wilcox.test(x = PO4_t_test_data$PO4_C, y = PO4_t_test_data$PO4_RGO,
            alternative = "two.sided", paired = FALSE)

t.test(x = PO4_t_test_data$PO4_C, y = PO4_t_test_data$PO4_RGO,
       alternative = "two.sided", paired = FALSE)

# Averaged concentrations
PO4_plotting %>%
  group_by(CB_type) %>%
  summarise(
    mean_PO4 = mean(mean_PO4_byepiweek),
    size_n = n(),
    se_PO4 = sd(mean_PO4_byepiweek)/sqrt(size_n)
  )


#### TSS ----

# TSS concentrations by Epiweek and catch basin type
TSS_byepiweek <- ggplot(data = subset(waterqual_clean_wide, !is.na(TSS)), aes(x = Epiweek, y = TSS, color = CB_type, fill = CB_type)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(0.9), lwd = 0.25) +
  geom_boxplot(outlier.size = 0.25, lwd = 0.25, width = 0.75, outlier.shape = 19, position = position_dodge(0.9)) +  
  labs(x = "Epiweek", y = "Total Suspended Solids (mg/L)") +
  scale_color_manual(values = c("black", "black"), labels = legend_labs, name = "Catch Basin Type") +
  scale_fill_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  figtheme

TSS_byepiweek

# Plot average TSS concentration by Epiweek and catch basin type
TSS_plotting <- waterqual_clean_wide %>%
  filter(!is.na(TSS)) %>%
  group_by(Epiweek, CB_type) %>%
  summarise(
    mean_TSS_byepiweek = mean(TSS)
  )

TSS_avg_byepiweek <- ggplot(TSS_plotting, aes(x = CB_type, y = mean_TSS_byepiweek, fill = CB_type)) +
  labs(x = "Catch Basin Type", y = "Average TSS\nConcentration (mg/L)by Epiweek") +
  scale_fill_manual(values = c("grey27", "green4")) +
  geom_boxplot() +
  theme_classic()
TSS_avg_byepiweek

# Nonparametric t-test
TSS_t_test_data <- pivot_wider(TSS_plotting,
                               names_from = CB_type,
                               names_prefix = "TSS_",
                               values_from = mean_TSS_byepiweek)

# Checking normality - NH3 in CCBs
shapiro.test(TSS_t_test_data$TSS_C)
qqnorm(TSS_t_test_data$TSS_C)
qqline(TSS_t_test_data$TSS_C)

# Checking normality - NH3 in RGOs
shapiro.test(TSS_t_test_data$TSS_RGO)
qqnorm(TSS_t_test_data$TSS_RGO)
qqline(TSS_t_test_data$TSS_RGO)

# Paired Wilcoxon test
wilcox.test(x = TSS_t_test_data$TSS_C, y = TSS_t_test_data$TSS_RGO,
            alternative = "two.sided", paired = FALSE)

t.test(x = TSS_t_test_data$TSS_C, y = TSS_t_test_data$TSS_RGO,
       alternative = "two.sided", paired = FALSE)

# Averaged concentrations
TSS_plotting %>%
  group_by(CB_type) %>%
  summarise(
    mean_TSS = mean(mean_TSS_byepiweek),
    size_n = n(),
    se_TSS = sd(mean_TSS_byepiweek)/sqrt(size_n)
  )


NH3_N_PO4_TSS_plots_byepiweek <- NH3_N_byepiweek / PO4_byepiweek / TSS_byepiweek +
  plot_layout(guides = "collect", axis_titles = "collect") &
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", plot.tag.position = c(0.04, 1.01), plot.tag = element_text(size = 9, face = "bold", family = "HelveticaNeueforSAS")) 
NH3_N_PO4_TSS_plots_byepiweek

ggsave(
  plot = NH3_N_PO4_TSS_plots_byepiweek,
  filename = "figures/manuscript/Exp2_water_quality_all.jpeg",
  device ="jpeg",
  units = "cm",
  height = 20, width = 8.5, dpi = 300
)

ggsave(
  plot = NH3_N_PO4_TSS_plots_byepiweek,
  filename = "figures/manuscript/Exp2_water_quality_all.pdf",
  device ="pdf",
  units = "cm",
  height = 20, width = 8.5, dpi = 300
)

#### Raw winglength data ----
raw_wing <- read.csv("data/Rain gardens/pupal sampling/starvation_resistance/Aurora_2015_starvation_resistance_20241007.csv")

# Subset Cx. pipiens
wing_clean <- raw_wing %>%
  subset(Species == "Cx. pipiens")

# Merge water quality data to winglength data
wing_water_merge <- merge(x = wing_clean,
                          y = waterqual_clean_wide,
                          by.x = c("CBID", "DateColl", "Epiweek", "CB_Class"),
                          by.y = c("CB_ID_G", "Coll_Date", "Epiweek", "CB_type"))


#### Read in raw weather data ----
weather_raw <- read.csv("data/weather/Weather_daily_lags_Apr2013toOct2015.csv")

weatherClean <- weather_raw %>%
  select(DATE:NSA_PRCP_mm, NS_PRCP_mm_sum_Lag1:NSA_minT_d0_d1)

# Merge wing length and NH3 data to weather data
wingWaterWeather <- merge(x = wing_water_merge,
                          y = weatherClean,
                          by.x = "DateColl",
                          by.y = "DATE",
                          all.x = FALSE)

# Drop obs with missing winglength data
wingWaterWeatherFINAL <- wingWaterWeather %>%
  drop_na(Winglength_mm)

# Save final df with water quality, winglength, and weather data
write.csv(wingWaterWeatherFINAL,"data/Rain gardens/wingWaterWeatherFINAL.csv", row.names = FALSE)

# Read in combined data
wingWaterWeatherFINAL <- read.csv("data/Rain gardens/wingWaterWeatherFINAL.csv")

#### Correlation btwn wing length and water quality ----

## Females
females_wing_water <- wingWaterWeatherFINAL %>%
  subset(Sex == "F")

# Average female wing length of all females emerging from a unique catch basin by unique sampling date
females_wing_water_byCB <- females_wing_water %>%
  group_by(CBID, DateColl) %>%
  summarise(
    CB_Class = first(CB_Class),
    avg_wing = mean(Winglength_mm),
    NH3_N = first(NH3_N),
    PO4 = first(PO4),
    TSS = first(TSS),
    total_n = n()
  )

# Males
males_wing_water <- wingWaterWeatherFINAL %>%
  subset(Sex == "M")

# Average male wing length of all males emerging from a unique catch basin by unique samplign date
males_wing_water_byCB <- males_wing_water %>%
  group_by(CBID, DateColl) %>%
  summarise(
    CB_Class = first(CB_Class),
    avg_wing = mean(Winglength_mm),
    NH3_N = first(NH3_N),
    PO4 = first(PO4),
    TSS = first(TSS),
    total_n = n()
  )

## Correlations btwn wing length and NH3
# Scatterplot - females
females_NH3_scatter <- ggplot(females_wing_water_byCB, aes(x = NH3_N, y = avg_wing, color = CB_Class)) +
  geom_point() +
  scale_color_manual(values = c("grey27", "green4")) +
  labs(x = "Ammonia Nitrogen Concentration", y = "Wing Length (mm)", title = "Females", 
       subtitle = "Spearman's corr = 0.0243; p = 0.9145") +
  theme_classic()
females_NH3_scatter

# All catch basins - females
cor.test(x = females_wing_water_byCB$avg_wing,
         y = females_wing_water_byCB$NH3_N,
         method = "spearman")

# Scatterplot - males
males_NH3_scatter <- ggplot(males_wing_water_byCB, aes(x = NH3_N, y = avg_wing, color = CB_Class)) +
  geom_point() +
  scale_color_manual(values = c("grey27", "green4")) +
  labs(x = "Ammonia Nitrogen Concentration", y = "Wing Length (mm)", title = "Males", 
       subtitle = "Spearman's corr = 0.2288; p = 0.2822") +
  theme_classic()
males_NH3_scatter

# All catch basins - males
cor.test(x = males_wing_water_byCB$avg_wing,
         y = males_wing_water_byCB$NH3_N,
         method = "spearman")

## Correlations btwn wing length and PO4
# Scatterplot - females
females_PO4_scatter <- ggplot(females_wing_water_byCB, aes(x = PO4, y = avg_wing, color = CB_Class)) +
  geom_point() +
  scale_color_manual(values = c("grey27", "green4")) +
  labs(x = "PO4 Concentration", y = "Wing Length (mm)", title = "Females", 
       subtitle = "Spearman's corr = -0.3473; p = 0.1137") +
  theme_classic()
females_PO4_scatter

# All catch basins - females
cor.test(x = females_wing_water_byCB$avg_wing,
         y = females_wing_water_byCB$PO4,
         method = "spearman")

# Scatterplot - males
males_PO4_scatter <- ggplot(males_wing_water_byCB, aes(x = PO4, y = avg_wing, color = CB_Class)) +
  geom_point() +
  scale_color_manual(values = c("grey27", "green4")) +
  labs(x = "PO4 Concentration", y = "Wing Length (mm)", title = "Males", 
       subtitle = "Spearman's corr = -0.2504; p = 0.2368") +
  theme_classic()
males_PO4_scatter

# All CBs - males
cor.test(x = males_wing_water_byCB$avg_wing,
         y = males_wing_water_byCB$PO4,
         method = "spearman")

## Correlations btwn wing length and TSS
# Scatterplot - females
females_TSS_scatter <- ggplot(females_wing_water_byCB, aes(x = TSS, y = avg_wing, color = CB_Class)) +
  geom_point() +
  scale_color_manual(values = c("grey27", "green4")) +
  labs(x = "TSS Concentration", y = "Wing Length (mm)", title = "Females", 
       subtitle = "Spearman's corr = -0.1033; p = 0.6462") +
  theme_classic() 
females_TSS_scatter

# All catch basins - females
cor.test(x = females_wing_water_byCB$avg_wing,
         y = females_wing_water_byCB$TSS,
         method = "spearman")

# Scatterplot - males
males_TSS_scatter <- ggplot(males_wing_water_byCB, aes(x = TSS, y = avg_wing, color = CB_Class)) +
  geom_point() +
  scale_color_manual(values = c("grey27", "green4")) +
  labs(x = "TSS Concentration", y = "Wing Length (mm)", title = "Males", 
       subtitle = "Spearman's corr = 0.0670; p = 0.7555") +
  theme_classic() 
males_TSS_scatter

# All catch basins - males
cor.test(x = males_wing_water_byCB$avg_wing,
         y = males_wing_water_byCB$TSS,
         method = "spearman")

wing_water_corrs <- ((females_NH3_scatter/males_NH3_scatter) | (females_PO4_scatter/males_PO4_scatter) | (females_TSS_scatter/males_TSS_scatter)) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")
wing_water_corrs


# Because there are no correlations btwn water quality and wing length, water quality is not included as a covariate in the wing length analyses
