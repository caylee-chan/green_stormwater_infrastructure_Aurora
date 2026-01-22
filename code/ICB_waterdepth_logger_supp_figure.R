# Logger data in 10 CCBs and 10 ICBs
# Caylee Chan
# Updated 22 Jan 2026
# Notes:


# Libraries:
library(dplyr)
library(ggplot2)
library(tidyr)
library(extrafont)

# Setwd
setwd("C:/Users/cayle.LAPTOP-QMLQMN4A/OneDrive - University of Illinois - Urbana/green_stormwater_infrastructure_Aurora/green_stormwater_infrastructure_Aurora")


# Set theme
figtheme <- theme(
  strip.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = 'black', linewidth = 0.3), 
  axis.title.x = element_text(size = 10, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.title.y = element_text(size = 10, color="black", face = "bold", family = "HelveticaNeueforSAS"), 
  axis.text.x = element_text(size = 10, color="black", family = "HelveticaNeueforSAS"), 
  axis.text.y = element_text(size = 10, color="black", family = "HelveticaNeueforSAS"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  legend.text = element_text(size = 10, family = "HelveticaNeueforSAS"), 
  legend.position = "top",
  legend.title = element_text(size = 10, color="black", family = "HelveticaNeueforSAS"))

# Set legend colors
legend_colors <- c("grey27", "green4")
legend_labs <- c("Conventional", "Infiltration")

# Raw data
logger_14 <- read.csv("data/Onset records/2014/water depth/All_2014_IPCA_T_Depth_15mininterval_Andrew_15Aug2025.csv")

str(logger_14$Date)
logger_14$Date <- as.Date(logger_14$Date, "%d-%b-%y")

#### Water temp ----

# Pivot longer - temp
logger_14_long_temp <- logger_14 %>%
  pivot_longer(
    cols = matches("_Temp_C$"), # Water temp columns
    names_to = "CB_ID", # Make water temp column names values of CB_ID
    values_to = "watertemp_C", # Convert values in water temp columns to watertemp_C
  ) %>%
  subset(!is.na(watertemp_C)) %>% # Remove NA water temp vals
  mutate(CB_ID = sub("_Temp_C", "", CB_ID)) %>% # Remove "_Temp_C" from CB_ID col name
  select("Date", "Hour", "Time", "CB_ID", "watertemp_C") %>%
  mutate(CB_type = ifelse(CB_ID %in% c("CB234", "CB249", "CB260", "CB263", "CB271", "CB287", "CB292", "CB505", "CB514", "CB534"), "ICB", "CCB"))  # Create col for catch basin type

# Calculate average temp for each day in CCBs and ICBs - by catch basin ID
avg_temp_per_day_by_CB_ID <- logger_14_long_temp %>%
  group_by(Date, CB_ID, CB_type) %>%
  summarise(
    avg_temp_C_byID = mean(watertemp_C, na.rm = TRUE)
  )

# Calculate average temp for each day in CCBs and ICBs - by catch basin TYPE
avg_temp_per_day_by_CB_TYPE <- avg_temp_per_day_by_CB_ID %>%
  group_by(Date, CB_type) %>%
  summarise(
    avg_temp_C_byTYPE = mean(avg_temp_C_byID, na.rm = TRUE)
  )

# Mean temp in ICBs and CCBs
mean(avg_temp_per_day_by_CB_TYPE$avg_temp_C_byTYPE[avg_temp_per_day_by_CB_TYPE$CB_type == "ICB"], na.rm = TRUE)
mean(avg_temp_per_day_by_CB_TYPE$avg_temp_C_byTYPE[avg_temp_per_day_by_CB_TYPE$CB_type == "CCB"], na.rm = TRUE)

str(avg_temp_per_day_by_CB_TYPE$CB_type)
avg_temp_per_day_by_CB_TYPE$CB_type <- factor(avg_temp_per_day_by_CB_TYPE$CB_type)


# Plot date on x-axis, water temp on y-axis, colored by catch basin 
water_temp <- ggplot(avg_temp_per_day_by_CB_TYPE, aes(x = Date, y = avg_temp_C_byTYPE, na.rm = TRUE, color = CB_type)) +
  geom_line(linewidth = 1) + #geom_point(size = 1.5) +
  scale_color_manual(values = legend_colors, labels = legend_labs, name = "Catch Basin Type") +
  labs(x = "Date", y = "Water Temperature (C)") +
  scale_y_continuous(limits = c(14, 26), expand = c(0,0)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  figtheme

water_temp


#### Water depth ----

# Pivot longer - depth
logger_14_long_depth <- logger_14 %>%
  pivot_longer(
    cols = matches("_Sensor_Depth_m"), # Sensor depth columns
    names_to = "CB_ID", # Make sensor depth column names values of CB_ID
    values_to = "sensordepth_m", # Convert values in sensor depth cols to sensordepth_m
  ) %>%
  subset(!is.na(sensordepth_m)) %>% # Remove NA sensor depth vals
  mutate(CB_ID = sub("_Sensor_Depth_m", "", CB_ID)) %>% # Remove "_Sensor_Depth_m" from CB_ID col name
  select("Date", "Hour", "Time", "CB_ID", "sensordepth_m") %>%
  mutate(CB_type = ifelse(CB_ID %in% c("CB234", "CB249", "CB260", "CB263", "CB271", "CB287", "CB292", "CB505", "CB514", "CB534"), "ICB", "CCB")) # Create col for catch basin type

# Calculate average depth for each day in CCBs and ICBs - by catch basin ID
avg_depth_per_day_by_CB_ID <- logger_14_long_depth %>%
  group_by(Date, CB_ID, CB_type) %>%
  summarise(
    avg_depth_m_byID = mean(sensordepth_m, na.rm = TRUE)
  )

# Calculate average depth for each day in CCBs and ICBs - by catch basin TYPE
avg_depth_per_day_by_CB_TYPE <- avg_depth_per_day_by_CB_ID %>%
  group_by(Date, CB_type) %>%
  summarise(
    avg_depth_m_byTYPE = mean(avg_depth_m_byID, na.rm = TRUE)
  )

# Average sensor depth for ICBs and CCBs
mean(avg_depth_per_day_by_CB_TYPE$avg_depth_m_byTYPE[avg_depth_per_day_by_CB_TYPE$CB_type == "ICB"], na.rm = TRUE)
mean(avg_depth_per_day_by_CB_TYPE$avg_depth_m_byTYPE[avg_depth_per_day_by_CB_TYPE$CB_type == "CCB"], na.rm = TRUE)

str(avg_depth_per_day_by_CB_TYPE$CB_type)
avg_depth_per_day_by_CB_TYPE$CB_type <- factor(avg_depth_per_day_by_CB_TYPE$CB_type)


# Plot date on x-axis, sensor depth on y-axis, colored by catch basin type
sensordepth <- ggplot(avg_depth_per_day_by_CB_TYPE, aes(x = Date, y = avg_depth_m_byTYPE, na.rm = TRUE, color = CB_type)) +
  geom_line(linewidth = 0.7)  +
  scale_color_manual(values = legend_colors, labels = legend_labs, name = "Catch\nBasin\nType") +
  labs(x = "Date", y = "Sensor Depth (m)") +
  scale_y_continuous(breaks = seq(0, 0.25, 0.05), expand = c(0,0)) +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  figtheme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sensordepth

ggsave(sensordepth,
       filename = "figures/manuscript/sensordepth.jpeg",
       device = "jpeg",
       units = "mm", 
       height = 120, width = 88, dpi = 300)


