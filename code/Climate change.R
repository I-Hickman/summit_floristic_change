

################ CLIMATE ##########
library(usethis)
library(devtools)
library(microclima)
library(NicheMapR)
library(dplyr)
library(future)
library(maps)
library(RNCEP)
library(future.apply)
library(cowplot)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tibble)
library(beepr)
library(tidyr)
suppressMessages(library(lubridate))
library(cropgrowdays)


#Read in file
sites <- read.csv("data/Summit site details.csv")


############# BASELINE #####################
#Download baseline climate data (1961-1990)
#Run loop through all columns
plan(multisession, workers = 4)
result_list <- future.apply::future_lapply(future.seed = NULL, 1:nrow(sites), function(x) {
  nm <- cropgrowdays::get_multi_silodata(
    latitude = sites$Lat[x],
    longitude = sites$Long[x],
    Sitename = sites$Site[x],
    START = "19610101", FINISH = "19900401",
    email = "irishickman11@hotmail.com") 
  out <- as.data.frame(nm)
  out
})  
result_list2 <- do.call(rbind, result_list)

#Save
write.csv(result_list2, "data/Baseline climate data.csv", row.names = FALSE)
#Select columns
baseline_clim <- result_list2[,c("Sitename", "date_met", "year", "day", "mint", "maxt", "rain" )]
#rename day column to julian day
baseline_clim <- baseline_clim %>% rename(julian_day = day)
#Separate the date_met column into year, month and day
baseline_clim2 <- baseline_clim %>% 
  separate(date_met, into = c("year", "month", "day"), sep = "-") %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day))
glimpse(baseline_clim2)

##### 1. Annual precipitation ######
## Step 1) Calculate the annual precipitation
baseline_precip <- baseline_clim2 %>% 
  group_by(Sitename, year) %>%
  summarise(sum_rain = sum(rain, na.rm = TRUE))

## Step 2) Calculate the mean annual precipitation for each site
baseline_mean_precip <- baseline_precip %>% 
  group_by(Sitename) %>%
  summarise(mean_rain = mean(sum_rain, na.rm = TRUE))


##### 2. Warming ######
## Step 1) get the mean annual maxt, and mint
baseline_maxt <- baseline_clim2 %>% 
  group_by(Sitename, year) %>%
  summarise(mean_maxt = mean(maxt, na.rm = TRUE))
baseline_mint <- baseline_clim2 %>% 
  group_by(Sitename, year) %>%
  summarise(mean_mint = mean(mint, na.rm = TRUE))

## Step 2) Calculate the mean annual precipitation
baseline_mean_annual_maxt <- baseline_maxt %>% 
  group_by(Sitename) %>%
  summarise(mean_maxt = mean(mean_maxt, na.rm = TRUE))

baseline_mean_annual_mint <- baseline_mint %>% 
  group_by(Sitename) %>%
  summarise(mean_mint = mean(mean_mint, na.rm = TRUE))

#### Dataframes
#Precipitation
baseline_mean_precip$time <- "baseline"
#Max temp
baseline_mean_annual_maxt$time <- "baseline"
#Min temp
baseline_mean_annual_mint$time <- "baseline"


############# CURRENT PERIOD #####################
#Download baseline climate data (1961-1990)
#Run loop through all columns
plan(multisession, workers = 4)
result_list <- future.apply::future_lapply(future.seed = NULL, 1:nrow(sites), function(x) {
  nm <- cropgrowdays::get_multi_silodata(
    latitude = sites$Lat[x],
    longitude = sites$Long[x],
    Sitename = sites$Site[x],
    START = "20030101", FINISH = "20220401",
    email = "irishickman11@hotmail.com") 
  out <- as.data.frame(nm)
  out
})  
result_list2 <- do.call(rbind, result_list)

#Save
write.csv(result_list2, "data/Current climate data.csv", row.names = FALSE)
#Select columns
current_clim <- result_list2[,c("Sitename", "date_met", "year", "day", "mint", "maxt", "rain" )]
#rename day column to julian day
current_clim <- current_clim %>% rename(julian_day = day)
#Separate the date_met column into year, month and day
current_clim2 <- current_clim %>% 
  separate(date_met, into = c("year", "month", "day"), sep = "-") %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day))


##### 1. Annual precipitation ######
## Step 1) Calculate the annual precipitation
current_precip <- current_clim2 %>% 
  group_by(Sitename, year) %>%
  summarise(sum_rain = sum(rain, na.rm = TRUE))

## Step 2) Calculate the mean annual precipitation for each site
current_mean_precip <- current_precip %>% 
  group_by(Sitename) %>%
  summarise(mean_rain = mean(sum_rain, na.rm = TRUE))


##### 2. Warming ######
## Step 1) get the mean annual maxt, and mint
current_maxt <- current_clim2 %>% 
  group_by(Sitename, year) %>%
  summarise(mean_maxt = mean(maxt, na.rm = TRUE))
current_mint <- current_clim2 %>% 
  group_by(Sitename, year) %>%
  summarise(mean_mint = mean(mint, na.rm = TRUE))

## Step 2) Calculate the mean annual precipitation
current_mean_annual_maxt <- current_maxt %>% 
  group_by(Sitename) %>%
  summarise(mean_maxt = mean(mean_maxt, na.rm = TRUE))

current_mean_annual_mint <- current_mint %>% 
  group_by(Sitename) %>%
  summarise(mean_mint = mean(mean_mint, na.rm = TRUE))

#### Dataframes
#Precipitation
current_mean_precip$time <- "current"
#Max temp
current_mean_annual_maxt$time <- "current"
#Min temp
current_mean_annual_mint$time <- "current"


#Join all dataframes together
precip <- rbind(baseline_mean_precip, current_mean_precip)
maxt <- rbind(baseline_mean_annual_maxt, current_mean_annual_maxt)
mint <- rbind(baseline_mean_annual_mint, current_mean_annual_mint)

########### CHANGE FROM BASELINE ########
#### 1. precipitation #####
#Calculate the change in precipitation
precip_change <- precip %>% 
  group_by(Sitename) %>%
  summarise(change_precip = mean(mean_rain[time == "current"], na.rm = TRUE) - mean(mean_rain[time == "baseline"], na.rm = TRUE))

#### 2. Warming #####
#Calculate the change in maxt
maxt_change <- maxt %>% 
  group_by(Sitename) %>%
  summarise(change_maxt = mean(mean_maxt[time == "current"], na.rm = TRUE) - mean(mean_maxt[time == "baseline"], na.rm = TRUE))

#Calculate the change in mint
mint_change <- mint %>% 
  group_by(Sitename) %>%
  summarise(change_mint = mean(mean_mint[time == "current"], na.rm = TRUE) - mean(mean_mint[time == "baseline"], na.rm = TRUE))

#Join all dataframes together
climate_change <- left_join(precip_change, maxt_change, by = "Sitename") %>%
  left_join(mint_change, by = "Sitename")

#Change Sitename to "Site"
climate_change <- climate_change %>% rename(Site = Sitename)

#Join elevation
climate_change <- left_join(climate_change, sites, by = "Site")


#######Plot change from baseline over elevation
library(paletteer)
library(ggplot2)
library(paletteer)

# Extract specific colors from the palette
colors <- paletteer::paletteer_d("nationalparkcolors::Acadia")

# Manually assign the extracted colors to the categories
color_precip <- colors[4]
color_maxt <- colors[1]
color_mint <- colors[3]

# Example usage in ggplot
# Precipitation
a <- ggplot(climate_change, aes(Elevation, change_precip)) +
  geom_point(size = 2, color = color_precip) +
  geom_smooth(method = "lm", se = FALSE, color = color_precip) +
  labs(x = "Elevation (m)", y = "Change in precipitation (mm)") +
  theme_classic() +
  ylim(NA, 0) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

# Max temperature
b <- ggplot(climate_change, aes(Elevation, change_maxt)) +
  geom_point(size = 2, color = color_maxt) +
  geom_smooth(method = "lm", se = FALSE, color = color_maxt) +
  labs(x = "Elevation (m)", y = "Change in max temperature (°C)") +
  theme_classic() +
  ylim(0, 1.5) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

# Min temperature
c <- ggplot(climate_change, aes(Elevation, change_mint)) +
  geom_point(size = 2, color = color_mint) +
  geom_smooth(method = "lm", se = FALSE, color = color_mint) +
  labs(x = "Elevation (m)", y = "Change in min temperature (°C)") +
  theme_classic() +
  ylim(0, 1.5) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))


#Stick together
library(gridExtra)
plot <- grid.arrange(a, b, c, ncol = 3)

#save
ggsave("figures/climate_change.png", plot, width = 12, height = 4, units = "in", dpi = 300)



