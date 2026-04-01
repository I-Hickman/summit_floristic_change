##Vegan package
#For an overview of codyn, see:
#Hallett LM, Jones SK, MacDonald AAA, Jones MB, Flynn DFB, Ripplinger J, Slaughter P, Gries C, Collins SL (2016) codyn: An R package of community dynamics metrics. Methods in Ecology and Evolution, 7(10):1146–1151. https://doi.org/10.1111/2041-210X.12569

#For a description of the newer spatial methods in codyn v2.x:
#Avolio ML, Carroll IT, Collins SL, Houseman GR, Hallett LM, Isbell F, Koerner SE, Komatsu KJ, Smith MD, Wilcox KR (2019) A comprehensive approach to analyzing community dynamics using rank abundance curves. Ecosphere, 10(10):e02881. https://doi.org/10.1002/ecs2.2881

#see for method write up:
#Cleland, Elsa E., Scott L. Collins, Timothy L. Dickson, Emily C. Farrer, Katherine L. Gross, Lau-
#reano A. Gherardi, Lauren M. Hallett, et al. (2013) "Sensitivity of grassland plant community
#composition to spatial vs. temporal variation in precipitation." Ecology 94, no. 8: 1687-96

library(codyn)
library(tidyverse)
library(ggplot2)
library(dplyr)

Florsitics <-read.csv("data/Raw floristic data.csv")
site_deets <- read.csv("data/Summit site details.csv")

#Convert to long form
TO_long <- Florsitics %>% 
  group_by(Year, Site, Transect, Aspect) %>% 
  gather(key = "Species", value = "Score", 5:107) 

#First need to get the average abundance of species across aspects
sp_abundnace_site <- TO_long %>% 
  #filter(Score != "0") %>% 
  group_by(Year, Site, Species) %>%
  dplyr::summarise(avg = mean(Score, na.rm = FALSE))

sp_abundnace_site_asp <- TO_long %>% 
  #filter(Score != "0") %>% 
  group_by(Year, Site, Aspect, Species) %>%
  dplyr::summarise(avg = mean(Score, na.rm = FALSE))


#######Community dynamics ########

###
# Step 1. Rank shift #####
###

###we use this to plot it on the graph!!
Rank_shift <- rank_shift(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site"
)


#Get mean and confidence intervals for each year
Rank_shift_mean <- Rank_shift %>%
  ungroup() %>%
  group_by(year_pair) %>%
  summarise(
    mean = mean(MRS, na.rm = TRUE),
    sd = sd(MRS, na.rm = TRUE),
    n = sum(!is.na(MRS)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

# Plot
# Extract the specific colors from the palette
colors <- paletteer::paletteer_d("ltc::trio3")

# Manually assign the extracted colors to the categories
custom_colors <- c("appearance" = colors[3], "disappearance" = colors[1], "total" = colors[2])

#Plot
RSM_plot <- ggplot(Rank_shift_mean, aes(x = year_pair, y = mean)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1, position = position_dodge(width = 0.5)) +
  labs(x = "Year-pair",
       y = "Mean rank shift") +
  theme_minimal() +
  ylim(0,15) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  #make text larger
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))

#make a new column for "Year" to use in the plot
Rank_shift_mean$Year = c("2012", "2017", "2022")
Rank_shift_mean$RAC = "MRS"
#Delete the year pair column
Rank_shift_mean$year_pair = NULL


#Save
#ggsave("figures/MRS_plot.png", Turnover_plot, width = 10, height = 5, units = "in", dpi = 300)



####
# Step 2: Rank Abundance Curve ####
###
#Change—tracks changes of a replicate through time
#Calculates changes in species richness, evenness, species’ ranks, gains, and losses for each replicate

#Community level
RAC <- RAC_change(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site")
  #reference.time = 2004)


#plot
RAC2 <- RAC %>% 
  ungroup() %>%
  group_by(Year, Year2, Site) %>%
  gather(key = "RAC", value = "Score", 4:8)

RAC2$Year2 = as.factor(RAC2$Year2)

#Get mean and confidence intervals for each year
RAC_mean <- RAC2 %>%
  ungroup() %>%
  group_by(Year2, RAC) %>%
  summarise(
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    n = sum(!is.na(Score)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

#Change the column name "Year2" to "Year"
RAC_mean <- RAC_mean %>% rename(Year = Year2)
Turnover <- read.csv("data/Turnover output.csv")

Turnover2 <- Turnover %>% 
  ungroup() %>%
  group_by(Year, Site) %>%
  gather(key = "Turnover", value = "Score", 3:5)

Turnover2 <- Turnover2 %>% filter(Turnover != "appearance" & Turnover != "disappearance")
Turnover2 <- Turnover2 %>% rename(RAC = Turnover)

Turnover_mean <- Turnover2 %>%
  ungroup() %>%
  group_by(Year, RAC) %>%
  summarise(
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    n = sum(!is.na(Score)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

Turnover_mean$Year = as.factor(Turnover_mean$Year)
RAC_mean$Year2 = as.factor(RAC_mean$Year2)

#Join
RAC_mean2 <- full_join(RAC_mean, Turnover_mean)
#RAC_mean3 <- full_join(RAC_mean2, Rank_shift_mean)

#Remove the richness_change rows
RAC_mean5 <- RAC_mean2 %>% filter(RAC != "richness_change")
RAC_mean6 <- RAC_mean5 %>% filter(RAC != "evenness_change")

#Plot
RAC_plot <- ggplot(RAC_mean6, aes(x = Year, y = mean, fill = RAC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1, position = position_dodge(width = 0.5)) +
  facet_wrap(~ RAC, labeller = as_labeller(c(
                                             rank_change = "Rank change", 
                                             gains = "Species gains", 
                                             losses = "Species losses",
                                             total = "Total turnover"))) +
  paletteer::scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  labs(x = "Year",
       y = " ",
       fill = "RAC") +
  theme_minimal() +
  theme_classic() +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

# Reorder the RAC factor
RAC_mean5 <- RAC_mean5 %>%
  mutate(RAC = factor(RAC, levels = c("gains", "losses", "total", "evenness_change", "rank_change")))

RAC_plot <- ggplot(RAC_mean5, aes(x = Year, y = mean, fill = RAC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1, position = position_dodge(width = 0.5)) +
  facet_wrap(~ RAC, labeller = as_labeller(c(
    gains = "Species gains", 
    losses = "Species losses",
    total = "Total turnover",
    evenness_change = "Evenness change", 
    rank_change = "Rank change"))) +
  paletteer::scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  labs(x = "Year",
       y = " ",
       fill = "RAC") +
  theme_minimal() +
  theme_classic() +
  ylim(-0.1, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

#Save
ggsave("figures/RAC_plot.png", RAC_plot, width = 7, height = 7, units = "in", dpi = 300)




### Site level
#Plot RAC$rank_change and facet by site
RAC$Year2 <- as.factor(RAC$Year2)
RAC_site <- left_join(RAC, site_deets, by = "Site")

#Add in turnover
#change Turnover2 column "Year" to "Year2"
Turnover$Year2 <- as.factor(Turnover$Year)
RAC_site2 <- left_join(RAC_site, Turnover, by = c("Year2", "Site"))

ggplot(RAC_site, aes(x = Year2, y = gains, group = Elevation)) +
  geom_line() +
  facet_wrap(~ Elevation) +
  theme_minimal() +
  theme_classic() +
  ylim(0, 0.5) +
  labs(x = "Year",
       y = "Rank change")

ggplot(RAC_site, aes(x = Year2, y = losses, group = Elevation)) +
  geom_line() +
  facet_wrap(~ Elevation) +
  theme_minimal() +
  theme_classic() +
  ylim(0, 0.5) +
  labs(x = "Year",
       y = "Species losses")

ggplot(RAC_site, aes(x = Year2, y = gains, group = Elevation)) +
  geom_line() +
  facet_wrap(~ Elevation) +
  theme_minimal() +
  theme_classic() +
  ylim(0, 0.5) +
  labs(x = "Year",
       y = "Species gains")

ggplot(RAC_site, aes(x = Year2, y = total, group = Elevation)) +
  geom_line() +
  facet_wrap(~ Elevation) +
  theme_minimal() +
  theme_classic() +
  ylim(0, 0.5) +
  labs(x = "Year",
       y = "Total turnover")



#Extract the specific colors from the palette
colors <- paletteer::paletteer_d("ltc::trio3")

# Manually assign the extracted colors to the categories
custom_colors <- c("gains" = colors[1], "losses" = colors[2], "total" = colors[3])


ggplot() +
  geom_line(data = RAC_site2, aes(x = Year2, y = losses, group = Elevation, color = "losses")) +
  geom_line(data = RAC_site2, aes(x = Year2, y = gains, group = Elevation, color = "gains")) +
  geom_line(data = RAC_site2, aes(x = Year2, y = total, group = Elevation, color = "total")) +
  facet_wrap(~ Elevation) +
  theme_minimal() +
  theme_classic() +
  ylim(0, 1) +
  labs(x = "Year",
       y = "Species turnover",
       color = "Turnover") +
  scale_color_manual(values = custom_colors)

#Save
ggsave("figures/Turnover_site_plot.png", width = 7, height = 5, units = "in", dpi = 300)


####
# Step 3: Rate of change ####
###

#tells us the rate of change for each site/fire history
Rate_change <- rate_change(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site"
)

rate_change_interval <- rate_change_interval(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site"
)


#plot the data
Rank_change_plot <- ggplot(data = rate_change_interval, aes(x =interval, y = distance)) + 
  geom_jitter(width = 0.05, height = 0, color = 'black', size = 2, alpha = 0.6) + 
  geom_smooth(method = "lm", colour = "black") +
  labs(x = "Time interval",
       y = "Euclidean distance") +
  theme_minimal() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))


#Save
ggsave("figures/Time_lag_plot.png", Rank_change_plot, width = 7, height = 5, units = "in", dpi = 300)

####
# Rate of change LM model #####
####

require(DHARMa)
require(lme4)
require(emmeans)
require(performance)
require(report)
require(ggplot2)
require(tidyverse)
library(cowplot)
library(arm)
library(modelsummary)

#Fit the models
m1 = lm(distance ~ interval, data = rate_change_interval)

### Step 5) Model validation
# Normality and homogeneity of variances
res = simulateResiduals(m1)
op = par(mfrow = c(1, 2), pty = 's')
plotQQunif(res)
plotResiduals(res, rate_change_interval$interval)
par(op) # Normality and homogeneity look good for m1 


# Check predictions
check_predictions(m1)
dev.off()


### Step 6) Model interpretation
# Model summary
summary(m1)
#none are significant

#### Step 7) Make a coefficient plot
#Formatting
m1_labels <- c(
  'interval' = 'Time interval',
  '(Intercept)' = 'Intercept')

b <- list(geom_vline(xintercept = 0, color = 'black', linetype = 'dashed'))

#Make plot
m1_plot <- modelplot(m1, coef_map = m1_labels, background = b) +
  labs(x = 'Coefficients') +
  theme_cowplot() 



#### Plotting for each site ######


RAC

### Site level
#Plot RAC$rank_change and facet by site
RAC$Year2 <- as.factor(RAC$Year2)
RAC_site <- left_join(RAC, site_deets, by = "Site")

#Add in turnover
#change Turnover2 column "Year" to "Year2"
Turnover$Year2 <- as.factor(Turnover$Year)
RAC_site2 <- left_join(RAC_site, Turnover, by = c("Year2", "Site"))
new_ele <- read.csv("data/Site evelvations.csv")
RAC_site2 <- left_join(RAC_site2, new_ele, by = "Site")

#Extract the specific colors from the palette
colors <- paletteer::paletteer_d("ltc::trio3")

# Manually assign the extracted colors to the categories
custom_colors <- c("gains" = colors[1], "losses" = colors[2], "total" = colors[3])


ggplot() +
  geom_line(data = RAC_site2, aes(x = Year2, y = losses, group = Elevation_new, color = "losses")) +
  geom_line(data = RAC_site2, aes(x = Year2, y = gains, group = Elevation_new, color = "gains")) +
  geom_line(data = RAC_site2, aes(x = Year2, y = total, group = Elevation_new, color = "total")) +
  facet_wrap(~ Elevation_new) +
  theme_minimal() +
  theme_classic() +
  ylim(0, 1) +
  labs(x = "Year",
       y = "Species turnover",
       color = "Turnover") +
  scale_color_manual(values = custom_colors)

#Save
ggsave("figures/Turnover_site_plot.png", width = 7, height = 5, units = "in", dpi = 300)



#Rank change
ggplot() +
  geom_line(data = RAC_site2, aes(x = Year2, y = rank_change, group = Elevation_new)) +
  facet_wrap(~ Elevation_new) +
  theme_minimal() +
  theme_classic() +
  ylim(0, 0.5) +
  labs(x = "Year",
       y = "Rank change")

#Save
ggsave("figures/Rank_change_site_plot.png", width = 7, height = 5, units = "in", dpi = 300)




####  Change in composition ######

multivariate_change <- multivariate_change(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site"
)


#Join together multivariate_change with 
RAC_df <- RAC_site2 %>%
  select(Site, Year.x, Year2, gains, losses, total, rank_change)

#change Year.x to Year
RAC_df <- RAC_df %>%
  rename(Year = Year.x)

#join with multivariate_change
multivariate_change$Year2 <- as.factor(multivariate_change$Year2)
multivariate_change$Year <- as.factor(multivariate_change$Year)
RAC_df$Year <- as.factor(RAC_df$Year)

RAC_df <- left_join(RAC_df, multivariate_change, by = c("Year", "Year2"))

## Check correlation between variables
library(GGally)
library(dplyr)
library(ggplot2)

colors <- paletteer::paletteer_d("ltc::trio3")

custom_colors2 <- c("gains" = colors[1], "losses" = colors[2], 
                    "total" = colors[3], "rank_change" = colors[4], 
                    "composition_change" = colors[5])


RAC_df %>%
  dplyr::select(-c(Site, Year, Year2, composition_change)) %>%
  GGally::ggpairs(
    columnLabels = c("Gains", "Losses", "Total", "Rank Change", "Composition Change"),
    upper = list(continuous = wrap("cor", size = 5)),
    lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.5)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.5))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 7)
  ) +
  theme_minimal() +
  theme_classic()


#save the plot
ggsave("figures/Correlation_plot.png", width = 8, height = 7, units = "in", dpi = 300)
