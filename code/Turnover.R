

##########################
#Turnover over time#
##########################
library(codyn)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(brms)

##TABLE OF CONTENTS
####
####
#### Step 1. Calculate turnover
#### Step 2. GLMM over time
###
##
##


##########################
# Step 1. Calculate turnover####
##########################

Florsitics <-read.csv("data/Raw floristic data.csv")

#Convert to long form
TO_long <- Florsitics %>% 
  group_by(Year, Site, Transect, Aspect) %>% 
  gather(key = "Species", value = "Score", 5:107) 

#First need to get the average abundance of species across aspects
sp_abundnace_site <- TO_long %>% 
  #filter(Score != "0") %>% 
  group_by(Year, Site, Species) %>%
  dplyr::summarise(avg = mean(Score, na.rm = FALSE))


#Species turnover
Turnover_total <- turnover(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site",
  metric = "total"
)

Turnover_app <- turnover(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site",
  metric = "appearance"
)

Turnover_disap <- turnover(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site",
  metric = "disappearance"
)


#join all together and add in site info
Turnover_join <- left_join(Turnover_total, Turnover_disap, by=c("Year", "Site")) %>% 
  left_join(Turnover_app, by=c("Year", "Site")) %>% 
  select(Year, Site, disappearance, appearance, total)

#Save
write.csv(Turnover_join, "data/Turnover output.csv", row.names = FALSE)

##########################################


# Step 2. Plot turnover ####
Turnover <- read.csv("data/Turnover output.csv")
Turnover2 <- Turnover %>% 
  ungroup() %>%
  group_by(Year, Site) %>%
  gather(key = "Turnover", value = "Score", 3:5)

Turnover2$Year = as.factor(Turnover2$Year)

#Get mean and confidence intervals for each year
Turnover_mean <- Turnover2 %>%
  ungroup() %>%
  group_by(Year, Turnover) %>%
  summarise(
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    n = sum(!is.na(Score)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

#Plot turnover and facet by Turnover
# Extract the specific colors from the palette
colors <- paletteer::paletteer_d("ltc::trio3")

# Manually assign the extracted colors to the categories
custom_colors <- c("appearance" = colors[3], "disappearance" = colors[1], "total" = colors[2])

#Plot
Turnover_plot <- ggplot(Turnover_mean, aes(x = Year, y = mean, color = Turnover)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Turnover, labeller = as_labeller(c(appearance = "Appearance", disappearance = "Disappearance", total = "Total Turnover"))) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Year",
       y = "Mean",
       color = "Turnover") +
  theme_minimal() +
  ylim(0, 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  #make text larger
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))
  
#Save
ggsave("figures/Turnover_plot.png", Turnover_plot, width = 10, height = 5, units = "in", dpi = 300)


#### ENDDD ######




##########################
# Step 3. GLMM ####
##########################
library(lme4)
library(mgcv)
library()
library(mgcv)
library(DHARMa)
library(performance)
library(adespatial)
library(MASS)
library(report) 
library(ggeffects)
library(cowplot)
library(sjPlot)
library(ggeffects)
library(glmmTMB)
library(webshot)
library(vegan)
library(broom.mixed)
library(jtools)
library(effects)
library(gridExtra)
library(ggplot2)
library(boot)
library(scales)
library(dplyr)
library(tidyverse)

##### Step 1: read in data, explore and scale
#Read in data
Turnover <- read.csv("data/Turnover output.csv")
glimpse(Turnover)
summary(Turnover)

#Plot
plot(Turnover$total~Turnover$Year)
hist(Turnover$total)
hist(Turnover$appearance)
hist(Turnover$disappearance)

#Transpose total, appearance and disappearance
Turnover2 <- Turnover %>% 
  ungroup() %>%
  group_by(Year, Site) %>%
  gather(key = "Turnover", value = "Score", 3:5)

#scale
Turnover2$Year = as.factor(Turnover2$Year)

#############################################
###     Total turnover model #####
########################################

TotalMod <-  brms::brm(Score ~ Turnover*Year + (1 | Site),
                       family = gaussian,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       control = list(adapt_delta = 0.8, max_treedepth=15),
                       data = Turnover2)

######## Step 1: Validate model 
summary(TotalMod)
plot(TotalMod)
pp_check(TotalMod) 
prior_summary(TotalMod)

##### Step 2: Fixed effects coefficient plot
Total_fixed <- as.data.frame(fixef(TotalMod))
Total_fixed2 <- rownames_to_column(Total_fixed, var = "Variable")
Total_fixed2$Variable <- factor(Total_fixed2$Variable,
                                 levels = c("SR.std:Year.std", "SR.std", "Intercept"))

custom_labels <- c("SR.std" = "SR",
                   "SR.std:Year.std" = "SR:Year",
                   "Intercept" = "Intercept")

Total_fixed_plot <- ggplot(Total_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")



#### Step 5) Random effects

Total_rand <- as.data.frame(ranef(TotalMod))
Total_rand2 <- rownames_to_column(Total_rand, var = "Variable")
Total_rand2 <- Total_rand2[order(Total_rand2$Site.Estimate.Intercept), ]

Total_rand_plot <- ggplot(Total_rand2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

ggsave("output/Turnover/Bayesian/ Total turnover (Site random effects).png", 
       plot = Total_rand_plot, width = 10, height = 7.07) 



#### Step 4) Plot predictions
#Total turnover vs solar radation
mean(Turnover$SR) #918.3846
sd(Turnover$SR) #40.14328
range(Turnover$SR) #830.2638 986.1972
#New dataframe
newdat_total <- data.frame(Site = "average site",
                           SR = seq(830, 985, by = 5), 
                           Year.F = 0) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(TotalMod, newdat= newdat_total, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_total$mean <- mean
newdat_total$lci <- lci
newdat_total$uci <- uci


##Plot predictions
Total_plot_SR <- ggplot(newdat_total, aes(x = SR, y = mean)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_point(data = Turnover, aes(x = SR, y = total), 
             shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  geom_line(size = 0.75) +
  theme_cowplot()+
  ylim(0,0.8) +
  ylab(expression("Total turnover"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 




### CHANGES OVER TIME
#What is 2004 and 2022 predictions against solar radiation
mean(Turnover$Year) #2019.5
sd(Turnover$Year) #2.511236
#Year.std = (Year - 2019.5)/2.511236
#2017 = -0.9955257
#2022 = 0.9955257

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
## 1. 2004 prediction
newdat_total_17 <- data.frame(Site = "average site",
                              SR = seq(830, 985, by = 5), 
                              Year.std = -0.9955257) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(TotalMod, newdat= newdat_total_17, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_total_17$mean <- mean
newdat_total_17$lci <- lci
newdat_total_17$uci <- uci


## 2. 2022 prediction
newdat_total_22 <- data.frame(Site = "average site",
                              SR = seq(830, 985, by = 5), 
                              Year.std = 0.9955257) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(TotalMod, newdat= newdat_total_22, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_total_22$mean <- mean
newdat_total_22$lci <- lci
newdat_total_22$uci <- uci


##Plot predictions
total_plot <- ggplot(newdat_total_22, aes(x = SR, y = mean)) + 
  geom_path(color = "darkorange") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkorange") +
  geom_path(data = newdat_total_17, aes(x = SR, y = mean), color = "purple") +
  geom_ribbon(data = newdat_total_17, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "purple") +
  theme_cowplot()+
  ylim(0,0.8) +
  ylab(expression("Total turnover"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 


#PLot all togeher
plot_total <- grid.arrange(Total_fixed_plot,
                           Total_plot_SR,
                           total_plot, ncol = 3)

ggsave("output/Turnover/Bayesian/Turnover plot.png", 
       plot = plot_forb2, width = 15, height = 3.5)





#############################################
###    #Appearances turnover model  #####
########################################

AppearMod <- brms::brm(appearance ~ SR.std + SR.std:Year.std + (1 | Site),
                       family = gaussian,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       control = list(adapt_delta = 0.8, max_treedepth=10),
                       data = Turnover)

######## Step 1: Validate model 
summary(AppearMod)
plot(AppearMod)
pp_check(AppearMod) 

##### Step 2: Fixed effects coefficient plot
Appear_fixed <- as.data.frame(fixef(AppearMod))
Appear_fixed2 <- rownames_to_column(Appear_fixed, var = "Variable")
Appear_fixed2$Variable <- factor(Appear_fixed2$Variable,
                                levels = c("SR.std:Year.std", "SR.std", "Intercept"))


Appear_fixed_plot <- ggplot(Appear_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")



#### Step 5) Random effects

Appear_rand <- as.data.frame(ranef(AppearMod))
Appear_rand2 <- rownames_to_column(Appear_rand, var = "Variable")
Appear_rand2 <- Appear_rand2[order(Appear_rand2$Site.Estimate.Intercept), ]

Appear_rand_plot <- ggplot(Appear_rand2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

ggsave("output/Turnover/Bayesian/ Appearances turnover (Site random effects).png", 
       plot = Appear_rand_plot, width = 10, height = 7.07) 



#### Step 4) Plot predictions
# Appearances vs solar radation
mean(Turnover$SR) #918.3846
sd(Turnover$SR) #40.14328
range(Turnover$SR) #830.2638 986.1972
#New dataframe
newdat_appear <- data.frame(Site = "average site",
                           SR = seq(830, 985, by = 5), 
                           Year.std = 0) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(AppearMod, newdat= newdat_appear, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_appear$mean <- mean
newdat_appear$lci <- lci
newdat_appear$uci <- uci


##Plot predictions
Appear_plot_SR <- ggplot(newdat_appear, aes(x = SR, y = mean)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_point(data = Turnover, aes(x = SR, y = appearance), 
             shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  geom_line(size = 0.75) +
  theme_cowplot()+
  ylim(0,0.8) +
  ylab(expression("Appearances"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 



### CHANGES OVER TIME
## 1. 2004 prediction
newdat_appear_17 <- data.frame(Site = "average site",
                              SR = seq(830, 985, by = 5), 
                              Year.std = -0.9955257) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(AppearMod, newdat= newdat_appear_17, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_appear_17$mean <- mean
newdat_appear_17$lci <- lci
newdat_appear_17$uci <- uci


## 2. 2022 prediction
newdat_appear_22 <- data.frame(Site = "average site",
                              SR = seq(830, 985, by = 5), 
                              Year.std = 0.9955257) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(AppearMod, newdat= newdat_appear_22, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_appear_22$mean <- mean
newdat_appear_22$lci <- lci
newdat_appear_22$uci <- uci


##Plot predictions
appear_plot <- ggplot(newdat_appear_22, aes(x = SR, y = mean)) + 
  geom_path(color = "darkorange") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkorange") +
  geom_path(data = newdat_appear_17, aes(x = SR, y = mean), color = "purple") +
  geom_ribbon(data = newdat_appear_17, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "purple") +
  theme_cowplot()+
  ylim(0,0.8) +
  ylab(expression("Appearances"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 


#PLot all togeher
plot_appearance <- grid.arrange(Appear_fixed_plot,
                           Appear_plot_SR,
                           appear_plot, ncol = 3)

ggsave("output/Turnover/Bayesian/Appearances plot.png", 
       plot = plot_appearance, width = 15, height = 3.5)






#############################################
###    #Disappearances turnover model  #####
########################################

DisappMod <-brms::brm(disappearance ~ SR.std + SR.std:Year.std + (1 | Site),
                        family = gaussian,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        control = list(adapt_delta = 0.8, max_treedepth=10),
                        data = Turnover)
  

######## Step 1: Validate model 
summary(DisappMod)
plot(DisappMod)
pp_check(DisappMod) 

##### Step 2: Fixed effects coefficient plot
Disappear_fixed <- as.data.frame(fixef(DisappMod))
Disappear_fixed2 <- rownames_to_column(Disappear_fixed, var = "Variable")
Disappear_fixed2$Variable <- factor(Disappear_fixed2$Variable,
                                 levels = c("SR.std:Year.std", "SR.std", "Intercept"))


Disappear_fixed_plot <- ggplot(Disappear_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")



#### Step 5) Random effects

Disappear_rand <- as.data.frame(ranef(DisappMod))
Disappear_rand2 <- rownames_to_column(Disappear_rand, var = "Variable")
Disappear_rand2 <- Appear_rand2[order(Disappear_rand2$Site.Estimate.Intercept), ]

Disappear_rand_plot <- ggplot(Disappear_rand2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")
Disappear_rand_plot

ggsave("output/Turnover/Bayesian/Disappearances turnover (Site random effects).png", 
       plot = Disappear_rand_plot, width = 10, height = 7.07) 



#### Step 4) Plot predictions
# 1. Disappearances vs solar radation
#New dataframe
newdat_disappear <- data.frame(Site = "average site",
                            SR = seq(830, 985, by = 5), 
                            Year.std = 0) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(DisappMod, newdat= newdat_disappear, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_disappear$mean <- mean
newdat_disappear$lci <- lci
newdat_disappear$uci <- uci


##Plot predictions
Disappear_plot_SR <- ggplot(newdat_disappear, aes(x = SR, y = mean)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_point(data = Turnover, aes(x = SR, y = disappearance), 
             shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  geom_line(size = 0.75) +
  theme_cowplot()+
  ylim(0,0.8) +
  ylab(expression("Disappearances"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 



### CHANGES OVER TIME

## 1. 2004 prediction
newdat_disappear_17 <- data.frame(Site = "average site",
                               SR = seq(830, 985, by = 5), 
                               Year.std = -0.9955257) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(DisappMod, newdat= newdat_disappear_17, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_disappear_17$mean <- mean
newdat_disappear_17$lci <- lci
newdat_disappear_17$uci <- uci


## 2. 2022 prediction
newdat_disappear_22 <- data.frame(Site = "average site",
                               SR = seq(830, 985, by = 5), 
                               Year.std = 0.9955257) %>% 
  dplyr::mutate(SR.std = (SR - 918.3846)/40.14328)
#predictions
predictions<- posterior_epred(DisappMod, newdat= newdat_disappear_22, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_disappear_22$mean <- mean
newdat_disappear_22$lci <- lci
newdat_disappear_22$uci <- uci


##Plot predictions
disappear_plot <- ggplot(newdat_disappear_22, aes(x = SR, y = mean)) + 
  geom_path(color = "darkorange") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkorange") +
  geom_path(data = newdat_disappear_17, aes(x = SR, y = mean), color = "purple") +
  geom_ribbon(data = newdat_disappear_17, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "purple") +
  theme_cowplot()+
  ylim(0,0.8) +
  ylab(expression("Disappearances"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 


#PLot all togeher
plot_disappearance <- grid.arrange(Disappear_fixed_plot,
                                Disappear_plot_SR,
                                disappear_plot, ncol = 3)

ggsave("output/Turnover/Bayesian/Appearances plot.png", 
       plot = plot_appearance, width = 15, height = 3.5)





##### PLOT ALL LIFEFORM PLOTS ####

plot_all <- grid.arrange(plot_total, 
                         plot_appearance,
                         plot_disappearance,
                         ncol = 1)

ggsave("output/Turnover/Bayesian/Turnover over time:SR plot.png", 
       plot = plot_all, width = 15, height = 8)


###############################
##### APPENDICES ######
###############################

#Disappearances
preds <- as.data.frame(fitted(DisappMod))
plot(preds$Estimate ~ DisappMod$data$disappearance)
abline(0, 1, col= 'red')

#Graminoid
preds <- as.data.frame(fitted(AppearMod))
plot(preds$Estimate ~ AppearMod$data$appearance)
abline(0, 1, col= 'red')

#Shrub
preds <- as.data.frame(fitted(TotalMod))
plot(preds$Estimate ~ TotalMod$data$total)
abline(0, 1, col= 'red')




####END




