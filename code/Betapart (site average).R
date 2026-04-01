##########################
#Turnover over time#
##########################
library(codyn)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(brms)
library(betapart)
library(tidyr)
library(ggplot2)
library(paletteer)
library(Rcpp)
library(bayesplot)
library(cowplot)
library(gridExtra)
library(tibble)
library(beepr)
library(bayestestR)
library(GGally)


##########################
# Step 1. Calculate turnover####
##########################


Florsitics <-read.csv("data/Raw floristic data.csv")

#Convert to long form
TO_long <- Florsitics %>% 
  group_by(Year, Site, Transect, Aspect) %>% 
  gather(key = "Species", value = "Score", 5:107) 

Florsitics <- TO_long %>% 
  #filter(Score != "0") %>% 
  group_by(Year, Site, Species) %>%
  dplyr::summarise(avg = mean(Score, na.rm = FALSE)) %>% 
  spread(key = Species, value = avg)

##### TURN INTO PA DATA IF NEEDED
Florsitics[, 3:105] <- ifelse(Florsitics[, 3:105] > 0, 1, 0)


####  beta.temp. turnover
time1 <- Florsitics %>% filter(Year == 2004)
time2 <- Florsitics %>% filter(Year == 2012)
time2_2 <- Florsitics %>% filter(Year == 2012)
time3 <- Florsitics %>% filter(Year == 2017)
time4 <- Florsitics %>% filter(Year == 2022)

#Get sites
site2004 <- Florsitics %>% filter(Year == 2004) %>% select(Site)
site12_22 <- Florsitics %>% filter(Year == 2012) %>% select(Site)

#Remove the year column
time1 <- time1[,-1]
time2 <- time2[,-1]
time2_2 <- time2_2[,-1]
time3 <- time3[,-1]
time4 <- time4[,-1]

#Remove sites from time2_2 that werent in 2004
time2_2 <- time2_2[time2_2$Site %in% time1$Site,]
unique(time2_2$Site)

#Make the the site column into the row names
rownames(time1) <- time1$Site
rownames(time2) <- time2$Site
rownames(time2_2) <- time2_2$Site
rownames(time3) <- time3$Site
rownames(time4) <- time4$Site

# Remove the 'Site' column
time1 <- time1[ , !(names(time1) %in% c("Site"))]
time2 <- time2[ , !(names(time2) %in% c("Site"))]
time2_2 <- time2_2[ , !(names(time2_2) %in% c("Site"))]
time3 <- time3[ , !(names(time3) %in% c("Site"))]
time4 <- time4[ , !(names(time4) %in% c("Site"))]

##### adespatial package #####
## looks for exceptions change in composition
site2004 
site12_22

library(adespatial)

dif0<- TBI(time1,time2_2, method = "sorensen", 
           pa.tr = TRUE, nperm = 999, BCD = TRUE, test.BC = TRUE,
           test.t.perm = TRUE)
dif1<- TBI(time2,time3, method = "sorensen", 
           pa.tr = TRUE, nperm = 999, BCD = TRUE, test.BC = TRUE,
           test.t.perm = TRUE)
dif2<- TBI(time3,time4, method = "sorensen", 
           pa.tr = TRUE, nperm = 999, BCD = TRUE, test.BC = TRUE,
           test.t.perm = TRUE)

#Look at summaries
dif0$BCD.mat 
#(D = (B + C) / den) = total dissimilarities
#species losses (B/den) and species gains (C/den)

#Get raw p-values
dif0$p.TBI
dif1$p.TBI
dif2$p.TBI

#adjusted for multiple testing 
dif0$p.adj
dif1$p.adj
dif2$p.adj


##### Step 1) Year 1 #####
#Columns B and C indicate which of the D values are 
#associated with large B (losses) or arge C values (gains)
dif0_df <- data.frame(Site = site2004$Site,
                      dif = dif0$TBI,
                      dif.p = dif0$p.TBI)

dif0_df2 <- cbind(dif0_df, dif0$BCD.mat)

#split Site into 2 columns by "_" and join with Site details
Site_details <- read.csv("data/Site evelvations.csv")

dif0_df3 <- dif0_df2 %>% 
  separate(Site, c("Site", "Transect"), sep = "_") %>% 
  left_join(Site_details, by = "Site")


##### Step 2) Year 2 #####
#Columns B and C indicate which of the D values are 
#associated with large B (losses) or arge C values (gains)
dif1_df <- data.frame(Site = site12_22$Site,
                      dif = dif1$TBI,
                      dif.p = dif1$p.TBI)

dif1_df2 <- cbind(dif1_df, dif1$BCD.mat)

#Join with Site details
dif1_df3 <- dif1_df2 %>% 
  separate(Site, c("Site", "Transect"), sep = "_") %>% 
  left_join(Site_details, by = "Site")

##### Step 3) Year 3 #####
#Columns B and C indicate which of the D values are 
#associated with large B (losses) or large C values (gains)
dif2_df <- data.frame(Site = site12_22$Site,
                      dif = dif2$TBI,
                      dif.p = dif2$p.TBI)

dif2_df2 <- cbind(dif2_df, dif2$BCD.mat)

#Join site details
dif2_df3 <- dif2_df2 %>% 
  separate(Site, c("Site", "Transect"), sep = "_") %>% 
  left_join(Site_details, by = "Site")



##### Step 4) Join al years together ######
dif0_df3$Year <- "2012"
dif1_df3$Year <- "2017"
dif2_df3$Year <- "2022"

#Join all years together
dif_all <- full_join(dif0_df3, dif1_df3)
dif_all2 <- full_join(dif_all, dif2_df3)
glimpse(dif_all2)


#### SPECIES GAINS #####
#Calculate the mean and CI for species gains
Gains_mean <- dif_all2 %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(
    mean = mean(`C/(2A+B+C)`, na.rm = TRUE),
    sd = sd(`C/(2A+B+C)`, na.rm = TRUE),
    n = sum(!is.na(`C/(2A+B+C)`)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

#Plot 
Gains_plot <- ggplot(Gains_mean, aes(x = Year, y = mean, fill = Year)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1, position = position_dodge(width = 0.5)) +
  paletteer::scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  labs(x = "Year",
       y = "Species gains",
       fill = "RAC") +
  theme_minimal() +
  theme_classic() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

Gains_plot <- Gains_plot +
  annotate("text", x = Inf, y = Inf, label = "(A)", vjust = 1, hjust = 1, size = 5)

#Save plot
ggsave("figures/Gains_plot.png", Gains_plot, width = 10, height = 6, units = "in", dpi = 300)


#### SPECIES LOSSES
#Calcuate the mean and CI for species gains
Losses_mean <- dif_all2 %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(
    mean = mean(`B/(2A+B+C)`, na.rm = TRUE),
    sd = sd(`B/(2A+B+C)`, na.rm = TRUE),
    n = sum(!is.na(`B/(2A+B+C)`)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

#Plot 
Losses_plot <- ggplot(Losses_mean, aes(x = Year, y = mean, fill = Year)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1, position = position_dodge(width = 0.5)) +
  paletteer::scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  labs(x = "Year",
       y = "Species losses") +
  theme_minimal() +
  theme_classic() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))


Losses_plot <- Losses_plot +
  annotate("text", x = Inf, y = Inf, label = "(B)", vjust = 1, hjust = 1, size = 5)

#Save plot
#Join plots together

plot <- gridExtra::grid.arrange(Gains_plot, Losses_plot,Rank_change_plot,
                                ncol = 2)
ggsave("figures/Turnover_plot3.png", plot, width = 8, height = 5, units = "in", dpi = 300)


######### Species losses and gain against elevation
Losses_ele_plot <- ggplot(data = dif_all2, aes(x = Elevation_new, y = `B/(2A+B+C)`))+
  geom_point(
    color = "grey",
    size = 2,
    alpha = 0.7)+
  geom_smooth(method="lm"
              #make line black
              , color = "black"
  )+
  ylab('Species losses') +
  xlab('Elevation (m)') +
  cowplot::theme_cowplot() + 
  ylim(0, 1) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

gains_ele_plot <- ggplot(data = dif_all2, aes(x = Elevation_new, y = `C/(2A+B+C)`))+
  geom_point(
    color = "grey",
    size = 2,
    alpha = 0.7)+
  geom_smooth(method="lm"
              #make line black
              , color = "black"
  )+
  ylab('Species gains') +
  xlab('Elevation (m)') +
  cowplot::theme_cowplot() + 
  ylim(0, 1) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

RC_ele_plot <- ggplot(data = RAC_rank_change, aes(x = Elevation_new, y = rank_change))+
  geom_point(
    color = "grey",
    size = 2,
    alpha = 0.7)+
  geom_smooth(method="lm"
              #make line black
              , color = "black"
  )+
  ylab('Rank change') +
  xlab('Elevation (m)') +
  cowplot::theme_cowplot() + 
  ylim(0, 0.5) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

gains_ele_plot <- gains_ele_plot +
  annotate("text", x = Inf, y = Inf, label = "(A)", vjust = 1, hjust = 1, size = 5)
Losses_ele_plot <- Losses_ele_plot +
  annotate("text", x = Inf, y = Inf, label = "(B)", vjust = 1, hjust = 1, size = 5)
RC_ele_plot <- RC_ele_plot +
  annotate("text", x = Inf, y = Inf, label = "(C)", vjust = 1, hjust = 1, size = 5)



#Save plot
#Join plots together
plot <- gridExtra::grid.arrange(gains_ele_plot, Losses_ele_plot,RC_ele_plot,
                                ncol = 1)
ggsave("figures/Turnover_elev_plot2.png", plot, width = 5, height = 8, units = "in", dpi = 300)



#########
#######
# Modeling dissimilarity #######

#### Step 1) Prepare dataset for modelling  #####
#D = (B + C) / den)

dif_all2

#Plot 
Dis_plot <- ggplot(dif_all2, aes(x = Year, y = `D=(B+C)/(2A+B+C)`, fill = Year)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.5)) +
  labs(x = "Year",
       y = "Dissimilarity") +
  theme_minimal() +
  theme_classic() +
  ylim(0, 1) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))

#Dis_plot <- Dis_plot +
 # annotate("text", x = Inf, y = Inf, label = "(A)", vjust = 1, hjust = 1, size = 5)
#Save plot
ggsave("figures/Gains_plot.png", Dis_plot, width = 10, height = 6, units = "in", dpi = 300)


#Join Fire data (all sites), Site evelvations and Summit CGDD csv
fire <- read.csv("data/Time since burnt.csv")
Elev <- read.csv("data/Site evelvations.csv")
CGDD <- read.csv("data/Summit CGDD.csv")
#remove 2004
fire <- fire[fire$Year != 2004,]
CGDD <- CGDD[CGDD$Year != 2004,]


#join fire, site elevations and CGDD data with dif_all2
dif_all3 <- merge(dif_all2, fire, by = c("Site", "Year"))
dif_all3 <- merge(dif_all3, Elev, by = "Site")
dif_all3 <- merge(dif_all3, CGDD, by = c("Site", "Year"))

#check
glimpse(dif_all3)

#Rename columns `B/(2A+B+C)` = Losses, `C/(2A+B+C)` = Gains, `D=(B+C)/(2A+B+C)` = Dissimilarity
colnames(dif_all3)[colnames(dif_all3) == "B/(2A+B+C)"] <- "Losses"
colnames(dif_all3)[colnames(dif_all3) == "C/(2A+B+C)"] <- "Gains"
colnames(dif_all3)[colnames(dif_all3) == "D=(B+C)/(2A+B+C)"] <- "Dissimilarity"
#rename columns Elevation_new.x
colnames(dif_all3)[colnames(dif_all3) == "Elevation_new.x"] <- "Elevation"
#remove Elevation_new.y 
dif_all3 <- dif_all3[, !grepl("Elevation_new.y", colnames(dif_all3))]

#check
glimpse(dif_all3)

##
###
## Step 2) GDM ####
###
##

glimpse(dif_all3)

### Step 1: Check for correlation between variables (<0.7) #######
dif_all3 %>%
  dplyr::select(-c(Site, Year, Transect, dif, dif.p, Change, 
                   Time_since_grazing,Grazing.ceased)) %>%
  GGally::ggpairs() 

### Step 2: Standardise variables & look at distribution ######
dif_all3$Year.F= as.factor(dif_all3$Year)
dif_all3$TSF.std= as.vector(scale(dif_all3$Time_since_burnt)) 
dif_all3$gdd.std = as.vector(scale(dif_all3$cum.gdd))
dif_all3$ele.std = as.vector(scale(dif_all3$Elevation))

#What is the distribution of the response variable
hist(dif_all3$Dissimilarity)

### Step 3: Create model ####
DisGLM <- brms::brm(Dissimilarity ~ 
                        TSF.std*gdd.std + ele.std +
                        (1 | Site) + (1|Year.F),
                        family = Beta(link_phi = "log"),
                      iter = 4000,  
                      seed = 123,
                      cores = 4,
                      chains = 4,
                      control = list(adapt_delta = 0.999, max_treedepth=15), 
                      data = dif_all3)

### Step 5: Validate model ####
#Check model
summary(DisGLM) #summary stats - look at Rhat and make sure <1.05
plot(DisGLM) #check Convergence of the model by looking at the Markov chains
pp_check(DisGLM) #check model predictions to observations
bayes_R2(DisGLM) # check the R2

### Step 6: Fixed effects coefficient plot ####
Dis_fixed <- as.data.frame(fixef(DisGLM))
Dis_fixed2 <- rownames_to_column(Dis_fixed, var = "Variable")
Dis_fixed2$Variable <- factor(Dis_fixed2$Variable,
                                 levels = c("TSF.std:gdd.std",
                                            "gdd.std", "TSF.std",
                                            "ele.std", "Intercept"))

custom_labels <- c( "ele.std" = "Elevation",
                    "gdd.std" = "CGD",
                    "TSF.std:gdd.std" = "TSF:CGD",
                    "TSF.std" = "TSF",
                    "Intercept" = "Intercept")

#Get the 80% credible intervals from model
ci_eti_fixed_89 <- ci(DisGLM, method = "ETI", ci = 0.89, effects = "fixed")
ci_eti_fixed_95 <- ci(DisGLM, method = "ETI", ci = 0.95, effects = "fixed")

#Rename CI low and CI high columns 
ci_eti_fixed_89 <- ci_eti_fixed_89 %>%
  rename(Q5.5 = CI_low, Q94.5 = CI_high) #Calcuates the 89% CI
ci_eti_fixed_95 <- ci_eti_fixed_95 %>%
  rename(Q2.5 = CI_low, Q97.5 = CI_high) #Calculates the 95% CI

#join the 80% and 95% credible intervals
Dis_fixed_ef <- left_join(ci_eti_fixed_89, ci_eti_fixed_95, by = "Parameter")

# add in estimate column 
Dis_fixed_ef$Estimate <- Dis_fixed2$Estimate
Dis_fixed_ef$Variable <- Dis_fixed2$Variable


#Plot
Dis_fixed <- ggplot(Dis_fixed_ef, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q5.5, xmax = Q94.5, height = 0), linewidth = 0.7) +  # 89% Confidence intervals)
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0), linewidth = 0.3) +  # 95% Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

Dis_fixed <- Dis_fixed +
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)

Dis_fixed

### Step 7: Random effects ####

# 1. Plot estimates of the standard deviation and CI of random effects 
#to see which one explains the most in the unexplained variation 
summary(DisGLM) #pull out the estimate and CI's from the model summary and create a new DF
dis_ran_DF <- data.frame(
  variable = c("Site", "Year"),
  estimate = c(0.36, 0.54),
  lci = c(0.11, 0.01),
  uci = c(0.69, 2.41)
)

Dis_rand_all <- ggplot(dis_ran_DF, aes(x = estimate, y = variable)) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = lci, xmax = uci, height = 0), linewidth = 0.3) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")

Dis_rand_all <- Dis_rand_all +
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)

ggsave("output/Forb/Forb random effects plot.png", 
       plot = Forb_rand_all, width = 7, height = 3.5) 

#### Random effects: site and year (appendix) ######## HAVENT DONE!!!
#SITE
Forb_rand <- as.data.frame(ranef(mForbGLM, groups = "Site"))
Forb_rand <- rownames_to_column(Forb_rand, var = "Variable")
Forb_rand <- Forb_rand[order(Forb_rand$Site.Estimate.Intercept), ]

Forb_rand_site <- ggplot(Forb_rand, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")


#YEAR
Forb_rand2 <- as.data.frame(ranef(mForbGLM, groups = "Year.F"))
Forb_rand2 <- rownames_to_column(Forb_rand2, var = "Variable")
Forb_rand2 <- Forb_rand2[order(Forb_rand2$Year.F.Estimate.Intercept), ]

Forb_rand_year <- ggplot(Forb_rand2, aes(x = Year.F.Estimate.Intercept, y = reorder(Variable, Year.F.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Year.F.Q2.5.Intercept, xmax = Year.F.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")


#Save both plots
plot <- grid.arrange(Forb_rand_site, Forb_rand_year,
                     ncol = 2)
#ggsave("output/Forb/Forb model output (random effects).png", 
       #plot = plot, width = 10, height = 3.5) 



### Step 8: Plot marginal predictions ####

#### Step 8.1: Time since fire X CGD ####

#Fire values - high and low
mean(dif_all3$Time_since_burnt) #21.94444
sd(dif_all3$Time_since_burnt) #25.81577
range(dif_all3$Time_since_burnt) #4 83
##LOW fire frequency (at 83 years unburnt - transform (83 -mean)/SD = 2.365049)
(83 - 21.94444)/25.81577 
#HIGH fire frequency (at 4 years unburnt - transform (4 -mean)/SD = -0.6950961)
(4 -21.94444)/25.81577 

#CGD values
mean(dif_all3$cum.gdd) #57109.19
sd(dif_all3$cum.gdd) #12562.7
range(dif_all3$cum.gdd) #35015.9 77193.0

### LOW FIRE PREDICTIONS = 83 years##
#Create new dataframe
newdat_DIS_F_CGD_Low <- data.frame(Site = "average site",
                                    Year.F = "average year",
                                    cum.gdd = seq(35000, 77200, by = 100), 
                                    TSF.std = 2.365049, #Predicting at low fire freq
                                    ele.std = 0) %>% 
  dplyr::mutate(gdd.std = (cum.gdd - 57109.19)/12562.7)
#Predictions
predictions<- posterior_epred(DisGLM, newdat= newdat_DIS_F_CGD_Low, re_formula = NULL, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_DIS_F_CGD_Low$mean <- mean
newdat_DIS_F_CGD_Low$lci <- lci
newdat_DIS_F_CGD_Low$uci <- uci


#### HIGH FIRE FREG PREDICTIONS = 4 years ##
newdat_DIS_F_CGD_High <- data.frame(Site = "average site",
                                     Year.F = "average year",
                                     cum.gdd = seq(35000, 77200, by = 100), 
                                     TSF.std = -0.6950961, #Predicting at high fire freq
                                     ele.std = 0) %>% 
  dplyr::mutate(gdd.std = (cum.gdd - 57109.19)/12562.7)
#Predictions
predictions<- posterior_epred(DisGLM, newdat= newdat_DIS_F_CGD_High, re_formula = NULL, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_DIS_F_CGD_High$mean <- mean
newdat_DIS_F_CGD_High$lci <- lci
newdat_DIS_F_CGD_High$uci <- uci


##Plot predictions
Forb2_plot_Fire_CGD <- ggplot(newdat_DIS_F_CGD_High, aes(x = cum.gdd, y = mean)) + 
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_DIS_F_CGD_Low, aes(x = cum.gdd, y = mean), color = "blue") +
  geom_ribbon(data = newdat_DIS_F_CGD_Low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  #Put in data points 
  #geom_point(data = Forbs2, aes(x = Time_since_burnt, y = Cover), 
  #shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  theme_cowplot()+
  ylim(0,1) +
  ylab(expression("Dissimilarity"))+
  xlab(expression("Cumulative summer growing degrees (°C)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

Forb2_plot_Fire_CGD <- Forb2_plot_Fire_CGD +
  annotate("text", x = Inf, y = Inf, label = "(c)", vjust = 1, hjust = 1, size = 6)

Forb2_plot_Fire_CGD



#### Step 8.2: Elevation ####

mean(Forbs2$ele) #1752.115
sd(Forbs2$ele) #96.6135
range(Forbs2$ele) #1574.000 1982.326
#Create new dataframe
newdat_forb4 <- data.frame(Site = "average site",
                           Year.F = "average year",
                           ele = seq(1500, 2000, by = 10), 
                           gdd.std = 0,
                           TSF.std = 0, 
                           total = 1) %>% 
  dplyr::mutate(ele.std = (ele - 1752.115)/96.6135)
#Predictions
predictions<- posterior_epred(mForbGLM, newdat= newdat_forb4, allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_forb4$mean <- mean
newdat_forb4$lci <- lci
newdat_forb4$uci <- uci


##Plot predictions
Forb2_plot_ele <- ggplot(newdat_forb4, aes(x = ele, y = mean)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  #geom_point(data = Forbs2, aes(x = ele, y = Cover), 
  #shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  #geom_line(size = 0.75) +
  theme_cowplot()+
  ylim(0,1) +
  ylab(expression(" "))+
  xlab(expression("Elevation (m)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

Forb2_plot_ele <- Forb2_plot_ele+
  annotate("text", x = Inf, y = Inf, label = "(d)", vjust = 1, hjust = 1, size = 6)

Forb2_plot_ele


#### Step 8.3: Plot all together #####
plot_forb2 <- grid.arrange(Forb2_fixed, Forb_rand_all,
                           Forb2_plot_Fire_CGD,
                           Forb2_plot_ele,
                           ncol = 2)
#Save
ggsave("output/Forb/Forb GLM.png", 
       plot = plot_forb2, width = 11, height = 9)



