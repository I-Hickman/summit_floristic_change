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


##########################
# Step 1. Calculate turnover####
##########################

Florsitics <-read.csv("data/Raw floristic data.csv")

####  beta.temp. turnover
#Merge the site and transect together with a "_" in between
Florsitics$Site <- paste(Florsitics$Site, Florsitics$Transect, sep = "_")

#Select year = 2004
Florsitics <- Florsitics[,-c(3,4)]

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

#Calcuate between time 1 and time 2  
turnover1 <-  beta.temp(time1, time2, index.family="sorensen")
#Calcuate between time 2 and time 3
turnover2 <-  beta.temp(time2, time3, index.family="sorensen")
#Calcuate between time 3 and time 4
turnover3 <-  beta.temp(time3, time4, index.family="sorensen")


##### adespatial package #####
## looks for exceptions change in composition
site2004 
site12_22

library(adespatial)

dif0<- TBI(time1,time2_2, method = "%difference", 
           pa.tr = TRUE, nperm = 999, BCD = TRUE, test.BC = TRUE,
           test.t.perm = TRUE)
dif1<- TBI(time2,time3, method = "%difference", 
            pa.tr = TRUE, nperm = 999, BCD = TRUE, test.BC = TRUE,
            test.t.perm = TRUE)
dif2<- TBI(time3,time4, method = "%difference", 
           pa.tr = TRUE, nperm = 999, BCD = TRUE, test.BC = TRUE,
           test.t.perm = TRUE)

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
#Calcuate the mean and CI for species gains
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
plot <- gridExtra::grid.arrange(Gains_plot, Losses_plot,
                                ncol = 1)
ggsave("figures/Turnover_plot2.png", plot, width = 5, height = 7, units = "in", dpi = 300)



##### Species losses and gain against elevation
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

####### Rank change 
RAC2 <- left_join(RAC, Site_evelvations)
RAC_rank_change <- RAC2[, c("Site", "Elevation_new", "Year2", "rank_change")]

Rank_mean <- RAC2 %>%
  ungroup() %>%
  group_by(Year2) %>%
  summarise(
    mean = mean(rank_change, na.rm = TRUE),
    sd = sd(rank_change, na.rm = TRUE),
    n = sum(!is.na(rank_change)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

#Plot 
Rank_change_plot <- ggplot(Rank_mean, aes(x = Year2, y = mean, fill = Year2)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1, position = position_dodge(width = 0.5)) +
  paletteer::scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  labs(x = "Year",
       y = "Rank change") +
  theme_minimal() +
  theme_classic() +
  ylim(0, 0.5) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14))


Rank_change_plot <- Rank_change_plot +
  annotate("text", x = Inf, y = Inf, label = "(C)", vjust = 1, hjust = 1, size = 5)

#Save plot
#Join plots together
plot <- gridExtra::grid.arrange(Gains_plot, Losses_plot,Rank_change_plot,
                                ncol = 2)
ggsave("figures/Turnover_plot3.png", plot, width = 8, height = 5, units = "in", dpi = 300)




#Rank change vs elevation
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




### Average turnover for each site
Losses_mean_site <- dif_all2 %>%
  ungroup() %>%
  group_by(Site, Elevation_new) %>%
  summarise(
    mean = mean(`B/(2A+B+C)`, na.rm = TRUE),
    sd = sd(`B/(2A+B+C)`, na.rm = TRUE),
    n = sum(!is.na(`B/(2A+B+C)`)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()

Gains_mean_site <- dif_all2 %>%
  ungroup() %>%
  group_by(Site, Elevation_new) %>%
  summarise(
    mean = mean(`C/(2A+B+C)`, na.rm = TRUE),
    sd = sd(`C/(2A+B+C)`, na.rm = TRUE),
    n = sum(!is.na(`C/(2A+B+C)`)),
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) %>%
  ungroup()




##### Test what processs is driving community change
Change_df <-RAC_df %>% 
  select(Site, Year2, composition_change)

Change_df <- left_join(Change_df, Site_details, by = "Site")




