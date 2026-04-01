####
####
####
# Step 4: Community stability and variance ratio #####
####
####
####


library(codyn)
library(tidyverse)
library(ggplot2)
library(dplyr)

Florsitics <-read.csv("data/Raw floristic data.csv")
Site_F <- read.csv("data/Fire history (broad).csv")

#Convert to long form
TO_long <- Florsitics %>% 
  group_by(Year, Site, Transect, Aspect) %>% 
  gather(key = "Species", value = "Score", 5:107) 

#First need to get the average abundance of species across aspects
sp_abundnace_site <- TO_long %>% 
  #filter(Score != "0") %>% 
  group_by(Year, Site, Species) %>%
  dplyr::summarise(avg = mean(Score, na.rm = FALSE))



########Community stability and variance ratio ########

plan(multisession, workers = 4)
variance.ratio <- variance_ratio(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site",
  bootnumber = 10000, #recommended null model repeat 
  average.replicates = FALSE,
  level = 0.95
)


Community.stability <- community_stability(
  df = sp_abundnace_site,
  time.var = "Year",
  abundance.var = "avg",
  replicate.var = "Site"
)

#Join Variance ratio and community stability
Community_stab_metrics <- Community.stability %>% 
  left_join(variance.ratio, by = "Site")
#Add in fire history
Community_stab_metrics <- left_join(Community_stab_metrics, Site_F, by = "Site")
#Plot
ggplot(data = Community_stab_metrics, aes(x =VR, y = stability)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ylab("Community stability") +
  xlab("Variance ratio")


##### Plot with loess function ####

## With loess function
cor_plot1 <- Community_stab_metrics %>%
  ggplot(aes(x = VR, y = stability, colour = Fire.history)) +
  geom_point(size = 2, alpha = 1) +
  #labels for each point
  stat_smooth(method = loess, 
              colour = "black", 
              na.rm = TRUE,
              linewidth = 0.9,
              se = TRUE) +
  #ggpubr::stat_cor(#aes(label = after_stat(r.label)), 
   #                  method = "spearman",
    #                 label.y = 10, label.x = 2.5,
     #                colour = "black",
      #               cor.coef.name = "rho") +
  theme_classic() +
  ylim(0,10) +
  labs(x = "Variance ratio", 
       y = "Community stability",
       color = "Fire history") +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) +
  paletteer::scale_fill_paletteer_d("nationalparkcolors::Acadia")
cor_plot1


##### Synchrony ####
Synchrony_L <- codyn::synchrony(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site",
  metric = "Loreau"
)

Synchrony_G <- codyn::synchrony(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site",
  metric = "Gross"
)

#Join community stability metrics to synchony
synchrony_metrics_G2 <- left_join(Synchrony_G, Community_stab_metrics, by = "Site")
synchrony_metrics_L2 <- left_join(Synchrony_L, Community_stab_metrics, by = "Site")

#plot
ggplot(data = synchrony_metrics_G2, aes(x =synchrony, y = stability)) + 
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm") +
  ylab("Community stability") +
  xlab("Synchrony (Gross)") +
  ylim(0,10)

ggplot(data = synchrony_metrics_G2, aes(x =synchrony, y = stability)) + 
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm") +
  ylab("Community stability") +
  xlab("Synchrony (Loreau)") +
  ylim(0,10)

##### Correlation
## With loess function
cor_plot2 <- synchrony_metrics_G2 %>%
  ggplot(aes(x = synchrony, y = stability, colour = Fire.history)) +
  geom_point(size = 2, alpha = 1) +
  #labels for each point
  stat_smooth(method = loess, 
              colour = "black", 
              na.rm = TRUE,
              linewidth = 0.9,
              se = TRUE) +
  ggpubr::stat_cor(#aes(label = after_stat(r.label)), 
                    method = "spearman",
                   label.y = 10, label.x = 0.5,
                  colour = "black",
                 cor.coef.name = "rho") +
  theme_classic() +
  ylim(0,10) +
  labs(x = "Synchrony (Gross)", 
       y = "Community stability",
       color = "Fire history") +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) 
cor_plot2







############### cyclic_shift() function ####

#Significance testing for community stability

cyclic_shift()

cyclic_shift <- cyclic_shift(
  df = sp_abundnace_site,
  time.var = "Year",
  species.var = "Species",
  abundance.var = "avg",
  replicate.var = "Site", 
  FUN = cov,
  bootnumber = 10
)

confint.cyclic_shift() #returns confidence intercals for cyclic_shift




#### MODEL ######

#Plot with R2
cor_plot <- wind_data3 %>%
  ggplot(aes(x = Tvel, y = perc_threshold_events, colour = DS)) +
  geom_point(size = 2, alpha = 1) +
  stat_smooth(method = lm, colour = "black", linewidth = 0.9) + # se = FALSE) +
  stat_cor(aes(label = after_stat(r.label)), 
           method = "spearman",
           label.y = 0.4, label.x = 4,
           colour = "#d95f02",
           cor.coef.name = "rho") +
  stat_cor(aes(label = after_stat(p.label)), 
           method = "spearman",
           label.y = 0.35, label.x = 4,
           colour = "#d95f02") +
  theme(axis.text.y = element_blank()) +
  theme_cowplot() +
  #ylim(-0.03,0.4) +
  labs(x = "Terminal velocity (m/s)", 
       y = "Proportion of vertical wind events > threshold",
       color = "Dispersal syndrome") +
  theme_cowplot()+
  ylab(expression(atop("Proportion of vertical wind", paste("events above threshold"))))+ 
  xlab(expression("Terminal velocity (m/s)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) 
cor_plot

#Plot without correlation 
cor_plot2 <- wind_data3 %>%
  ggplot(aes(x = Tvel, y = perc_threshold_events, colour = DS)) +
  geom_point(size = 2, alpha = 1) +
  stat_smooth(method = gam, 
              colour = "black", 
              na.rm = TRUE,
              linewidth = 0.9, 
              se = FALSE) +
  theme(axis.text.y = element_blank()) +
  theme_cowplot() +
  ylim(0,0.4) +
  labs(x = "Terminal velocity (m/s)", 
       y = "Proportion of vertical wind events > threshold",
       color = "Dispersal syndrome") +
  theme_cowplot()+
  ylab(expression(atop("Proportion of vertical wind", paste("events above threshold"))))+ 
  xlab(expression("Terminal velocity (m/s)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) 
cor_plot2

## With loess function
cor_plot3 <- wind_data3 %>%
  ggplot(aes(x = Tvel, y = perc_threshold_events, colour = DS)) +
  geom_point(size = 2, alpha = 1) +
  stat_smooth(method = loess, 
              colour = "black", 
              na.rm = TRUE,
              linewidth = 0.9,
              se = FALSE) +
  theme(axis.text.y = element_blank()) +
  theme_cowplot() +
  ylim(0,0.4) +
  labs(x = "Terminal velocity (m/s)", 
       y = "Proportion of vertical wind events > threshold",
       color = "Dispersal syndrome") +
  theme_cowplot()+
  ylab(expression(atop("Proportion of vertical wind", paste("events above threshold"))))+ 
  xlab(expression("Terminal velocity (m/s)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) 
cor_plot3

#Save
ggsave("output/Vertical wind events model/Vertical wind correlation plot (without numbers).png", 
       plot = cor_plot2, width = 8, height = 3.5) 