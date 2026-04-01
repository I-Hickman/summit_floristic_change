Temp_change_from_baseline <- read.csv("Data/Temp_change_from_baseline.csv")

temp <- Temp_change_from_baseline %>% 
  gather(key = "Elevation", value = "Temp", 2:9) 


#### Step 2: MODEL (I have removed 'Time since grazing')
mTemp <- brms::brm(Temp ~ 
                         + Year +
                         ( 1 | Elevation),
                       family = gaussian(),
                       iter = 2000, #increase if model not converging 
                       cores = 4,
                       chains = 4,
                       control = list(adapt_delta = 0.8, max_treedepth=15), # if you get a transition warning message - increase adapt_delta to 0.9 (was 0.8)
                       data = temp)


# 1. Change in max summer temp

range(temp$Year) #2000 2022
#Create new dataframe
newdat <- data.frame(Elevation = "average elevation",
                     Year = seq(2000, 2022, by = 10)
#Predictions
predictions<- posterior_epred(mTemp, newdat= newdat, re.form = NA)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
newdat_Shrub3$mean <- mean
newdat_Shrub3$lci <- lci
newdat_Shrub3$uci <- uci


##Plot predictions
temp <- ggplot(temp, aes(x = Year, y = Temp)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_point(data = temp, aes(x = Average_change_max_summer_temp, y = Cover), 
             shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  #geom_line(size = 0.75) +
  theme_cowplot()+
  ylim(0,0.8) +
  ylab(expression("Shrub proportional cover"))+
  xlab(expression("Δ in summer max temp (°)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 