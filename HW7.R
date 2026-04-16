#read in data
ghg <- read.csv("/cloud/project/activity07/Deemer_GHG_Data.csv")
ETdata <- read.csv("/cloud/project/activity07/ETdata.csv")
#load packages 
library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)

# log transform methane fluxes
ghg$log.ch4 <- log(ghg$ch4+1)
#transform other metrics 
ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)
ghg$log.surfArea <- log(ghg$surface.area)


unique(ghg$Region)

# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)
# binary variable for alpine region
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)
# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)

#original model 
mod.full <- lm(log.ch4 ~ airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV, data=ghg)
#check assumptions
res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

# shapiro-wilks test
shapiro.test(res.full)

#residuals
plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)

# isolate continuous model variables into data frame:

reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)
# run stepwise
full.step <- ols_step_forward_aic(mod.full)
# view table
full.step 
# check full model
full.step$model

# plot AIC over time
plot(full.step )

#create improved model
hydro_mod.full <- lm(log.ch4 ~ airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV +log.surfArea, data=ghg)

plot(hydro_mod.full)
summary(hydro_mod.full)

#check assumptions 
hydro_res.full <- rstandard(hydro_mod.full)
hydro_fit.full <- fitted.values(hydro_mod.full)

# check normality with qq plot
qqnorm(hydro_res.full, pch=19, col="grey50")
qqline(hydro_res.full)

# shapiro-wilks test
shapiro.test(hydro_res.full)

#residuals 
plot(hydro_fit.full, hydro_res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:

hydro_reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip, 
                       ghg$log.surfArea)
# make a correlation matrix 
chart.Correlation(hydro_reg.data, histogram=TRUE, pch=19)

# run stepwise
hydro_full.step <- ols_step_forward_aic(hydro_mod.full)
# view table
hydro_full.step 

#check full model
hydro_full.step$model

# plot AIC over time
plot(hydro_full.step )

#count NAs
sum(is.na(ghg$chlorophyll.a))
sum(is.na(ghg$hydropower))
sum(is.na(ghg$runoff))

summary(hydro_mod.full)
summary(mod.full)

table_out<- summary(hydro_mod.full)$coefficents

plot(hydro_mod.full)

write.csv(table_out, "/cloud/project/ch4_mod.csv", row.names = TRUE)

# Homework Questions ------------------------------------------------------
ETdat <- read.csv("/cloud/project/activity07/ETdata.csv")
ETdat

unique(ETdat$crop)

library(lubridate)
library(ggplot2)
library(forecast)
library(dplyr)

#Setting up a Time Series 
#characterize almond water use in the region across multiple fields

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

#set up regular intervals 
# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)

# plot decomposition
plot(almond_dec)

almondTrend <- almond_dec$trend
almondSeason <- almond_dec$seasonal

#autocorrelation - measure of correlation from one observation 
#to the observations before and after
#look at lag periods
acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

#autoregressive models
pacf.plot <- pacf(na.omit(almond_ts))

#omit NAs
almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

#try running 4th order
model4 <- arima(almond_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4

# calculate fit
AR_fit1 <- almond_y - residuals(model1) 
AR_fit4 <- almond_y - residuals(model4)
#plot data
plot(almond_y)
# plot fit
points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

#forecast future data
newAlmond <- forecast(model4)
newAlmond

#make dataframe for plotting
newAlmondF <- data.frame(newAlmond)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newAlmondF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = almond, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(almond$date[1]),newAlmondF$dateF[24])+  # Plotting original data
  geom_line(data = newAlmondF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newAlmondF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")
