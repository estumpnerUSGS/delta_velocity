#retrieval, imputation, and summary of instantaneous USGS velocity data
#USGS 11455350 CACHE SLOUGH A RYER ISLAND
#USGS 11455385 CACHE SLOUGH AB RYER ISLAND FERRY NR RIO VISTA CA

library(dataRetrieval)
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(imputeTS)
library(zoo)

#data retrieval pull through NWIS web services

siteNumbers <- c("11455350")  #2007-10-01 	2019-04-01 ,# 2018-07-03 	2022-07-17

startDate <- ""

endDate <- ""

parameterCd <- "72255"

uvRYI <- readNWISuv(siteNumbers,parameterCd,startDate,endDate)

#convert to PST

hrs <- 8 * 60 * 60

uvRYI$dateTime_PST <- as.POSIXct(uvRYI$dateTime, format = "%Y-%m-%d %H:%M:%S") - hrs

uvRYI <- subset(uvRYI, select = c(4, 7))

uvRYI <- uvRYI %>%
  rename(velocity_ft_s = "X_72255_00000")

#write out to .csv

write_csv(uvRYI, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvRYI_vel.csv")

#import .csv

uvRYI <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvRYI_vel.csv")

Tgodinfn<-function(tint,xint) {

  # apply Godin filter
  # it includes three passes in time domain, with windows 24, 24, and 25 hrs long
  # will assume 15-min data here, so
  # for first pass, 48 points have to be skipped and last 44
  # window is 12-1-11 hrs in length, i.e. use previous 12 hours and
  # subsequent 11 hrs to define result at 13th hour.

  xfilt<-rep_len(NA,length.out=length(xint))

  # do first pass, with window of 12,1,11 hrs
  for (i in seq(from=49,to=(length(xint)-44))) {
    xfilt[i]=sum(xint[(i-48):(i+44)])/93
  }

  xfilt1<-xfilt
  # now 2nd pass, same approach but window is 11,1,12 hrs
  for (i in seq(from=45, to=(length(xint)-48))) {
    xfilt1[i]=sum(xfilt[(i-44):(i+48)])/93;
  }
  xfilt=xfilt1;
  # now 3rd pass with 12-1-12 window
  for (i in seq(from=49,to=(length(xint)-48))) {
    xfilt[i]=sum(xfilt1[(i-48):(i+48)])/97;
  }
  # Now throw away 36 hrs at each end
  xfilt[1:(36*4)]<-NA
  isize=length(xfilt)
  xfilt[(isize-36*4):isize]<-NA

  return(xfilt)

}

#apply godin filter------------------------
uvRYI$vel_tf <- Tgodinfn(uvRYI$dateTime_PST, uvRYI$velocity_ft_s)

#tidal flow = velocity - net velocity(tidally-filtered)

uvRYI$tidal_vel <- uvRYI$velocity_ft_s - uvRYI$vel_tf

#plots - take a minute to load
#plot(uvRYI$dateTime_PST, uvRYI$velocity_ft_s)
#lines(uvRYI$dateTime_PST, uvRYI$tidal_vel, col = 3)
#lines(uvRYI$dateTime_PST, uvRYI$vel_tf, col = 2)
#legend("topright", legend = c("velocity", "tidal_vel", "net_vel"), col = c(1:3), lwd = 2)

#calculate amplitude- first break df in two
uvRYI_a <- subset(uvRYI, select = c(1, 2, 3))
uvRYI_b <- subset(uvRYI, select = c(2, 4))

#calcualte 30 hr moving window
uvRYI_b$roll_max <- rollapply(uvRYI_b$tidal_vel, 120, max, by = 1, partial = TRUE)
uvRYI_b$roll_min <- rollapply(uvRYI_b$tidal_vel, 120, min, by = 1, partial = TRUE)

#add high and low amplitude and divide by two
uvRYI_b$amp <- (uvRYI_b$roll_max + abs(uvRYI_b$roll_min))/2

#plots - take a minute to load
#plot(uvRYI_b$dateTime_PST, uvRYI_b$tidal_vel, col = 1, ylim = c(-5, 5), xlim = c(as.POSIXct('2020-10-01 17:30:00', format="%Y-%m-%d %H:%M:%S"), as.POSIXct('2020-10-15 17:30:00', format="%Y-%m-%d %H:%M:%S")))
#lines(uvRYI_b$dateTime_PST, uvRYI_b$roll_max, col = 2)
#lines(uvRYI_b$dateTime_PST, uvRYI_b$roll_min, col = 3)
#lines(uvRYI_b$dateTime_PST, uvRYI_b$amp, col = 5)
#legend("top", inset = c(-0.45, 0), legend = c("tidal_vel", "env_max", "env_min", "amp"), col = c(1:3, 5), lwd = 2)

#merge two dfs back together

uvRYI_c <- merge(uvRYI_a, uvRYI_b, by = 'dateTime_PST')

#subset

uvRYI_c <- subset(uvRYI_c, select = c(1:4, 7))

#calculate ratio of mean amplitude: mean tidal velocity

uvRYI_c$ratio <- uvRYI_c$amp/uvRYI_c$vel_tf

#>1 tidal influence, <1 river influence

#plot(uvRYI_c$dateTime_PST, uvRYI_c$ratio, ylim = c(-1000, 1000))

write_csv(uvRYI_c, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYI_vel.csv")

#downstep to daily - add n column

dvRYI <- uvRYI_c %>%
  mutate(date = as.Date(dateTime_PST)) %>%
  group_by(date) %>%
  summarise(mean_ratio = mean(ratio), mean_tidal_vel = mean(tidal_vel), mean_net_vel = mean(vel_tf), max_abs_velocity = max(abs(velocity_ft_s), na.rm = TRUE), min_vel_ft_s = min(velocity_ft_s, na.rm = TRUE), max_vel_ft_s = max(velocity_ft_s, na.rm = TRUE), mean_vel_ft_s = mean(velocity_ft_s, na.rm = TRUE), mean_amp = mean(amp), n_vel = n())

continous.dates <- data.frame (x = 1:4201, date = seq(as.Date('2007-10-01'),as.Date('2019-04-01'), by='day'))

RYI_vel_daily <- merge(continous.dates, dvRYI, by = "date", all = TRUE)

plot(RYI_vel_daily$date, RYI_vel_daily$n_vel)

plot(RYI_vel_daily$date, RYI_vel_daily$mean_vel_ft_s)

#if n_flow is NA replace with 0

RYI_vel_daily$n_vel[is.na(RYI_vel_daily$n_vel)] <- 0

sum(RYI_vel_daily$n_vel<=91)#78 days with <95% of flow measurements

#add column to identify if flow data is measured or will be imputed

RYI_vel_daily$group <- ifelse(RYI_vel_daily$n_vel>= 91, "measure", "impute")

#if n_value is <91, change mean, max, and min velocity to NA

RYI_vel_daily$mean_final_vel <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$mean_vel_ft_s, NA)

RYI_vel_daily$max_final_vel <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$max_vel_ft_s, NA)

RYI_vel_daily$min_final_vel <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$min_vel_ft_s, NA)

RYI_vel_daily$max_abs_final_vel <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$max_abs_velocity, NA)

RYI_vel_daily$final_tidal_vel <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$mean_tidal_vel, NA)

RYI_vel_daily$final_net_vel <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$mean_net_vel, NA)

RYI_vel_daily$final_ratio <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$mean_ratio, NA)

RYI_vel_daily$final_amp <- ifelse(RYI_vel_daily$n_vel>= 91, RYI_vel_daily$mean_amp, NA)

#impute missing values

RYI_vel_daily$final_mean_vel_final <- na_ma(RYI_vel_daily$mean_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYI_vel_daily$final_max_vel_final <- na_ma(RYI_vel_daily$max_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYI_vel_daily$final_min_vel_final <- na_ma(RYI_vel_daily$min_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYI_vel_daily$final_max_abs_final <- na_ma(RYI_vel_daily$max_abs_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYI_vel_daily$final_ratio_final <- na_ma(RYI_vel_daily$final_ratio, k = 7, weighting = "exponential", maxgap = Inf)

RYI_vel_daily$final_tidal_vel_final <- na_ma(RYI_vel_daily$final_tidal_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYI_vel_daily$final_net_vel_final <- na_ma(RYI_vel_daily$final_net_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYI_vel_daily$final_amp_final <- na_ma(RYI_vel_daily$final_amp, k = 7, weighting = "exponential", maxgap = Inf)

summary(RYI_vel_daily)

plot(RYI_vel_daily$date, RYI_vel_daily$final_mean_vel_final, ylim = c(-4, 6))
lines(RYI_vel_daily$date, RYI_vel_daily$final_max_vel_final, col = 2)
lines(RYI_vel_daily$date, RYI_vel_daily$final_min_vel_final, col = 3)
lines(RYI_vel_daily$date, RYI_vel_daily$final_max_abs_final, col = 4)

RYI_vel_daily <- RYI_vel_daily %>%
  rename(mean_vel = final_mean_vel_final,
         max_vel = final_max_vel_final,
         min_vel = final_min_vel_final,
         max_abs_vel = final_max_abs_final,
         mean_tide_vel = final_tidal_vel_final,
         ratio_mean = final_ratio_final,
         net_vel_mean = final_net_vel_final,
         amp = final_amp_final)

RYI_vel_daily <- subset(RYI_vel_daily, select = c(1, 21:28))

#write out to .csv

write_csv(RYI_vel_daily, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYI_vel_daily.csv")

#import .csv

RYI_vel_daily <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYI_vel_daily.csv")

#repeat for RYF----------------------

siteNumbers <- c("11455385")  # 2018-07-03 	2022-07-17

startDate <- ""

endDate <- ""

parameterCd <- "72255"

uvRYF <- readNWISuv(siteNumbers,parameterCd,startDate,endDate)

#convert to PST

hrs <- 8 * 60 * 60

uvRYF$dateTime_PST <- as.POSIXct(uvRYF$dateTime, format = "%Y-%m-%d %H:%M:%S") - hrs

uvRYF <- subset(uvRYF, select = c(4, 7))

uvRYF <- uvRYF %>%
  rename(velocity_ft_s = "X_72255_00000")

#write out to .csv

write_csv(uvRYF, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvRYF_vel.csv")

#import .csv

uvRYF <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvRYF_vel.csv")

#apply godin filter

Tgodinfn<-function(tint,xint) {

  # apply Godin filter
  # it includes three passes in time domain, with windows 24, 24, and 25 hrs long
  # will assume 15-min data here, so
  # for first pass, 48 points have to be skipped and last 44
  # window is 12-1-11 hrs in length, i.e. use previous 12 hours and
  # subsequent 11 hrs to define result at 13th hour.

  xfilt<-rep_len(NA,length.out=length(xint))

  # do first pass, with window of 12,1,11 hrs
  for (i in seq(from=49,to=(length(xint)-44))) {
    xfilt[i]=sum(xint[(i-48):(i+44)])/93
  }

  xfilt1<-xfilt
  # now 2nd pass, same approach but window is 11,1,12 hrs
  for (i in seq(from=45, to=(length(xint)-48))) {
    xfilt1[i]=sum(xfilt[(i-44):(i+48)])/93;
  }
  xfilt=xfilt1;
  # now 3rd pass with 12-1-12 window
  for (i in seq(from=49,to=(length(xint)-48))) {
    xfilt[i]=sum(xfilt1[(i-48):(i+48)])/97;
  }
  # Now throw away 36 hrs at each end
  xfilt[1:(36*4)]<-NA
  isize=length(xfilt)
  xfilt[(isize-36*4):isize]<-NA

  return(xfilt)

}

#apply godin filter------------------------
uvRYF$vel_tf <- Tgodinfn(uvRYF$dateTime_PST, uvRYF$velocity_ft_s)

summary(uvRYF$vel_tf)

#tidal flow = velocity - net velocity(tidally-filtered)

uvRYF$tidal_vel <- uvRYF$velocity_ft_s - uvRYF$vel_tf

#plots - take a minute to load
#plot(uvRYF$dateTime_PST, uvRYF$velocity_ft_s)
#lines(uvRYF$dateTime_PST, uvRYF$tidal_vel, col = 3)
#lines(uvRYF$dateTime_PST, uvRYF$vel_tf, col = 2)
#legend("topright", legend = c("velocity", "tidal_vel", "net_vel"), col = c(1:3), lwd = 2)

#calculate amplitude- first break df in two
uvRYF_a <- subset(uvRYF, select = c(1, 2, 3))
uvRYF_b <- subset(uvRYF, select = c(2, 4))

#calcualte 30 hr moving window
uvRYF_b$roll_max <- rollapply(uvRYF_b$tidal_vel, 120, max, by = 1, partial = TRUE)
uvRYF_b$roll_min <- rollapply(uvRYF_b$tidal_vel, 120, min, by = 1, partial = TRUE)

#add high and low amplitude and divide by two
uvRYF_b$amp <- (uvRYF_b$roll_max + abs(uvRYF_b$roll_min))/2

#plots - take a minute to load
#plot(uvRYF_b$dateTime_PST, uvRYF_b$tidal_vel, col = 1, ylim = c(-5, 5), xlim = c(as.POSIXct('2020-10-01 17:30:00', format="%Y-%m-%d %H:%M:%S"), as.POSIXct('2020-10-15 17:30:00', format="%Y-%m-%d %H:%M:%S")))
#lines(uvRYF_b$dateTime_PST, uvRYF_b$roll_max, col = 2)
#lines(uvRYF_b$dateTime_PST, uvRYF_b$roll_min, col = 3)
#lines(uvRYF_b$dateTime_PST, uvRYF_b$amp, col = 5)
#legend("top", inset = c(-0.45, 0), legend = c("tidal_vel", "env_max", "env_min", "amp"), col = c(1:3, 5), lwd = 2)

#merge two dfs back together

uvRYF_c <- merge(uvRYF_a, uvRYF_b, by = 'dateTime_PST')

#subset

uvRYF_c <- subset(uvRYF_c, select = c(1:4, 7))

#calculate ratio of mean amplitude: mean tidal velocity

uvRYF_c$ratio <- uvRYF_c$amp/uvRYF_c$vel_tf

#>1 tidal influence, <1 river influence

#plot(uvRYF_c$dateTime_PST, uvRYF_c$ratio, ylim = c(-1000, 1000))

write_csv(uvRYF_c, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYF_vel.csv")
#downstep to daily - add n column

dvRYF <- uvRYF_c %>%
  mutate(date = as.Date(dateTime_PST)) %>%
  group_by(date) %>%
  summarise(mean_ratio = mean(ratio), mean_tidal_vel = mean(tidal_vel), mean_net_vel = mean(vel_tf), max_abs_velocity = max(abs(velocity_ft_s), na.rm = TRUE), min_vel_ft_s = min(velocity_ft_s, na.rm = TRUE), max_vel_ft_s = max(velocity_ft_s, na.rm = TRUE), mean_vel_ft_s = mean(velocity_ft_s, na.rm = TRUE), mean_amp = mean(amp), n_vel = n())

continous.dates <- data.frame (x = 1:1476, date = seq(as.Date('2018-07-03'),as.Date('2022-07-17'), by='day'))

RYF_vel_daily <- merge(continous.dates, dvRYF, by = "date", all = TRUE)

plot(RYF_vel_daily$date, RYF_vel_daily$n_vel)

plot(RYF_vel_daily$date, RYF_vel_daily$mean_vel_ft_s)

#if n_flow is NA replace with 0

RYF_vel_daily$n_vel[is.na(RYF_vel_daily$n_vel)] <- 0

sum(RYF_vel_daily$n_vel<=91)#24 days with <95% of flow measurements

#add column to identify if flow data is measured or will be imputed

RYF_vel_daily$group <- ifelse(RYF_vel_daily$n_vel>= 91, "measure", "impute")

#if n_value is <91, change mean, max, and min velocity to NA

RYF_vel_daily$mean_final_vel <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$mean_vel_ft_s, NA)

RYF_vel_daily$max_final_vel <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$max_vel_ft_s, NA)

RYF_vel_daily$min_final_vel <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$min_vel_ft_s, NA)

RYF_vel_daily$max_abs_final_vel <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$max_abs_velocity, NA)

RYF_vel_daily$final_tidal_vel <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$mean_tidal_vel, NA)

RYF_vel_daily$final_net_vel <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$mean_net_vel, NA)

RYF_vel_daily$final_ratio <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$mean_ratio, NA)

RYF_vel_daily$final_amp <- ifelse(RYF_vel_daily$n_vel>= 91, RYF_vel_daily$mean_amp, NA)

#impute missing values

RYF_vel_daily$final_mean_vel_final <- na_ma(RYF_vel_daily$mean_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYF_vel_daily$final_max_vel_final <- na_ma(RYF_vel_daily$max_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYF_vel_daily$final_min_vel_final <- na_ma(RYF_vel_daily$min_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYF_vel_daily$final_max_abs_final <- na_ma(RYF_vel_daily$max_abs_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYF_vel_daily$final_ratio_final <- na_ma(RYF_vel_daily$final_ratio, k = 7, weighting = "exponential", maxgap = Inf)

RYF_vel_daily$final_tidal_vel_final <- na_ma(RYF_vel_daily$final_tidal_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYF_vel_daily$final_net_vel_final <- na_ma(RYF_vel_daily$final_net_vel, k = 7, weighting = "exponential", maxgap = Inf)

RYF_vel_daily$final_amp_final <- na_ma(RYF_vel_daily$final_amp, k = 7, weighting = "exponential", maxgap = Inf)

summary(RYF_vel_daily)

plot(RYF_vel_daily$date, RYF_vel_daily$final_mean_vel_final, ylim = c(-4, 6))
lines(RYF_vel_daily$date, RYF_vel_daily$final_max_vel_final, col = 2)
lines(RYF_vel_daily$date, RYF_vel_daily$final_min_vel_final, col = 3)
lines(RYF_vel_daily$date, RYF_vel_daily$final_max_abs_final, col = 4)

RYF_vel_daily <- RYF_vel_daily %>%
  rename(mean_vel = final_mean_vel_final,
         max_vel = final_max_vel_final,
         min_vel = final_min_vel_final,
         max_abs_vel = final_max_abs_final,
         mean_tidal_vel = final_tidal_vel_final,
         ratio_mean = final_ratio_final,
         net_vel_mean = final_net_vel_final,
         amp = final_amp_final)

RYF_vel_daily <- subset(RYF_vel_daily, select = c(1, 21:28))

#write out to .csv

write_csv(RYF_vel_daily, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYF_vel_daily.csv")

#merge RYI and RYF datasets---------------------
#first rename columns

RYF_vel_daily <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYF_vel_daily.csv")

RYI_vel_daily <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYI_vel_daily.csv")

RYI_vel_daily <- RYI_vel_daily %>%
  mutate(RYI_mean_vel = mean_vel,
         RYI_max_vel = max_vel,
         RYI_min_vel = min_vel,
         RYI_max_abs_vel = max_abs_vel,
         RYI_mean_tide_vel = mean_tide_vel,
         RYI_ratio_mean = ratio_mean,
         RYI_net_vel_mean = net_vel_mean,
         RYI_amp = amp
         )

RYF_vel_daily <- RYF_vel_daily %>%
  mutate(RYF_mean_vel = mean_vel,
         RYF_max_vel = max_vel,
         RYF_min_vel = min_vel,
         RYF_max_abs_vel = max_abs_vel,
         RYF_mean_tide_vel = mean_tide_vel,
         RYF_ratio_mean = ratio_mean,
         RYF_net_vel_mean = net_vel_mean,
         RYF_amp = amp
  )

RYI_RYF <- full_join(RYI_vel_daily, RYF_vel_daily, na.rm = TRUE, by = "date")
summary(RYI_RYF)

#run when using read.csv function
RYI_RYF <- RYI_RYF[-c(2:9, 18:25)]

#create final vel vector with if else statement: if RYF is na then fill in with RYI

RYI_RYF$mean_vel <- with(RYI_RYF, ifelse(is.na(RYF_mean_vel), RYI_mean_vel, RYF_mean_vel))

RYI_RYF$min_vel <- with(RYI_RYF, ifelse(is.na(RYF_min_vel), RYI_min_vel, RYF_min_vel))

RYI_RYF$max_vel <- with(RYI_RYF, ifelse(is.na(RYF_max_vel), RYI_max_vel, RYF_max_vel))

RYI_RYF$max_abs_vel <- with(RYI_RYF, ifelse(is.na(RYF_max_abs_vel), RYI_max_abs_vel, RYF_max_abs_vel))

RYI_RYF$mean_tide_vel <- with(RYI_RYF, ifelse(is.na(RYF_mean_tide_vel), RYI_mean_tide_vel, RYF_mean_tide_vel))

RYI_RYF$ratio_mean <- with(RYI_RYF, ifelse(is.na(RYF_ratio_mean), RYI_ratio_mean, RYF_ratio_mean))

RYI_RYF$net_vel_mean <- with(RYI_RYF, ifelse(is.na(RYF_net_vel_mean), RYI_net_vel_mean, RYF_net_vel_mean))

RYI_RYF$amp <- with(RYI_RYF, ifelse(is.na(RYF_amp), RYI_amp, RYF_amp))

RYI_RYF <- rename(RYI_RYF, Date = date)

RYI_RYF <- subset(RYI_RYF, select = c(1, 18:25))

#create sign column

RYI_RYF <- RYI_RYF %>%
  mutate(sign=case_when(abs(min_vel) > max_vel ~ "-", abs(min_vel) < max_vel ~ "+"))

RYI_RYF$station <- "Cache"
#write to .csv

write_csv(RYI_RYF, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYI_RYF.csv")

#extra plots and ANOVA below

#import .csv

RYI_RYF <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYI_RYF.csv")

#below used for ANOVA
#merge YearAdj and season

Cache <- merge(RYI_RYF, raw_hydro_1975_2021, by = "Date")

Cache$log_out <- log(Cache$Outflow)

#plot velocity across years

boxplot(Cache$mean_vel ~ Cache$YearAdj)

boxplot(Cache$abs_min_vel ~ Cache$YearAdj)
boxplot(Cache$max_vel ~ Cache$YearAdj)
boxplot(Cache$min_vel ~ Cache$YearAdj)

#merge WY type and drought
lt_seasonal_2 <- subset(lt_seasonal, select = -c(3, 6:11))

Cache_2 <- merge(Cache, lt_seasonal_2, by = c("YearAdj", "Season"))



#ready for ANOVA

save_dir<-file.path("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/output")

#make Season, YearAdj, and Drought factors

Cache_3 <- Cache_2 %>% mutate(Season=factor(Season, levels=c("Winter", "Spring", "Summer", "Fall")),YearAdj=factor(YearAdj),Drought=factor(Drought, levels=c("W", "N", "D")))

#plot explore

plot(Cache_3$Outflow, Cache_3$abs_vel)
plot(Cache_3$log_out, Cache_3$final_vel, col = ifelse(Cache_3$Drought == "D", "gold2", "green3"))
legend("topleft", lty = c(1,2), legend = c("D", "N"), col = c("gold2", "green3"))

plot(Cache_3$log_out, Cache_3$final_vel, col = Cache_3$Season)
legend("top", cex = 0.5, lty = 1, legend = c("Fall", "Winter", "Spring", "Summer"), col = c(1:4))

#ggplot explore

# this works
ggplot(Cache_3, aes(x=log_out, y=final_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  ylab("vel (ft/s)")+
  theme_bw()

#color palette works in legend but not in figs
ggplot(Cache_3, aes(x=log_out, y=final_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  wyt_pal()+
  ylab("vel (ft/s)")+
  theme_bw()

ggplot(Cache_3, aes(x=log_out, y=final_max_abs_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="fixed")+
  ylab("max abs(vel (ft/s))")+
  theme_bw()

ggplot(Cache_3, aes(x=log_out, y=final_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  ylab("vel (ft/s)")+
  theme_bw()

ggplot(Cache_3, aes(x=log_out, y=final_max_abs_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~YearAdj, scales="fixed")+
  ylab("max abs vel (ft/s)")+
  theme_bw()


ggplot(Cache_3, aes(x=log_out, y=final_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="fixed")+
  ylab("abs vel (ft/s)")+
  theme_bw()

ggplot(Cache_3, aes(x=log_out, y=final_max_abs_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="fixed")+
  wyt_pal()+
  ylab("abs vel (ft/s)")+
  theme_bw()

#beeswarm

install.packages('ggbeeswarm')
library(ggbeeswarm)

ggplot(Cache_3, aes(x=log_out, y=final_vel, fill=YearType))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  ylab("vel (ft/s)")+
  theme_bw() +
  geom_quasirandom(              # jitters the points to show the distribution
    shape=21,                    # sets the shape of the points
    size=2,                      # changes the point size
    alpha=0.75,                  # makes the points 75% transparent
    color = "slategrey",         # changes the color (outline)
    fill = "steelblue") +        # changes the point fill
  stat_summary(                  # adds a statistical summary of the data
    fun = median,                # fun, fun.min, fun.max used to add the median
    fun.min = median,
    fun.max = median,
    geom = "crossbar",           # sets the shape of the median to a crossbar
    width = 0.8,                 # width of the median
    size = 0.8,                  # size (thickness) of the median
    color = "black") +           # color of the median
  xlab("log_outflow") +             # adds a custom label for the x-axis
  theme_classic(base_size=20)



#load Sam's function for year + season analysis

model_plotter_cache<-function(model, data){
  data<-data%>%
    mutate(Residuals=resid(model),
           Fitted=predict(model))

  p_hist<-ggplot(data, aes(x=Residuals))+
    geom_histogram()+
    xlab("cache (ft/s)")+
    theme_bw()

  p_res_fit<-ggplot(data, aes(x=Residuals, y=Fitted))+
    geom_point()+
    ylab("Predicted cache (ft/s)")+
    xlab("Residuals cache (ft/s)")+
    theme_bw()

  p_obs_fit<-ggplot(data, aes(x=final_max_abs_vel, y=Fitted))+
    geom_point()+
    geom_abline(slope=1, intercept=0, color="red")+
    ylab("Predicted cache (ft/s)")+
    xlab("Residuals cache (ft/s)")+
    theme_bw()

  vel<-(p_hist+plot_layout(ncol=1))+(p_res_fit+p_obs_fit+plot_layout(ncol=2))+plot_layout(nrow=2, widths=c(1, 0.5, 0.5))

  return(vel)
}

#load Sam's tukey function for year

tukey_plotter_cache_yr<-function(model, data, data_type, model_type){

  tuk<-emmeans(model, list(data=data_type, model=model_type))

  tuk_data<-as_tibble(cld(tuk$data, sort=FALSE, Letters = letters))%>%
    mutate(.group=str_remove_all(.group, fixed(" ")))%>%
    left_join(data%>%
                group_by(across(all_of(data_type)))%>%
                summarise(max_abs_vel=max(final_max_abs_vel), .groups="drop"),
              by=data_type)

  tuk_model<-as_tibble(cld(tuk$model, sort=FALSE, Letters = letters))%>%
    mutate(.group=str_remove_all(.group, fixed(" ")))%>%
    left_join(data%>%
                group_by(across(all_of(model_type)))%>%
                summarise(max_abs_vel=max(final_max_abs_vel), .groups="drop"),
              by=model_type)

  p_data<-ggplot(tuk_data, aes(x=.data[[data_type]], y=emmean, ymin=lower.CL, ymax=upper.CL, label=.group))+
    geom_boxplot(data=data, aes(x=.data[[data_type]], y=final_vel), inherit.aes = FALSE)+
    geom_pointrange(color="red", position=position_nudge(x=0.1))+
    geom_text(aes(y=final_max_abs_vel+(max(data$final_max_abs_vel)-min(data$final_max_abs_vel))/20), size=6)+
    ylab("cache (ft/s)")+
    theme_bw(base_size=16)

  p_model<-ggplot(tuk_model, aes(x=.data[[model_type]], y=emmean, ymin=lower.CL, ymax=upper.CL, label=.group))+
    geom_boxplot(data=data, aes(x=.data[[model_type]], y=final_vel), inherit.aes = FALSE)+
    geom_pointrange(color="red", position=position_nudge(x=0.1))+
    geom_text(aes(y=final_max_abs_vel+(max(data$final_max_abs_vel)-min(data$final_max_abs_vel))/20), angle=if_else(model_type=="YearAdj", 90, 0), hjust=if_else(model_type=="YearAdj", "left", NA_character_), vjust=0.25, size=6)+
    ylab("cache (ft/s)")+
    theme_bw(base_size=16)+
    {if(model_type=="YearAdj"){
      list(geom_tile(data=data,
                     aes(x=YearAdj, y=min(final_max_abs_vel)-(max(final_max_abs_vel)-min(final_max_abs_vel))/20,
                         fill=Drought, height=(max(final_max_abs_vel)-min(final_max_abs_vel))/20),
                     inherit.aes = FALSE),
           xlab("Year"),
           theme(axis.text.x=element_text(angle=90, vjust=0.5)),
           scale_y_continuous(expand = expansion(mult=c(0,0.1))),
           drt_color_pal_drought())
    }}

  vel<-p_data/p_model+plot_annotation(tag_levels="A")

  if(model_type=="YearAdj"){
    vel<-vel+plot_layout(heights = c(0.8, 1))
  }

  return(vel)
}

#Parameter: USGS Cache velocity
#ANOVA:  Metric ~ Drought/Wet + Season

#Data ??? seasonal averages from 2011-2021

#ANOVA:  Metric ~ factor(Year) + Season

#Data ??? seasonal averages from 2011-2021

#Plot vel by year

ggplot(Cache_3, aes(x=YearAdj, y=final_vel, color=YearType))+
  #geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  wyt_pal()+
  ylab("vel (ft/s)")+
  theme_bw()+
geom_quasirandom(              # jitters the points to show the distribution
  #shape=21,                    # sets the shape of the points
  size=2,                      # changes the point size
  alpha=0.75,                  # makes the points 75% transparent
  #color = "slategrey",         # changes the color (outline)
  #fill = "steelblue"
  ) +        # changes the point fill
  stat_summary(                  # adds a statistical summary of the data
    fun = median,                # fun, fun.min, fun.max used to add the median
    fun.min = median,
    fun.max = median,
    geom = "crossbar",           # sets the shape of the median to a crossbar
    width = 0.8,                 # width of the median
    size = 0.8,                  # size (thickness) of the median
    color = "black") +           # color of the median
  xlab("log_outflow") +             # adds a custom label for the x-axis
  theme_classic(base_size=20)

#Plot Cache velocity by Drought Index

ggplot(Cache_3, aes(x=Season, y=final_vel, fill=Drought))+
  geom_boxplot()+
  facet_wrap(~Season, scales="fixed")+
  drt_color_pal_drought()+
  ylab("vel (ft/s)")+
  theme_bw()

ggplot(Cache_3, aes(x=Season, y=final_max_abs_vel, fill=Drought))+
  geom_boxplot()+
  facet_wrap(~Season, scales="fixed")+
  drt_color_pal_drought()+
  ylab("vel (ft/s)")+
  theme_bw()

# ANOVA for velocity by YearAdj + Season

m_cache_year <- aov(final_max_abs_vel ~ factor(YearAdj) + Season, data = Cache_3)

#Estimated effects may be unbalanced with factor(YearAdj)*Season
m_cache_year

m_cache_year_Anova <- Anova(m_cache_year, type=2)

m_cache_year_Anova

# check residuals

#hist(resid(m_out_year))

#fitted values of the model

data.fit = fitted(m_cache_year)

#residuals of the model

data.res = resid(m_cache_year)

#look at residuals

hist(data.res)
hist(data.fit)

data.stdres = rstandard(m_cache_year)

#make qq plot and add line. first obtain normal prob plot for standardized resids

data.stdres = rstandard(m_cache_year)

qqnorm(data.stdres)

qqline(data.stdres)

#check assumptions using Sam's function

model_plotter_cache(m_cache_year, Cache_3)

#anova table

summary(m_cache_year)

#post - hoc using Sam's tukey function

cache_year_tukey <-tukey_plotter_cache_yr(m_cache_year, Cache_3, "Season", "YearAdj")

cache_year_tukey

ggsave(plot=cache_year_tukey, filename=file.path(save_dir, "cache_year_tukey.png"), device="png", height=12, width=15, units="in")

# ANOVA for Cache velocity by Drought + Season

m_cache_dr <- aov(final_max_abs_vel ~ Drought + Season, data = Cache_3)

summary(m_cache_dr)

m_cache_dr_Anova <- Anova(m_cache_dr, type = 2)

# check residuals

#hist(resid(ex_year))

#fitted values of the model

data.fit = fitted(m_cache_dr)

#residuals of the model

data.res = resid(m_cache_dr)

#look at residuals

hist(data.res)
hist(data.fit)

data.stdres = rstandard(m_cache_dr)

#make qq plot and add line. first obtain normal prob plot for standardized resids

data.stdres = rstandard(m_cache_dr)

qqnorm(data.stdres)

qqline(data.stdres)

#check assumptions using Sam's function

model_plotter_cache(m_cache_dr, Cache_3)

#anova table

summary(m_cache_dr)

#post - hoc using Sam's tukey function

cache_drought_tukey <-tukey_plotter_cache_yr(m_cache_dr, Cache_3, "Season", "Drought")

cache_drought_tukey

ggsave(plot=cache_drought_tukey, filename=file.path(save_dir, "cache_season_drought_model.png"), device="png", height=12, width=15, units="in")

#save all Anova output
anovas<-bind_rows(
  mutate(as_tibble(m_cache_year_Anova, rownames = "Parameter"), model="Year_Season"),
  mutate(as_tibble(m_cache_dr_Anova, rownames = "Parameter"), model="Season_Drought")
)%>%
  mutate(`Pr(>F)`=if_else(`Pr(>F)`<0.001, "< 0.001", as.character(round(`Pr(>F)`, 4))))%>%
  write_csv(file.path(save_dir, "cache_vel_anovas.csv"))

#save all Anova output
anovas<-bind_rows(
  mutate(as_tibble(m_CCH_year_Anova, rownames = "Parameter"), model="Year_Season"),
  mutate(as_tibble(m_CCH_dr_Anova, rownames = "Parameter"), model="Season_Drought")
)

anovas$Metric <- c('cache_vel')

anovas <- anovas %>% relocate(Metric, .before = Parameter)%>%
  mutate(`Pr(>F)`=if_else(`Pr(>F)`<0.001, "< 0.001", as.character(round(`Pr(>F)`, 4))))%>%
  write_csv(file.path(save_dir, "cache_vel_anovas.csv"))

#compare 2021 to prior years
#load data
#raw_out <- raw_hydro_1975_2021

#raw_out$log_Outflow = log(raw_out$Outflow)

#adding drought_20_21 and yeartype_20_21 columns

cache<-Cache_3%>%
  filter(!is.na(final_vel))%>%
  #left_join(lt_regional%>%distinct(YearAdj, SVIndex, YearType, Drought),by="YearAdj")%>%
  mutate(across(c(Drought, YearType), list(`20_21`=~case_when(YearAdj==2021 ~ "2021",YearAdj==2020 ~ "2020",TRUE ~ as.character(.x)))),across(c(YearType, YearType_20_21), ~factor(.x, levels=c("2020", "2021", "Critical", "Dry", "Below Normal", "Above Normal", "Wet"))),Season=factor(Season, levels=c("Winter", "Spring", "Summer", "Fall")))

# graph how 2021 compares to Drought, Normal, and Wet periods?

cache_vel<-ggplot(cache, aes(x=Drought_20_21, y=final_vel, fill=Drought))+geom_boxplot()+drt_color_pal_drought()+xlab("Drought")+ylab("final_vel")+theme_bw()

cache_vel

ggsave(plot=cache_vel, filename=file.path(save_dir, "cache_vel_drought_20_21.png"), device="png", height=4, width=5, units="in")


##extra stuff---------add a few columns------------------

RYI_RYF$pos_vel <- ifelse(RYI_RYF$final_vel > "0", "1", "0")
RYI_RYF$pos_vel <- as.numeric(RYI_RYF$pos_vel)
RYI_RYF$neg_vel <- ifelse(RYI_RYF$final_vel < "0", "1", "0")
RYI_RYF$neg_vel <- as.numeric(RYI_RYF$neg_vel)
RYI_RYF$Month <- format(as.Date(RYI_RYF$date), "%m")

RYI_RYF$Year <- format(as.Date(RYI_RYF$date), "%Y")

#quick view of time series---------------

plot(RYI_RYF$date, RYI_RYF$final_vel)

hist(RYI_RYF$final_vel)

write.csv(RYI_RYF, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/RYI_RYF.csv")

#assign adjusted water year and season columns
#run ANOVA as completed with dayflow stuff

install.packages("zoo")
library(zoo)

RYI_RYF$Month_Yr <- format(as.Date(RYI_RYF$date), "%Y-%m")

RYI_RYF_month <- RYI_RYF %>%
  group_by(Month_Yr) %>%
  summarize(mean_vel = mean(final_vel, na.rm=TRUE), sd_vel = sd(final_vel,na.rm=TRUE), n_flow = n(), pos_vel = sum(pos_vel), neg_vel = sum(neg_vel))

RYI_RYF_month$Month_Yr <- as.yearmon(RYI_RYF_month$Month_Yr)

#RYI_RYF_month$Month_Yr <- as.POSIXct(RYI_RYF_month$Month_Yr, format = "%Y-%m")

plot(RYI_RYF_month$Month_Yr, RYI_RYF_month$pos_vel)
lines(RYI_RYF_month$Month_Yr, RYI_RYF_month$neg_vel, col = "3")
#YearAdj is the adjusted calendar year (December-November, with December of the
#previous calendar year included with the following year

