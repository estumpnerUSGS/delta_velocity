#retrieval, imputation, and summary of instantaneous USGS velocity data
#USGS 11455420 SACRAMENTO R A RIO VISTA CA

library(dataRetrieval)
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(imputeTS)
library(zoo)

#data retrieval pull through NWIS web services

siteNumbers <- c("11455420")

startDate <- ""

endDate <- ""

parameterCd <- "72255"

uvSRV <- readNWISuv(siteNumbers,parameterCd,startDate,endDate)

#convert to PST

hrs <- 8 * 60 * 60

uvSRV$dateTime_PST <- as.POSIXct(uvSRV$dateTime, format = "%d/%m/%Y %H:%M:%S") - hrs

uvSRV_new <- subset(uvSRV, select = c(4, 7))

uvSRV_new <- uvSRV_new %>%
  rename(velocity_ft_s = "X_72255_00000")

write_csv(uvSRV_new, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvSRV_vel.csv")

SRV_vel <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvSRV_vel.csv")

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
SRV_vel$vel_tf <- Tgodinfn(uvSJJ$dateTime_PST, SRV_vel$velocity_ft_s)

#tidal flow = velocity - net velocity(tidally-filtered)

SRV_vel$tidal_vel <- SRV_vel$velocity_ft_s - SRV_vel$vel_tf

#plots - take a minute to load
#plot(SRV_vel$dateTime_PST, SRV_vel$velocity_ft_s)
#lines(SRV_vel$dateTime_PST, SRV_vel$tidal_vel, col = 3)
#lines(SRV_vel$dateTime_PST, SRV_vel$vel_tf, col = 2)

#calculate amplitude - first break df in two
SRV_vel_a <- subset(SRV_vel, select = c(1, 2, 3))
SRV_vel_b <- subset(SRV_vel, select = c(2, 4))

#calcualte 30 hr moving window
SRV_vel_b$roll_max <- rollapply(SRV_vel_b$tidal_vel, 120, max, by = 1, partial = TRUE)
SRV_vel_b$roll_min <- rollapply(SRV_vel_b$tidal_vel, 120, min, by = 1, partial = TRUE)

#add high and low amplitude and divide by two
SRV_vel_b$amp <- (SRV_vel_b$roll_max + abs(SRV_vel_b$roll_min))/2

#plots - take a minute to load
#plot(SRV_vel_b$dateTime_PST, SRV_vel_b$tidal_vel, col = 1, ylim = c(-5, 5)) #xlim = c(as.POSIXct('2020-10-01 17:30:00', format="%Y-%m-%d %H:%M:%S"), as.POSIXct('2020-10-15 17:30:00', format="%Y-%m-%d %H:%M:%S")))
#lines(SRV_vel_b$dateTime_PST, SRV_vel_b$roll_max, col = 2)
#lines(SRV_vel_b$dateTime_PST, SRV_vel_b$roll_min, col = 3)
#lines(SRV_vel_b$dateTime_PST, SRV_vel_b$amp, col = 5)

#merge two dfs back together

SRV_vel_c <- merge(SRV_vel_a, SRV_vel_b, by = 'dateTime_PST')

#subset

SRV_vel_c <- subset(SRV_vel_c, select = c(1:4, 7))

#calculate ratio of mean amplitude: mean tidal velocity

SRV_vel_c$ratio <- SRV_vel_c$amp/SRV_vel_c$vel_tf

#>1 tidal influence, <1 river influence

#plot(SRV_vel_c$dateTime_PST, SRV_vel_c$ratio, ylim = c(-1000, 1000))

write_csv(SRV_vel_c, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/SRV_vel.csv")

SRV_vel_c <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/SRV_vel.csv")

#downstep to daily - add n column
dvSRV <- SRV_vel_c %>%
  mutate(date = as.Date(dateTime_PST)) %>%
  group_by(date) %>%
  summarise(mean_ratio = mean(ratio), mean_tidal_vel = mean(tidal_vel), mean_net_vel = mean(vel_tf), max_abs_vel_ft_s = max(abs(velocity_ft_s), na.rm = TRUE), min_vel_ft_s = min(velocity_ft_s, na.rm = TRUE), max_vel_ft_s = max(velocity_ft_s, na.rm = TRUE), mean_vel_ft_s = mean(velocity_ft_s, na.rm = TRUE), mean_amp = mean(amp), n_vel = n())

continous.dates <- data.frame (x = 1:5418, date = seq(as.Date('2007-10-01'),as.Date('2022-07-31'), by='day'))

SRV_vel_daily <- merge(continous.dates, dvSRV, by = "date", all = TRUE)

plot(SRV_vel_daily$date, SRV_vel_daily$n_vel)

plot(SRV_vel_daily$date, SRV_vel_daily$mean_vel_ft_s)

#if n_flow is NA replace with 0

SRV_vel_daily$n_vel[is.na(SRV_vel_daily$n_vel)] <- 0

sum(SRV_vel_daily$n_vel<=91)#201 days with <95% of flow measurements

#add column to identify if flow data is measured or will be imputed

SRV_vel_daily$group <- ifelse(SRV_vel_daily$n_vel>= 91, "measure", "impute")

#if n_value is <91, change mean, max, and min velocity to NA

SRV_vel_daily$mean_final_vel <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$mean_vel_ft_s, NA)

SRV_vel_daily$max_final_vel <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$max_vel_ft_s, NA)

SRV_vel_daily$min_final_vel <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$min_vel_ft_s, NA)

SRV_vel_daily$max_abs_final_vel <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$max_abs_vel_ft_s, NA)

SRV_vel_daily$final_tidal_vel <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$mean_tidal_vel, NA)

SRV_vel_daily$final_net_vel <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$mean_net_vel, NA)

SRV_vel_daily$final_ratio <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$mean_ratio, NA)

SRV_vel_daily$final_amp <- ifelse(SRV_vel_daily$n_vel>= 91, SRV_vel_daily$mean_amp, NA)

#impute missing values

SRV_vel_daily$final_mean_vel_final <- na_ma(SRV_vel_daily$mean_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SRV_vel_daily$final_max_vel_final <- na_ma(SRV_vel_daily$max_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SRV_vel_daily$final_min_vel_final <- na_ma(SRV_vel_daily$min_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SRV_vel_daily$final_max_abs_final <- na_ma(SRV_vel_daily$max_abs_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SRV_vel_daily$final_ratio_final <- na_ma(SRV_vel_daily$final_ratio, k = 7, weighting = "exponential", maxgap = Inf)

SRV_vel_daily$final_tidal_vel_final <- na_ma(SRV_vel_daily$final_tidal_vel, k = 7, weighting = "exponential", maxgap = Inf)

SRV_vel_daily$final_net_vel_final <- na_ma(SRV_vel_daily$final_net_vel, k = 7, weighting = "exponential", maxgap = Inf)

SRV_vel_daily$final_amp_final <- na_ma(SRV_vel_daily$final_amp, k = 7, weighting = "exponential", maxgap = Inf)

#summary(SRV_vel_daily)

#plot(SRV_vel_daily$date, SRV_vel_daily$final_max_abs_final, ylim = c(-0.6, 2))
#lines(SRV_vel_daily$date, SRV_vel_daily$final_tidal_vel_final, col = 2)
#lines(SRV_vel_daily$date, SRV_vel_daily$final_net_vel_final, col = 3)

SRV_vel_daily <- SRV_vel_daily %>%
  rename(mean_vel = final_mean_vel_final,
         max_vel = final_max_vel_final,
         min_vel = final_min_vel_final,
         max_abs_vel = final_max_abs_final,
         mean_tide_vel = final_tidal_vel_final,
         ratio_mean = final_ratio_final,
         net_vel_mean = final_net_vel_final,
         amp = final_amp_final)

SRV_vel_daily <- subset(SRV_vel_daily, select = c(1, 21:28))

SRV_vel_daily <- rename(SRV_vel_daily, Date = date)

SRV_vel_daily$station <- "Rio Vista"

#create column to assign sign for max abs velocity column

SRV_vel_daily <- SRV_vel_daily %>%
  mutate(sign=case_when(abs(min_vel) > max_vel ~ "-", abs(min_vel) < max_vel ~ "+"))

write_csv(SRV_vel_daily, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/SRV_vel_daily.csv")
