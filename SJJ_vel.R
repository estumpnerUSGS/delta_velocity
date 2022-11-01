#retrieval, imputation, and summary of instantaneous USGS velocity data
#USGS 11337190 SAN JOAQUIN R A JERSEY POINT CA

library(dataRetrieval)
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(imputeTS)
library(zoo)

siteNumbers <- c("11337190")

startDate <- ""

endDate <- ""

parameterCd <- "72255"

uvSJJ <- readNWISuv(siteNumbers,parameterCd,startDate,endDate)

#convert to PST

hrs <- 8 * 60 * 60

uvSJJ$dateTime_PST <- as.POSIXct(uvSJJ$dateTime, format = "%Y-%m-%d %H:%M:%S") - hrs

uvSJJ <- subset(uvSJJ, select = c(4, 7))

uvSJJ <- uvSJJ %>%
  rename(velocity_ft_s = "X_72255_00000")

write.csv(uvSJJ, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvSJJ_vel.csv")

#import .csv

uvSJJ <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/uvSJJ_vel.csv")

#reformat df

uvSJJ <- subset(uvSJJ, select = c(2:3))

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
uvSJJ$vel_tf <- Tgodinfn(uvSJJ$dateTime_PST, uvSJJ$velocity_ft_s)

#tidal flow = velocity - net velocity(tidally-filtered)

uvSJJ$tidal_vel <- uvSJJ$velocity_ft_s - uvSJJ$vel_tf

#plots - take a minute to load
#plot(uvSJJ$dateTime_PST, uvSJJ$velocity_ft_s)
#lines(uvSJJ$dateTime_PST, uvSJJ$tidal_vel, col = 3)
#lines(uvSJJ$dateTime_PST, uvSJJ$vel_tf, col = 2)
#legend("topright", legend = c("velocity", "tidal_vel", "net_vel"), col = c(1:3), lwd = 2)

#calculate amplitude- first break df in two
uvSJJ_a <- subset(uvSJJ, select = c(1, 2, 3))
uvSJJ_b <- subset(uvSJJ, select = c(2, 4))

#calcualte 30 hr moving window
uvSJJ_b$roll_max <- rollapply(uvSJJ_b$tidal_vel, 120, max, by = 1, partial = TRUE)
uvSJJ_b$roll_min <- rollapply(uvSJJ_b$tidal_vel, 120, min, by = 1, partial = TRUE)

#add high and low amplitude and divide by two
uvSJJ_b$amp <- (uvSJJ_b$roll_max + abs(uvSJJ_b$roll_min))/2

#plots - take a minute to load
#plot(uvSJJ_b$dateTime_PST, uvSJJ_b$tidal_vel, col = 1, ylim = c(-5, 5), xlim = c(as.POSIXct('2020-10-01 17:30:00', format="%Y-%m-%d %H:%M:%S"), as.POSIXct('2020-10-15 17:30:00', format="%Y-%m-%d %H:%M:%S")))
#lines(uvSJJ_b$dateTime_PST, uvSJJ_b$roll_max, col = 2)
#lines(uvSJJ_b$dateTime_PST, uvSJJ_b$roll_min, col = 3)
#lines(uvSJJ_b$dateTime_PST, uvSJJ_b$amp, col = 5)
#legend("top", inset = c(-0.45, 0), legend = c("tidal_vel", "env_max", "env_min", "amp"), col = c(1:3, 5), lwd = 2)

#merge two dfs back together

uvSJJ_c <- merge(uvSJJ_a, uvSJJ_b, by = 'dateTime_PST')

#subset

uvSJJ_c <- subset(uvSJJ_c, select = c(1:4, 7))

#calculate ratio of mean amplitude: mean tidal velocity

uvSJJ_c$ratio <- uvSJJ_c$amp/uvSJJ_c$vel_tf

#>1 tidal influence, <1 river influence

#plot(uvSJJ_c$dateTime_PST, uvSJJ_c$ratio, ylim = c(-1000, 1000))

write_csv(uvSJJ_c, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/SJJ_vel.csv")

#downstep to daily - add n column

dvSJJ <- uvSJJ_c %>%
  mutate(date = as.Date(dateTime_PST)) %>%
  group_by(date) %>%
  summarise(mean_ratio = mean(ratio), mean_tidal_vel = mean(tidal_vel), mean_net_vel = mean(vel_tf), max_abs_velocity = max(abs(velocity_ft_s), na.rm = TRUE), min_vel_ft_s = min(velocity_ft_s, na.rm = TRUE), max_vel_ft_s = max(velocity_ft_s, na.rm = TRUE), mean_vel_ft_s = mean(velocity_ft_s, na.rm = TRUE), mean_amp = mean(amp), n_vel = n())

continous.dates <- data.frame (x = 1:4201, date = seq(as.Date('2007-10-01'),as.Date('2019-04-01'), by='day'))

SJJ_vel_daily <- merge(continous.dates, dvSJJ, by = "date", all = TRUE)

plot(SJJ_vel_daily$date, SJJ_vel_daily$n_vel)

plot(SJJ_vel_daily$date, SJJ_vel_daily$mean_vel_ft_s)

#if n_flow is NA replace with 0

SJJ_vel_daily$n_vel[is.na(SJJ_vel_daily$n_vel)] <- 0

sum(SJJ_vel_daily$n_vel<=91)#194 days with <95% of flow measurements

#add column to identify if flow data is measured or will be imputed

SJJ_vel_daily$group <- ifelse(SJJ_vel_daily$n_vel>= 91, "measure", "impute")

#if n_value is <91, change mean, max, and min velocity to NA

SJJ_vel_daily$mean_final_vel <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$mean_vel_ft_s, NA)

SJJ_vel_daily$max_final_vel <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$max_vel_ft_s, NA)

SJJ_vel_daily$min_final_vel <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$min_vel_ft_s, NA)

SJJ_vel_daily$max_abs_final_vel <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$max_abs_velocity, NA)

SJJ_vel_daily$final_tidal_vel <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$mean_tidal_vel, NA)

SJJ_vel_daily$final_net_vel <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$mean_net_vel, NA)

SJJ_vel_daily$final_ratio <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$mean_ratio, NA)

SJJ_vel_daily$final_amp <- ifelse(SJJ_vel_daily$n_vel>= 91, SJJ_vel_daily$mean_amp, NA)

#impute missing values

SJJ_vel_daily$final_mean_vel_final <- na_ma(SJJ_vel_daily$mean_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily$final_max_vel_final <- na_ma(SJJ_vel_daily$max_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily$final_min_vel_final <- na_ma(SJJ_vel_daily$min_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily$final_max_abs_final <- na_ma(SJJ_vel_daily$max_abs_final_vel, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily$final_ratio_final <- na_ma(SJJ_vel_daily$final_ratio, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily$final_tidal_vel_final <- na_ma(SJJ_vel_daily$final_tidal_vel, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily$final_net_vel_final <- na_ma(SJJ_vel_daily$final_net_vel, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily$final_amp_final <- na_ma(SJJ_vel_daily$final_amp, k = 7, weighting = "exponential", maxgap = Inf)

SJJ_vel_daily <- SJJ_vel_daily %>%
  rename(mean_vel = final_mean_vel_final,
         max_vel = final_max_vel_final,
         min_vel = final_min_vel_final,
         max_abs_vel = final_max_abs_final,
         mean_tide_vel = final_tidal_vel_final,
         ratio_mean = final_ratio_final,
         net_vel_mean = final_net_vel_final,
         amp = final_amp_final)


#subset

SJJ_vel_daily <- subset(SJJ_vel_daily, select = c(1, 21:28))

SJJ_vel_daily <- rename(SJJ_vel_daily, Date = date)

#create column to assign sign for max abs velocity column

SJJ_vel_daily <- SJJ_vel_daily %>%
  mutate(sign=case_when(abs(min_vel) > max_vel ~ "-", abs(min_vel) < max_vel ~ "+"))

SJJ_vel_daily$station <- "Jersey"

write_csv(SJJ_vel_daily, "C:/Users/estumpne/Documents/R/drought_synthesis_velocity/SJJ_vel_daily.csv")

#extra plots and ANOVA below

#import .csv

SJJ_vel_daily <- read_csv("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/SJJ_vel_daily.csv")

#merge YearAdj and season

Jersey <- merge(SJJ_vel_daily, raw_hydro_1975_2021, by = "Date")

lt_seas <- lt_seasonal[c(1, 4:5)]

lt_seas <- unique(lt_seas)

Jersey <- merge(Jersey, lt_seas, by = "YearAdj")

Jersey$log_out <- log(Jersey$Outflow)

#add wateryear to daily data frame

wtr_yr <- function(Date, start_month=10) {
  # Convert dates into POSIXlt
  dates.posix = as.POSIXlt(Date)
  # Year offset
  offset = ifelse(dates.posix$mon >= start_month - 1, 1, 0)
  # Water year
  adj.year = dates.posix$year + 1900 + offset
  # Return the water year
  adj.year
}

Jersey_WY <- Jersey %>%
  mutate(wtr_yr = wtr_yr(Date))

Jersey_WY <- Jersey_WY %>%
  group_by(wtr_yr) %>%
  mutate(wtr_day = (as.integer(difftime(Date,ymd(paste0(wtr_yr - 1 ,'-09-30')), units = "days"))))

Jersey_WY$log_export <- log(Jersey_WY$Export)

Jersey_WY %>%
  ggplot(aes(x=Date, y=log_export)) +
  geom_point(aes(x=Date, y=log_out))

plot(Jersey_WY$Date, Jersey_WY$log_export, ylim = c(-1, 15), xlab = "date", ylab = "cfs", main = "Jersey 2007 - 2021")
axis.Date(1, at=seq(min(Jersey_WY$Date), max(Jersey_WY$Date), by="year"), format="%Y", las = 2)
lines(Jersey_WY$Date, Jersey_WY$log_out, col = 2)
lines(Jersey_WY$Date, Jersey_WY$amp, ylim = c(-2,3), col = 3)
lines(Jersey_WY$Date, Jersey_WY$net_vel_mean, col = 4) +
lines(Jersey_WY$Date, Jersey_WY$mean_tide_vel, col = 5)
legend("top", inset = c(-0.45, 0), legend = c("log_exp", "log_out", "amp", "net_vel", "tide_vel"), col = c(1:5), lwd = 2)

Jersey_WY %>%
  ggplot(aes(x=wtr_day, y=amp, color = YearType)) +
  geom_smooth() +
  theme(text = element_text (size = 30)) +
  labs(y = "amp (ft)", xlab = "WY day", title = "Jersey Pt. amplitude 2007 - 2021")

Jersey_WY %>%
  ggplot(aes(x=wtr_day, y=net_vel_mean, color = YearType)) +
  geom_smooth() +
  theme(text = element_text (size = 30)) +
  labs(y = "net vel (cfs)", xlab = "WY day", title = "Jersey Pt. net vel. 2007 - 2021")

Jersey_WY %>%
  ggplot(aes(x=wtr_day, y=mean_tide_vel, color = YearType)) +
  geom_smooth() +
  theme(text = element_text (size = 30)) +
  labs(y = "tide vel (cfs)", xlab = "WY day", title = "Jersey Pt. tide vel. 2007 - 2021")

Jersey_WY %>%
  ggplot(aes(x=wtr_day, y=ratio_mean, color = YearType)) +
  geom_smooth() +
  theme(text = element_text (size = 30)) +
  labs(y = "ratio", xlab = "WY day", title = "Jersey Pt. ratio amp:net vel. 2007 - 2021")
#beeswarm 1 - max(abs(velocity))

Jersey_WY %>%
  ggplot(aes(x=YearType, y=amp)) +
  geom_quasirandom() +
  facet_wrap(~Season, scales = "fixed")+
  labs(y = "max abs vel (cfs)", x = "Year type", title = "SJJ Max abs velocity") +
  stat_summary(                  # adds a statistical summary of the data
    fun = median,                # fun, fun.min, fun.max used to add the median
    #fun.min = median,
    #fun.max = median,
    geom = "crossbar",           # sets the shape of the median to a crossbar
    width = 0.8,                 # width of the median
    size = 0.8,                  # size (thickness) of the median
    color = "black") +
  theme (text = element_text (size = 30))# color of the median


#beeswarm 2
#mean of tidal velocity amplitude/mean of net velocity
#>1 tide <1 river

Jersey_WY %>%
  ggplot(aes(x=YearType, y=ratio_mean, color = sign)) +
  coord_cartesian(ylim = c(-250, 250)) +
  geom_quasirandom() +
  facet_wrap(~Season, scales = "fixed")+
  labs(y = "max abs vel (cfs)", x = "Year type", title = "SJJ ratio of amplitude/net velocity") +
  geom_hline(yintercept=1, linetype='solid', col = 'black')+
  stat_summary(                  # adds a statistical summary of the data
    fun = median,                # fun, fun.min, fun.max used to add the median
    fun.min = median,
    fun.max = median,
    geom = "crossbar",           # sets the shape of the median to a crossbar
    width = 0.8,                 # width of the median
    size = 0.8,                  # size (thickness) of the median
    color = "black") +           # color of the median
  theme (text = element_text (size = 30))

#plot velocity across years

boxplot(Jersey$final_vel ~ Jersey$YearAdj)

#merge WY type and drought
lt_seasonal_2 <- subset(lt_seasonal, select = -c(3, 6:11))

Jersey_2 <- merge(Jersey, lt_seasonal_2, by = c("YearAdj", "Season"))



#ready for ANOVA

save_dir<-file.path("C:/Users/estumpne/Documents/R/drought_synthesis_velocity/output")

#make Season, YearAdj, and Drought factors

Jersey_3 <- Jersey_2 %>% mutate(Season=factor(Season, levels=c("Winter", "Spring", "Summer", "Fall")),YearAdj=factor(YearAdj),Drought=factor(Drought.x, levels=c("W", "N", "D")))

#plot explore

plot(Jersey_3$Outflow, Jersey_3$abs_vel)

plot(Jersey_3$log_out, Jersey_3$final_vel, col = ifelse(Jersey_3$Drought == "D", "gold2", "green3"))
legend("topleft", lty = c(1,2), legend = c("D", "N"), col = c("gold2", "green3"))

plot(Jersey_3$log_out, Jersey_3$max_abs_vel, col = ifelse(Jersey_3$Drought == "D", "gold2", "green3"))
legend("topleft", lty = c(1,2), legend = c("D", "N"), col = c("gold2", "green3"))

plot(Jersey_3$log_out, Jersey_3$final_vel, col = Jersey_3$Season)
legend("top", cex = 0.5, lty = 1, legend = c("Fall", "Winter", "Spring", "Summer"), col = c(1:4))

#ggplot explore

ggplot(Jersey_3, aes(x=log_out, y=final_vel, fill=Drought))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  drt_color_pal_drought()+
  ylab("vel (ft/s)")+
  theme_bw()

ggplot(Jersey_3, aes(x=log_out, y=max_abs_vel, fill=Drought))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  drt_color_pal_drought()+
  ylab("max abs vel (ft/s)")+
  theme_bw()

ggplot(Jersey_3, aes(x=log_out, y=max_abs_vel, fill=Drought))+
  geom_point(shape=21, color="black")+
  facet_wrap(~YearAdj, scales="fixed")+
  drt_color_pal_drought()+
  ylab("vel (ft/s)")+
  theme_bw()

ggplot(Jersey_3, aes(x=log_out, y=final_vel, fill=YearAdj))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  ylab("vel (ft/s)")+
  theme_bw()

ggplot(Jersey_3, aes(x=log_out, y=abs_vel, fill=Drought))+
  geom_point(shape=21, color="black")+
  facet_wrap(~YearAdj, scales="fixed")+
  drt_color_pal_drought()+
  ylab("abs vel (ft/s)")+
  theme_bw()

ggplot(Jersey_3, aes(x=log_out, y=abs_vel, fill=Drought))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  drt_color_pal_drought()+
  ylab("abs vel (ft/s)")+
  theme_bw()

#load Sam's function for year + season analysis

model_plotter_jersey<-function(model, data){
  data<-data%>%
    mutate(Residuals=resid(model),
           Fitted=predict(model))

  p_hist<-ggplot(data, aes(x=Residuals))+
    geom_histogram()+
    xlab("jersey (ft/s)")+
    theme_bw()

  p_res_fit<-ggplot(data, aes(x=Residuals, y=Fitted))+
    geom_point()+
    ylab("Predicted jersey (ft/s)")+
    xlab("Residuals jersey (ft/s)")+
    theme_bw()

  p_obs_fit<-ggplot(data, aes(x=final_vel, y=Fitted))+
    geom_point()+
    geom_abline(slope=1, intercept=0, color="red")+
    ylab("Predicted jersey (ft/s)")+
    xlab("Residuals jersey (ft/s)")+
    theme_bw()

  vel<-(p_hist+plot_layout(ncol=1))+(p_res_fit+p_obs_fit+plot_layout(ncol=2))+plot_layout(nrow=2, widths=c(1, 0.5, 0.5))

  return(vel)
}

#load Sam's tukey function for year

tukey_plotter_jersey_yr<-function(model, data, data_type, model_type){

  tuk<-emmeans(model, list(data=data_type, model=model_type))

  tuk_data<-as_tibble(cld(tuk$data, sort=FALSE, Letters = letters))%>%
    mutate(.group=str_remove_all(.group, fixed(" ")))%>%
    left_join(data%>%
                group_by(across(all_of(data_type)))%>%
                summarise(max_vel=max(final_vel), .groups="drop"),
              by=data_type)

  tuk_model<-as_tibble(cld(tuk$model, sort=FALSE, Letters = letters))%>%
    mutate(.group=str_remove_all(.group, fixed(" ")))%>%
    left_join(data%>%
                group_by(across(all_of(model_type)))%>%
                summarise(max_vel=max(final_vel), .groups="drop"),
              by=model_type)

  p_data<-ggplot(tuk_data, aes(x=.data[[data_type]], y=emmean, ymin=lower.CL, ymax=upper.CL, label=.group))+
    geom_boxplot(data=data, aes(x=.data[[data_type]], y=final_vel), inherit.aes = FALSE)+
    geom_pointrange(color="red", position=position_nudge(x=0.1))+
    geom_text(aes(y=max_vel+(max(data$final_vel)-min(data$final_vel))/20), size=6)+
    ylab("jersey (ft/s)")+
    theme_bw(base_size=16)

  p_model<-ggplot(tuk_model, aes(x=.data[[model_type]], y=emmean, ymin=lower.CL, ymax=upper.CL, label=.group))+
    geom_boxplot(data=data, aes(x=.data[[model_type]], y=final_vel), inherit.aes = FALSE)+
    geom_pointrange(color="red", position=position_nudge(x=0.1))+
    geom_text(aes(y=max_vel+(max(data$final_vel)-min(data$final_vel))/20), angle=if_else(model_type=="YearAdj", 90, 0), hjust=if_else(model_type=="YearAdj", "left", NA_character_), vjust=0.25, size=6)+
    ylab("jersey (ft/s)")+
    theme_bw(base_size=16)+
    {if(model_type=="YearAdj"){
      list(geom_tile(data=data,
                     aes(x=YearAdj, y=min(final_vel)-(max(final_vel)-min(final_vel))/20,
                         fill=Drought, height=(max(final_vel)-min(final_vel))/20),
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

#Parameter: USGS Jersey Velocity
#ANOVA:  Metric ~ Drought/Wet + Season

#Data – seasonal averages from 2011-2021

#ANOVA:  Metric ~ factor(Year) + Season

#Data – seasonal averages from 2011-2021

#Plot seasonal velocity by year

ggplot(Jersey_3, aes(x=YearAdj, y=final_vel, fill=Drought))+
  geom_point(shape=21, color="black")+
  facet_wrap(~Season, scales="free")+
  drt_color_pal_drought()+
  ylab("vel (ft/s)")+
  theme_bw()

#Plot Jersey velocity by Drought Index

ggplot(Jersey_3, aes(x=Season, y=final_vel, fill=Drought))+
  geom_boxplot()+
  facet_wrap(~Season, scales="free")+
  drt_color_pal_drought()+
  ylab("vel (ft/s)")+
  theme_bw()

# ANOVA for velocity by YearAdj + Season

m_jersey_year <- aov(final_vel ~ factor(YearAdj) + Season, data = Jersey_3)

#Estimated effects may be unbalanced with factor(YearAdj)*Season
m_jersey_year

m_jersey_year_Anova <- Anova(m_jersey_year, type=2)

m_jersey_year_Anova

# check residuals

#hist(resid(m_out_year))

#fitted values of the model

data.fit = fitted(m_jersey_year)

#residuals of the model

data.res = resid(m_jersey_year)

#look at residuals

hist(data.res)
hist(data.fit)

data.stdres = rstandard(m_jersey_year)

#make qq plot and add line. first obtain normal prob plot for standardized resids

data.stdres = rstandard(m_jersey_year)

qqnorm(data.stdres)

qqline(data.stdres)

#check assumptions using Sam's function

model_plotter_jersey(m_jersey_year, Jersey_3)

#anova table

summary(m_jersey_year)

#post - hoc using Sam's tukey function

jersey_year_tukey <-tukey_plotter_jersey_yr(m_jersey_year, Jersey_3, "Season", "YearAdj")

jersey_year_tukey

ggsave(plot=jersey_year_tukey, filename=file.path(save_dir, "jersey_year_tukey.png"), device="png", height=12, width=15, units="in")

# ANOVA for Jersey velocity by Drought + Season

m_jersey_dr <- aov(final_vel ~ Drought + Season, data = Jersey_3)

summary(m_jersey_dr)

m_jersey_dr_Anova <- Anova(m_jersey_dr, type = 2)

# check residuals

#hist(resid(ex_year))

#fitted values of the model

data.fit = fitted(m_jersey_dr)

#residuals of the model

data.res = resid(m_jersey_dr)

#look at residuals

hist(data.res)
hist(data.fit)

data.stdres = rstandard(m_jersey_dr)

#make qq plot and add line. first obtain normal prob plot for standardized resids

data.stdres = rstandard(m_jersey_dr)

qqnorm(data.stdres)

qqline(data.stdres)

#check assumptions using Sam's function

model_plotter_jersey(m_jersey_dr, Jersey_3)

#anova table

summary(m_jersey_dr)

#post - hoc using Sam's tukey function

jersey_drought_tukey <-tukey_plotter_jersey_yr(m_jersey_dr, Jersey_3, "Season", "Drought")

jersey_drought_tukey

ggsave(plot=jersey_drought_tukey, filename=file.path(save_dir, "jersey_season_drought_model.png"), device="png", height=12, width=15, units="in")

#save all Anova output
anovas<-bind_rows(
  mutate(as_tibble(m_jersey_year_Anova, rownames = "Parameter"), model="Year_Season"),
  mutate(as_tibble(m_jersey_dr_Anova, rownames = "Parameter"), model="Season_Drought")
)%>%
  mutate(`Pr(>F)`=if_else(`Pr(>F)`<0.001, "< 0.001", as.character(round(`Pr(>F)`, 4))))%>%
  write_csv(file.path(save_dir, "jersey_vel_anovas.csv"))

#save all Anova output
anovas<-bind_rows(
  mutate(as_tibble(m_jersey_year_Anova, rownames = "Parameter"), model="Year_Season"),
  mutate(as_tibble(m_jersey_dr_Anova, rownames = "Parameter"), model="Season_Drought")
)

anovas$Metric <- c('jersey_vel')

anovas <- anovas %>% relocate(Metric, .before = Parameter)%>%
  mutate(`Pr(>F)`=if_else(`Pr(>F)`<0.001, "< 0.001", as.character(round(`Pr(>F)`, 4))))%>%
  write_csv(file.path(save_dir, "jersey_vel_anovas.csv"))

#compare 2021 to prior years
#load data
#raw_out <- raw_hydro_1975_2021

#raw_out$log_Outflow = log(raw_out$Outflow)

#adding drought_20_21 and yeartype_20_21 columns

jersey<-Jersey_3%>%
  filter(!is.na(final_vel))%>%
  #left_join(lt_regional%>%distinct(YearAdj, SVIndex, YearType, Drought),by="YearAdj")%>%
  mutate(across(c(Drought, YearType), list(`20_21`=~case_when(YearAdj==2021 ~ "2021",YearAdj==2020 ~ "2020",TRUE ~ as.character(.x)))),across(c(YearType, YearType_20_21), ~factor(.x, levels=c("2020", "2021", "Critical", "Dry", "Below Normal", "Above Normal", "Wet"))),Season=factor(Season, levels=c("Winter", "Spring", "Summer", "Fall")))

# graph how 2021 compares to Drought, Normal, and Wet periods?

jersey_vel<-ggplot(jersey, aes(x=Drought_20_21, y=final_vel, fill=Drought))+geom_boxplot()+drt_color_pal_drought()+xlab("Drought")+ylab("final_vel")+theme_bw()

jersey_vel

ggsave(plot=jersey_vel, filename=file.path(save_dir, "jersey_vel_drought_20_21.png"), device="png", height=4, width=5, units="in")

