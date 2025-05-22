rm(list=ls())

library(R.matlab)
library(ggplot2)
library(poweRlaw)
library(extRemes)
library(lubridate)
library(dplyr)
library(gridExtra)

source('util.R')

mypath = getwd()

data_dir = paste0(dirname(mypath),'/data/')

# ===================================
# Read the Dcx data
# ===================================

cat("\n\n Reading Dcx Data \n\n")

var = 'Dcx'

filename = paste0(data_dir,var,'.mat')

df_dcx = read_mat_data(filename)

# provisional data

filename = paste0(data_dir,var,'_provisional.txt')

df_dcx_prov = read_prov_realtime_data(filename)

# realtime data

filename = paste0(data_dir,var,'_realtime.txt')

df_dcx_realt = read_prov_realtime_data(filename)

# merge the data frames

df_dcx_total = rbind(df_dcx, df_dcx_prov, df_dcx_realt)

# =====================================================
# start year and month of each cycle and cycle number
# =====================================================

cycle_month = c(   9,    2,    4,   10,    3,    9,    8,   12,   12)
cycle_year  = c(1933, 1944, 1954, 1964, 1976, 1986, 1996, 2008, 2019)


ncycle = length(cycle_month)

cycle_num = seq(from=(24-ncycle+2),to=24)

ncycle1 = ncycle + 1

cycle_day = rep(15, ncycle)

cycle_date = ISOdate(cycle_year, cycle_month, cycle_day, hour = rep(12, ncycle), tz = 'UTC')

# there is a total of 9 complete cycles.  We are now in the 10th one
# we will analyse one by one and also anlayze the entire Dcx data set


#  define some thresholds

evntThreshold <- -20.
evntMax       <- -100.
evntPOT       <- -20
xaxmin        <- 100

date_breaks <- "10 years"

# for standard Clauset analysis - use the same xcrit event value or if NULL uses the largest event for each cycle
# For the Carrington event used 850

inp.xcrit <- NULL # 565 # NULL
dt <- 10
yaxmax = 1
# number of bootstraps
boot = 1e4
# Number of bootstrap to show
bootShow = 10

continous <- TRUE

# =====================================================
# calculate storms for Dcx
# =====================================================

var_vec = c('Dcx')

#df_list = dfc_list = list()

dataRaw <- df_dcx_total

colnames(dataRaw) <- c('evnt_time', 'evnt')

# data frame for plotting raw data

df <- reshape::melt(dataRaw, id.var = "evnt_time")

# calculate the storms

data_storms = calc_storms(dataRaw = dataRaw, evntMax = evntMax, evntThreshold = evntThreshold)

# data frame for plotting storms

dfc <- reshape::melt(data_storms, id.var = "evnt_time")

df = data_storms

df_storm_dcx = dfc

# =====================================================
# figureS1
# =====================================================

fname = paste0(dirname(mypath),'/plots/figureS1.png')

png(file=fname, width = 800, height = 600)

par(mfrow = c(3,3))
for (ii in 1:(ncycle-1)) {

  ii1 = ii + 1

  date_start = ISOdate(cycle_year[ii], cycle_month[ii], cycle_day[ii], hour = 12, tz = 'UTC')
  date_end   = ISOdate(cycle_year[ii1], cycle_month[ii1], cycle_day[ii1], hour = 12, tz = 'UTC')
  data       = subset(df_storm_dcx, evnt_time >= date_start & evnt_time <= date_end)

  title = paste0('Cycle Number ', cycle_num[ii],' (',as.character(cycle_year[ii]),'-',as.character(cycle_year[ii1]),')')
  qqnorm(log(-data$value), cex.lab = 1.25, main = title)

}

data = df_storm_dcx
title = 'All Cycles'
qqnorm(log(-data$value), cex.lab = 1.25, main = title)

dev.off()

# =====================================================
# figureS2
# =====================================================

xmin_arr = rep(0, ncycle -1)

for (ii in 1:(ncycle-1)) {

  ii1 = ii + 1

  date_start = ISOdate(cycle_year[ii], cycle_month[ii], cycle_day[ii], hour = 12, tz = 'UTC')
  date_end   = ISOdate(cycle_year[ii1], cycle_month[ii1], cycle_day[ii1], hour = 12, tz = 'UTC')
  data       = subset(df_storm_dcx, evnt_time >= date_start & evnt_time <= date_end)

  data_estimates <- get_estimates(abs(data$value), continous = continous)
  threshold = round(data_estimates$estimates[2])

  xmin_arr[ii] = threshold

  cat("For Cycle ", ii, " KS Estimate for xmin: ", xmin_arr[ii], '\n')
}

fname = paste0(dirname(mypath),'/plots/figureS2.png')
png(file=fname, width = 800, height = 600)
par(mfrow = c(3,3), mar=c(2,4,3,2))

for (ii in 1:(ncycle-1)) {

  ii1 = ii + 1

  date_start = ISOdate(cycle_year[ii], cycle_month[ii], cycle_day[ii], hour = 12, tz = 'UTC')
  date_end   = ISOdate(cycle_year[ii1], cycle_month[ii1], cycle_day[ii1], hour = 12, tz = 'UTC')
  data       = subset(df_storm_dcx, evnt_time >= date_start & evnt_time <= date_end)
  title = paste0('Cycle Number ', cycle_num[ii],' (',as.character(cycle_year[ii]),'-',as.character(cycle_year[ii1]),') \n',
                 'KS Threshold ~ ', xmin_arr[ii])
  cat(ii, range(-data$value), '\n')
  
  rmin = 110
  rmax = round(max(-data$value)*0.7)
  nint = 40
  rmax = min(rmax, 200)
  mrl_plot <- mrlplot(-data$value, xlim = c(rmin, rmax), main = title)
  mrl_plot <- mrl_plot + geom_abline(intercept = 1, slope = 0.5, linetype = "dashed", color = "red")
  
}

data = df_storm_dcx
rmin = 110
rmax = round(max(-data$value)*0.9)
rmax = min(rmax, 200)
title = 'All Cycles'

mrlplot(-data$value, xlim = c(rmin, rmax), main = title)

dev.off()

