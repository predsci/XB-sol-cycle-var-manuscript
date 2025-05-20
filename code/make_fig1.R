rm(list=ls())

graphics.off()

library(R.matlab)
library(ggplot2)
library(lubridate)
library(dplyr)
library(gridExtra)
library(reshape2)

# source('/Users/pete/Dropbox/research/XB-2024/util.R')
source('util.R')

# ================================
# Plot Dcx/Dxt/Dst and SSN
# ================================

var_vec = c('Dcx', 'Dxt', 'Dst')

# This is in Michal's original, but then must run Rstudio from specific DIR
#mypath = getwd()

# mypath = "/Users/pete/Dropbox/NAPR08/oulu_data/XB-sol-cycle-var-manuscript"
mypath = getwd()

data_dir = paste0(dirname(mypath),'/data/')

nvar = length(var_vec)

df_list = list()

# ================================
# Read the data
# ================================

for (ii in 1:nvar) {

  var = var_vec[ii]

  cat("Reading ", var, ' Data \n')
  if (ii < nvar) {


    # past data
    filename = paste0(data_dir,var,'.mat')

    df = read_mat_data(filename)

    # provisional data

    filename = paste0(data_dir,var,'_provisional.txt')

    df_prov = read_prov_realtime_data(filename)

    # realtime data

    filename = paste0(data_dir,var,'_realtime.txt')

    df_realt = read_prov_realtime_data(filename)

    # merge the data frames

    df_total = rbind(df, df_prov, df_realt)

  } else {

    filename = paste0(data_dir, var,'.txt')

    df = read_prov_realtime_data(filename)

    # now read data from kyoto Starting from Jan 1, 2015 todate

    filename = paste0(data_dir, var,'.kyoto.iaga2002.txt')

    df_kyoto = read_kyoto_data(filename)

    df_total = rbind(df, df_kyoto)
  }

  df_list[[ii]] = df_total
}

# =====================================================
# calculate storms
# =====================================================

evntThreshold <- -20.
evntMax       <- -100.
evntPOT       <- -20
xaxmin        <- 100

df_storm_list = list()

for (ii in 1:nvar) {
  dataRaw = df_list[[ii]]
  colnames(dataRaw) <- c('evnt_time', 'evnt')
  data = calc_storms(dataRaw = dataRaw, evntMax = evntMax, evntThreshold = evntThreshold)
  df_storm_list[[ii]] = data
}


df_dcx = df_list[[1]]
df_dxt = df_list[[2]]
df_dst = df_list[[3]]

df_storm_dcx = df_storm_list[[1]]
df_storm_dxt = df_storm_list[[2]]
df_storm_dst = df_storm_list[[3]]

# =====================================================
# read the sun-spot data
# =====================================================

filename = paste0(data_dir,'/SN_ms_tot_V2.0.txt')
sn_data <- read_sn_txt(filename)

date_start <- df_dcx[1,'datetime']

year_start <- year(date_start)
month_start <- month(date_start)

colnames(sn_data) <- c('year', 'month', 'time', 'ssn', 'std', 'measurements')
sn_data <- subset(sn_data, year >= year_start & month >= month_start)

n_ssn = length(sn_data[,1])

# =====================================================
# start year and month of each cycle and cycle number
# =====================================================

cycle_month = c(   9,    2,    4,   10,    3,    9,    8,   12,   12)
cycle_year  = c(1933, 1944, 1954, 1964, 1976, 1986, 1996, 2008, 2019)
cycle_num   = c(17,18,19,20,21,22,23,24)

ncycle = length(cycle_month)

ncycle1 = ncycle + 1

cycle_day = rep(15, ncycle)

cycle_date = ISOdate(cycle_year, cycle_month, cycle_day, hour = rep(12, ncycle), tz = 'UTC')

ssn_date = ISOdate(sn_data[,'year'], sn_data[,'month'], rep(15, n_ssn), hour = rep(12, n_ssn), tz = 'UTC')

df_ssn = data.frame(datetime = ssn_date, ssn = sn_data[,'ssn'], std = sn_data[, 'std'])

# =====================================================
# Plot SSN data along with Dcx, Dxt and Dst
# =====================================================


# The Dcx strom magnitude is x2 the SSN number - this puts them on the same scale

ssnColor <- "#69b3a2"
valColor <- 'coral' #rgb(0.2, 0.6, 0.9, 1)

coeff <- 2

pl = list()

# First for Dcx

pl[[1]] <- ggplot() +
  geom_line(data = df_ssn, aes(x = datetime, y=ssn, color = ssnColor), size =2) +
  geom_segment(data = df_storm_dcx, aes(x = evnt_time, xend = evnt_time, y=0, yend= -evnt/ coeff, color = valColor)) + # Divide by 10 to get the same range than the temperature
  scale_x_datetime(name = '') +
  scale_y_continuous(

    # Features of the first axis
    name = "SSN Monthly ",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Dcx Storms [nT]")
  )+

  theme(
    axis.title.y.left = element_text(color = valColor, size=13),
    axis.title.y.right = element_text(color = ssnColor, size=13),
    legend.position = 'none'
  ) +

  ggtitle("Monthly SSN data and Dcx Storms") +
  annotate(geom ='text', x=(cycle_date[1:8]+60*60*24*1200), y= rep(200, 8), label = cycle_num, color = 'black')

# Second for Dxt

pl[[2]] <- ggplot() +
  geom_line(data = df_ssn, aes(x = datetime, y=ssn, color = ssnColor), size =2) +
  geom_segment(data = df_storm_dxt, aes(x = evnt_time, xend = evnt_time, y=0, yend= -evnt/ coeff, color = valColor)) + # Divide by 10 to get the same range than the temperature
  scale_x_datetime(name = '') +
  scale_y_continuous(

    # Features of the first axis
    name = "SSN Monthly ",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Dxt Storms [nT]")
  )+

  theme(
    axis.title.y.left = element_text(color = valColor, size=13),
    axis.title.y.right = element_text(color = ssnColor, size=13),
    legend.position = 'none'
  ) +

  ggtitle("Monthly SSN data and Dxt Storms")

# Third for Dst

pl[[3]] <- ggplot() +
  geom_line(data = df_ssn, aes(x = datetime, y=ssn, color = ssnColor), size =2) +
  geom_segment(data = df_storm_dst, aes(x = evnt_time, xend = evnt_time, y=0, yend= -evnt/ coeff, color = valColor)) + # Divide by 10 to get the same range than the temperature
  scale_x_datetime(name = '') +
  scale_y_continuous(

    # Features of the first axis
    name = "SSN Monthly ",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Dst Storms [nT]")
  )+

  theme(
    axis.title.y.left = element_text(color = valColor, size=13),
    axis.title.y.right = element_text(color = ssnColor, size=13),
    legend.position = 'none'
  ) +

ggtitle("Monthly SSN data and Dst Storms")

mp = grid.arrange(grobs = pl,nrow = 3, ncol = 1)

fname = paste0(dirname(mypath),'/plots/figure1.png')
ggsave(fname, mp, width = 8, height = 11, units = 'in')

cat('\nSaved Plot at: ', fname,'\n')
# =====================================================

