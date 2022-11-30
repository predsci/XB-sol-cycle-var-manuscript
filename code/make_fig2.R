rm(list=ls())

graphics.off()

library(R.matlab)
library(ggplot2)
library(lubridate)
library(dplyr)
library(gridExtra)
library(reshape2)

source('util.R')

var_vec = c('Dcx', 'Dxt', 'Dst')

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


# =====================================================
# Plot SSN data along with mean storm size and std or total storm
# =====================================================

n1 = ncycle - 1
vec = rep(NA, n1)
stat_table = data.frame(cycle = 1:n1, dcx_mean = vec, dxt_mean = vec, dst_mean = vec, dcx_std = vec, dxt_std = vec, dst_std = vec, dcx_max = vec, dxt_max = vec, dst_max = vec, ssn_max = vec, Dcx = vec, Dxt = vec, Dst = vec)


for (ii in 1:n1) {
  start_date = cycle_date[ii]
  end_date   = cycle_date[(ii+1)]
  ssn_subset = subset(df_ssn, datetime >= start_date & datetime <= end_date)
  tmp1= tmp2 = tmp3 = tmp4 = rep(0, nvar)
  for (jj in 1:nvar) {
    dataRaw = df_storm_list[[jj]]
    data    = subset(dataRaw, evnt_time >= start_date & evnt_time <= end_date)
    if(dim(data)[1] == 0) next
    tmp1[jj] = mean(-data[,'evnt'])
    tmp2[jj] = sd(data[,'evnt'])
    tmp3[jj] = max(-data[,'evnt'])
    tmp4[jj] = length(data[,'evnt']) /(as.numeric(end_date - start_date)) # this average gives number of storms per day and we will convert to per year
    tmp4[jj] = tmp4[jj] * 365.25
  }
  stat_table[ii, c('dcx_mean', 'dxt_mean', 'dst_mean')] = tmp1
  stat_table[ii, c('dcx_std', 'dxt_std', 'dst_std')] = tmp2
  stat_table[ii, c('dcx_max', 'dxt_max', 'dst_max')] = tmp3
  stat_table[ii, 'ssn_max'] = max(ssn_subset[, 'ssn'])
  stat_table[ii, c('Dcx', 'Dxt', 'Dst')] = tmp4
}

df_stat = reshape2::melt(stat_table[,c( 'ssn_max', 'Dcx', 'Dxt', 'Dst')], id.vars = 'ssn_max')
df_stat[df_stat[,'value'] == 0, 'value' ] <- NA

pl.stat = list()

mtit = 'Average Per Year # Storms of Each Cycle Vs. Maximum SSN of Cycle '
pl.stat[[1]] <- ggplot() +
  geom_point(data = df_stat, aes(x = ssn_max, y = value, color = variable), size = 4)  +
  stat_smooth(data = df_stat, method = "lm", aes(x = ssn_max, y = value, color = variable), geom = 'smooth', show.legend = FALSE, level = 0.95, alpha = 0.4) +
  scale_y_continuous(name = 'Average # Storms per year') +
  scale_x_continuous(name = 'Maximum SSN [per cycle]', breaks = seq(from=100, to = 300, by = 25)) + #ggtitle(mtit) +
  theme(axis.text.x=element_text(size=18), axis.text.y=element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18), axis.title=element_text(size=18), legend.position = 'none') 


x = stat_table[,'ssn_max']

fit = list()

for (var in var_vec) {
  y = stat_table[,var]
  fit[[var]] = lm(y ~ x)
}

df_max = stat_table[,c( 'ssn_max', 'dcx_max', 'dxt_max', 'dst_max')]

fit2 = list()

for (var in c('dcx_max', 'dxt_max', 'dst_max')) {
  y = df_max[,var]
  fit2[[var]] = lm(y ~ x)
}

df_max[df_max[,'dst_max'] == 0, 'dst_max'] <- NA

colnames(df_max) = c('ssn_max', 'Dcx', 'Dxt', 'Dst')

df_max = reshape2::melt(df_max, id.vars = 'ssn_max')


mtit = 'Strom Max vs SSN Max'
pl.stat[[2]] <- ggplot() +
  geom_point(data = df_max, aes(x = ssn_max, y = value, color = variable), size = 4)  +
  stat_smooth(data = df_max, method = "lm", aes(x = ssn_max, y = value, color = variable), 
              geom = 'smooth', show.legend = FALSE, level = 0.95, alpha = 0.4) +
  scale_y_continuous(name = 'Maximum Storm Magnitude') +
  scale_x_continuous(name = 'Maximum SSN [per cycle]', breaks = seq(from=100, to = 300, by = 25)) + #ggtitle(mtit) +
  theme(axis.text.x=element_text(size=18), axis.text.y=element_text(size=18), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18), axis.title=element_text(size=18)) 

mp = grid.arrange(grobs = pl.stat, ncol = 2, nrow = 1)

fname = paste0(dirname(mypath),'/plots/figure2.png')
ggsave(fname, mp, width = 12, height = 5, units = 'in')

cat('\nSaved Plot at: ', fname,'\n')
# =====================================================
