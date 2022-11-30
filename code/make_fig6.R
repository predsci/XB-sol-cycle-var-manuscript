rm(list=ls())
graphics.off()

library(R.matlab)
library(ggplot2)
library(lubridate)
library(poweRlaw)
library(extRemes)
library(dplyr)
library(gridExtra)

source('util.R')

# ===================================
# create the 10 -year prob of exceeding event 
# ===================================

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
# calculate storms 
# =====================================================

evntThreshold <- -20.
evntMax       <- -100.
evntPOT       <- -20
xaxmin        <- 100

dataRaw = df_dcx_total
colnames(dataRaw) <- c('evnt_time', 'evnt')
data = calc_storms(dataRaw = dataRaw, evntMax = evntMax, evntThreshold = evntThreshold)
df_storm = data

df_dcx = df_dcx_total

df_storm_dcx =df_storm 

# =====================================================
# start year and month of each cycle and cycle number
# =====================================================

cycle_month = c(   9,    2,    4,   10,    3,    9,    8,   12,   12)
cycle_year  = c(1933, 1944, 1954, 1964, 1976, 1986, 1996, 2008, 2019)

ncycle = length(cycle_month)

cycle_num = seq(from=(24-ncycle+2),to=24)

ncycle = length(cycle_month)

ncycle1 = ncycle + 1

cycle_day = rep(15, ncycle)

cycle_date = ISOdate(cycle_year, cycle_month, cycle_day, hour = rep(12, ncycle), tz = 'UTC')

# there is a total of 9 complete cycles.  We are now in the 10th one 
# we will analyse one by one and also anlayse the entire Dcx data set 

date_breaks <- "10 years"

inp.xcrit <- NULL
dt <- 10
yaxmax = 1
boot = 500
bootShow = 10

continous <- TRUE

input.var = "Dcx"

# calculate the 10 year probability of exceeding each of these events for each cycle and for all 

xcrit_vec =  c(400, 500, 600) 

ncrit = length(xcrit_vec)

pEventXyears_df = data.frame(years = rep(0, ncycle), median = rep(0, ncycle), low = rep(0, ncycle), high = rep(0, ncycle))

pEventXyears_df_list = list()

icount = 0

for (jj in 1:ncrit) {

	inp.xcrit <- xcrit_vec[jj]
	
	cat("Processing Xcrit of ", inp.xcrit, ' [nT]\n')
	for (ii in 1:ncycle) {

	  ii1m = ii - 1
		ii1  = ii + 1
		
		if (ii == 1) {
		  data = df_storm_dcx
		  date_start = data[1, "evnt_time"]
		  date_end = max(data[, "evnt_time"])
		  cycle_years = "ALL" 
		  cat('Start and End dates for cycle ALL\n')
		} else {
		  data = df_storm_dcx
		  date_start = ISOdate(cycle_year[ii1m] ,  cycle_month[ii1m],  cycle_day[ii1m], hour = 12, tz = "UTC")
		  date_end   = ISOdate(cycle_year[ii], cycle_month[ii], cycle_day[ii], hour = 12, tz = "UTC")
		  data = subset(df_storm_dcx, evnt_time >= date_start & evnt_time <= date_end)
		  cycle_years = paste0(as.character(cycle_year[ii1m]), "-", as.character(cycle_year[ii]))
		  cat('Start and End dates for cycle ', cycle_num[ii1m],'\n')
		}
		
		print(date_start)
		print(date_end)

		data_estimates <- get_estimates(abs(data$evnt), continous = continous)
		threshold = round(data_estimates$estimates[2])

		xmin = threshold

		tSpan <- (as.numeric(data$evnt_time[length(data$evnt_time)]) - as.numeric(data$evnt_time[1]))/(86400 * 365)

		##
		## Clauset Analysis
		##

		clauset.list <- clauset(data_estimates, boot = boot, bootShow = bootShow, xcrit = inp.xcrit, xmin = NULL, xaxmin = NULL, yaxmax = yaxmax, dt = dt, tSpan = tSpan, continous = continous)

		xcrit <- clauset.list$xcrit

		pEventXyears <- clauset.list$pEventXyears

		pEventXyears_df[ii, 2:4] = as.numeric(pEventXyears[1, 2:4])
		pEventXyears_df[ii, "years"] = cycle_years

	}
	
	pEventXyears_df_list[[jj]] = pEventXyears_df
}


# unpack add xcrit column and pack as one big data frame

pEventXyears_df_all = NULL
for (jj in 1:ncrit) {
  pEventXyears_df = pEventXyears_df_list[[jj]]
  pEventXyears_df$xcrit = rep(xcrit_vec[jj], ncycle)
  pEventXyears_df_all = rbind(pEventXyears_df_all, pEventXyears_df)
}


# ####################################################################################################################
# Append to this data frame the maximum SSN for each cycle, number of storms and average storm size
# ####################################################################################################################


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

ssn_date = ISOdate(sn_data[,'year'], sn_data[,'month'], rep(15, n_ssn), hour = rep(12, n_ssn), tz = 'UTC')

df_ssn = data.frame(datetime = ssn_date, ssn = sn_data[,'ssn'], std = sn_data[, 'std'])


n1 = ncycle - 1
vec = rep(NA, n1)
stat_table = data.frame(cycle = 1:n1, dcx_mean = vec, dxt_mean = vec, dst_mean = vec, dcx_std = vec, dxt_std = vec, dst_std = vec, dcx_total = vec, dxt_total = vec, dst_total = vec, ssn_max = vec, ssn_sum = vec, ssn_mean = vec)


for (ii in 1:n1) {
  start_date = cycle_date[ii]
  end_date   = cycle_date[(ii+1)]
  ssn_subset = subset(df_ssn, datetime >= start_date & datetime <= end_date)

  dataRaw = df_storm_dcx
  data    = subset(dataRaw, evnt_time >= start_date & evnt_time <= end_date)
  if(dim(data)[1] == 0) next
    tmp1 = mean(-data[,'evnt'])
    tmp2 = sd(data[,'evnt'])
    tmp3 = sum(-data[,'evnt'])
    
  stat_table[ii, c('dcx_mean', 'dxt_mean', 'dst_mean')] = tmp1
  stat_table[ii, c('dcx_std', 'dxt_std', 'dst_std')] = tmp2
  stat_table[ii, c('dcx_total', 'dxt_total', 'dst_total')] = tmp3
  stat_table[ii, 'ssn_max'] = max(ssn_subset[, 'ssn'])
  stat_table[ii, 'ssn_sum'] = sum(ssn_subset[, 'ssn'])
  stat_table[ii, 'ssn_mean'] = mean(ssn_subset[, 'ssn'])
}


# we need to remove the 'ALL' from the data frame before we add the ssn data 
# will give it a different name

pEventXyears_ssn_df = subset(pEventXyears_df_all , years != 'ALL')

pEventXyears_ssn_df$xcrit   = rep(as.character(xcrit_vec), each = n1)
pEventXyears_ssn_df$ssn_sum = rep(stat_table[, 'ssn_sum'], ncrit)
pEventXyears_ssn_df$ssn_max = rep(stat_table[, 'ssn_max'], ncrit)
pEventXyears_ssn_df$ssn_mean = rep(stat_table[, 'ssn_mean'], ncrit)


mtit = 'Probability of Event Vs. Maximum SSN'

pl.ssn = ggplot(pEventXyears_ssn_df) + geom_point(aes(x = ssn_max, y = median/100., color = xcrit),  size = 4) + 
  geom_errorbar(aes(x = ssn_max, ymin = low/100., ymax = high/100., color = xcrit), alpha = 0.8, width = 0.2) +
  scale_y_continuous(name = 'Probability') +
  scale_x_continuous(name = 'Maximum SSN [per cycle]', breaks = seq(from=100, to = 300, by = 25)) + ggtitle(mtit) + 
  guides(color = guide_legend(title = "Xcrit [nT]"),size = guide_legend(14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14), axis.title.y=element_text(size=18), axis.title=element_text(size=14)) +
  geom_smooth(data = pEventXyears_ssn_df, method = 'lm', formula=(y~x), aes(x= ssn_max, y = median/100, color = xcrit), stat = 'smooth', span = 1)

fname = paste0(dirname(mypath),'/plots/figure6.png')
ggsave(fname, pl.ssn, width = 12, height = 7, units = 'in')

cat('\nSaved Plot at: ', fname,'\n')

#  =====================================================


