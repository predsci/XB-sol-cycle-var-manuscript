rm(list=ls())
graphics.off()

library(R.matlab)
library(ggplot2)
library(lubridate)
library(dplyr)
library(poweRlaw)
library(extRemes)

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


continous <- TRUE

input.var = "Dcx"

# calculate the 10 year probability of exceeding each of these events for each cycle and for all 

xcrit_vec =  c(400, 500, 600) 

ncrit = length(xcrit_vec)


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


##########################################
# calculate return period using one cycle data
##########################################


returnlevel_df = NULL 
tmp_df = data.frame(xcrit = xcrit_vec, return_period = rep(0, ncrit),year = rep(0, ncrit), ssn= rep(0, ncrit))
icount = 0

for (ii in 1:(ncycle-1)) {

	ii1 = ii + 1
	if (ii == ncycle) {
		data = df_storm_dcx
		date_start = data[1, "evnt_time"]
		date_end = max(stdata[, "evnt_time"])
		data_ssn = df_ssn
		cycle_years = "ALL"
	} else {
		data = df_storm_dcx
		date_start = ISOdate(cycle_year[ii], cycle_month[ii], cycle_day[ii], hour = 12, tz = "UTC")
		date_end = ISOdate(cycle_year[ii1], cycle_month[ii1], cycle_day[ii1], hour = 12, tz = "UTC")
		data = subset(df_storm_dcx, evnt_time >= date_start & evnt_time <= date_end)
		data_ssn = subset(df_ssn, datetime >= date_start & datetime <= date_end)
		cycle_years = paste0(as.character(cycle_year[ii]), "-", as.character(cycle_year[ii1]))
		
	}

	print("Start and End dates")
	print(date_start)
	print(date_end)

	# data estimates

	data_estimates <- get_estimates(abs(data$evnt), continous = continous)

	threshold = round(data_estimates$estimates[2])

	# =====================================================
	# POT analysis - returns the return levels with confidence intervals
  # =====================================================

	fitGP <- dst_extRemes(data = data, threshold = threshold)

	rlpoints <- getrlpoints(fitGP)

	ci_df <- getcidf(fitGP, alpha = 0.5, rperiods = 2:1000)
  
	# now for this cycle interpolate find the return periods for the events of magnitude we are interested in
	for (jj in 1:ncrit) {
		xcrit = xcrit_vec[jj]
		val = approx(x = ci_df[,'Median'], y=ci_df[,'Return-Period'], xout = xcrit)$y
		tmp_df[jj, 'return_period'] = round(val, digits = 2) #ci_df[ind,'Return-Period']
		tmp_df[jj, 'year'] = cycle_years
		tmp_df[jj, 'ssn']  = max(data_ssn$ssn)
	}
	
	returnlevel_df  = rbind(returnlevel_df, tmp_df)

}

returnlevel_df[, 'xcrit'] = as.character(returnlevel_df[, 'xcrit'])

# =====================================================
# plot the return period vs max ssn
# =====================================================

mtit = paste0('Return Period vs. Max SSN of each cycle')
pl.rl.ssn.logy = ggplot(returnlevel_df) + geom_point(aes(x = ssn, y = return_period, color = xcrit),  size = 4) + 
  scale_y_continuous(name = 'Return Period [years]', trans = "log10") +
  scale_x_continuous(name = 'Maximum SSN [per cycle]', breaks = seq(from=100, to = 300, by = 25)) + ggtitle(mtit) + 
  guides(color = guide_legend(title = "Xcrit [nT]"),size = guide_legend(14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18), axis.title=element_text(size=14)) 


fname = paste0(dirname(mypath),'/plots/figure7.png')
ggsave(fname, pl.rl.ssn.logy, width = 12, height = 7, units = 'in')

cat('\nSaved Plot at: ', fname,'\n')
# =====================================================

