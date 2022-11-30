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

xcrit_vec = c(400, 500, 600) 

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


##########################################
# Plot all in one panel
#########################################

pEventXyears_df_all = NULL
for (jj in 1:ncrit) {
  pEventXyears_df_all = rbind(pEventXyears_df_all, pEventXyears_df_list[[jj]])
}

mtit = paste0('10 Year Probability of an Event Exceeding ', xcrit_vec[1],' ,',xcrit_vec[2],' ,',  xcrit_vec[3],' [nT]' )

pEventXyears_df_all$x = rep(c(0, 1:(ncycle-1)), ncrit)
pEventXyears_df_all$xcrit = as.character(rep(xcrit_vec, each = ncycle))

# labels for x-axis
x_labels = unique(pEventXyears_df_all$years) 


pl.all = ggplot(pEventXyears_df_all) + geom_point(aes(x = x, y = median/100., color = xcrit),  size = 4) + 
  geom_errorbar(aes(x = x, ymin = low/100., ymax = high/100., color = xcrit), alpha = 0.8, width = 0.1) +
  scale_y_continuous(name = 'Probability') +
  scale_x_continuous(name = 'Years', breaks = 0:(ncycle-1), labels = x_labels) + ggtitle(mtit) +
  guides(color = guide_legend(title = "Xcrit [nT]"),size = guide_legend(14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14), axis.title.y=element_text(size=18), axis.title=element_text(size=14))


fname = paste0(dirname(mypath),'/plots/figure5.png')
ggsave(fname, pl.all, width = 12, height = 7, units = 'in')

cat('\nSaved Plot at: ', fname,'\n')
# =====================================================


