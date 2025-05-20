rm(list=ls())

graphics.off()

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

storm_table = data.frame(years = rep(0, ncycle), nstorms = rep(0, ncycle), nstorms_evd = rep(0, ncycle), threshold = rep(0, ncycle))

# =====================================================
# Create titles for plots
# =====================================================

input.var = 'Dcx'

ytit1 <- make_ytit1(input.var=input.var)

ytit11 <- make_ytit1(input.var=input.var)

# =====================================================
# Single cycle Clauset analysis for Dcx
# =====================================================

ytit2 = "P(X>xmax)"

pl.clauset = list()

icount = 0


xmin_arr = rep(0, ncycle -1)

for (ii in 1:(ncycle-1)) {

  ii1 = ii + 1

  date_start = ISOdate(cycle_year[ii], cycle_month[ii], cycle_day[ii], hour = 12, tz = 'UTC')
  date_end   = ISOdate(cycle_year[ii1], cycle_month[ii1], cycle_day[ii1], hour = 12, tz = 'UTC')
  data       = subset(df_storm_dcx, evnt_time >= date_start & evnt_time <= date_end)

  #cat('Start and End dates for cycle ', cycle_num[ii],'\n')
  # print(date_start)
  #  print(date_end)


  data_estimates <- get_estimates(abs(data$value), continous = continous)
  threshold = round(data_estimates$estimates[2])

  # order the data and set threshold as being say 20% of the data

  data_ordered <-  -data$value[order(-data$value)]
  ind <- round(length(data_ordered) * 0.5)
  xmin_arr[ii] = data_ordered[ind]
  cat(ind, length(data_ordered), '\n')
  cat("For Cycle ", ii, " KS Estimate for xmin: ", xmin_arr[ii], '\n')
}


for (ii in 1:(ncycle-1)) {

  # if (ii == (ncycle-1))
  #   next

  ii1 = ii + 1

  date_start = ISOdate(cycle_year[ii], cycle_month[ii], cycle_day[ii], hour = 12, tz = 'UTC')
  date_end   = ISOdate(cycle_year[ii1], cycle_month[ii1], cycle_day[ii1], hour = 12, tz = 'UTC')
  data       = subset(df_storm_dcx, evnt_time >= date_start & evnt_time <= date_end)


  data_estimates <- get_estimates(abs(data$value), continous = continous)
  threshold = round(data_estimates$estimates[2])

  # Do the analysis with the subjective threshold


  xmin = xmin_arr[ii]

  # Duration of data in seconds
  tSpan <-
    (as.numeric(data$evnt_time[length(data$evnt_time)]) - as.numeric(data$evnt_time[1])) /
    (86400 * 365)

  ##
  ## Clauset Analysis
  ##

  clauset.list <-
    clauset(
      data_estimates,
      boot = boot,
      bootShow = bootShow,
      xcrit = inp.xcrit,
      xmin = xmin,
      xaxmin = NULL,
      yaxmax = yaxmax,
      dt = dt,
      tSpan = tSpan,
      continous = continous
    )

  df_emp_ccdf <- clauset.list$df_emp_ccdf

  list_pl_ccdf <- clauset.list$pl_ccdf

  list_ln_ccdf <- clauset.list$ln_ccdf

  list_exp_ccdf <- clauset.list$exp_ccdf

  ## pHat_df is retrieved

  pHat_df <- clauset.list$pHat_df

  xcrit <- clauset.list$xcrit

  ##
  ## End Clauset Analysis
  ##

  title_clauset = paste0('Cycle ',cycle_num[ii],' (',as.character(cycle_year[ii]),'-',as.character(cycle_year[ii1]),')')

  ## Plots of the traces

  xaxmin = min(min(abs(df_emp_ccdf$x)), xmin)

  yaxmin = 1e-3
  yaxmax = 1

  icount = icount + 1
  pl.clauset[[icount]] <-
    ggplot() + scale_y_continuous(
      name = ytit2,
      trans = "log",
      breaks = c(1e-3, 1e-2, 0.1, 1),
    ) + scale_x_continuous(
      name = ytit1,
      trans = "log",
      breaks = seq(from=100, to = 575, by = 50)) +
    ggtitle(title_clauset) +
    coord_cartesian(xlim = c(99, 580), ylim=c(yaxmin, yaxmax))

  for (i in 1:bootShow) {

    boot_exp_ccdf = list_exp_ccdf[[i]]
    boot_ln_ccdf  = list_ln_ccdf[[i]]
    boot_pl_ccdf  = list_pl_ccdf[[i]]

    pl.clauset[[icount]] <-
      pl.clauset[[icount]] + geom_line(data = boot_exp_ccdf,
                                       aes(
                                         x = value,
                                         y = p_x,
                                         colour = "Exponential"), alpha = 0.4)

    pl.clauset[[icount]] <-
      pl.clauset[[icount]] + geom_line(data = boot_ln_ccdf,
                                       aes(
                                         x = value,
                                         y = p_x,
                                         colour = "Log-Normal"), alpha = 0.4)
    pl.clauset[[icount]] <-
      pl.clauset[[icount]] + geom_line(data = boot_pl_ccdf,
                                       aes(
                                         x = value,
                                         y = p_x,
                                         colour = "Power-Law"), alpha = 0.4)

    pl.clauset[[icount]] <-
      pl.clauset[[icount]] + geom_point(
        data = df_emp_ccdf,
        aes(x = x, y = y, colour = "Observations"),
        color = "black",
        alpha = 0.5,
        size = 2
      )

    if (i == 1)
      pl.clauset[[icount]] <-
      pl.clauset[[icount]] + scale_colour_manual(values = c(
        `Power-Law` = "red",
        `Log-Normal` = "blue",
        Exponential = "green"
      )) + scale_alpha_manual(values = c(
        `Power-Law` = 0.2,
        `Log-Normal` = 0.2,
        Exponential = 0.2
      ))
    pl.clauset[[icount]] <-
      pl.clauset[[icount]] + theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.direction = "horizontal"
      )
    pl.clauset[[icount]] <-
      pl.clauset[[icount]] + geom_vline(
        data = NULL,
        xintercept = xmin,
        color = "grey2",
        size = 1,
        linetype = "dashed"
      )

  }

  ## The density plot

  icount = icount + 1
  pl.clauset[[icount]] <-
    ggplot() + geom_density(
      data = pHat_df,
      aes(
        x = prob,
        stat(scaled),
        colour = model,
        fill = model
      ),
      adjust = 1.5,
      alpha = 0.2
    ) + scale_x_continuous(name = ytit2,
                           limits = c(0, 1)) + scale_y_continuous(name = "Density") + scale_colour_manual(values = c(
                             `Power-Law` = "red",
                             `Log-Normal` = "blue",
                             Exponential = "green"
                           )) + scale_alpha_manual(values = c(
                             `Power-Law` = 0.2,
                             `Log-Normal` = 0.2,
                             Exponential = 0.2
                           )) + scale_fill_manual(values = c(
                             `Power-Law` = "red",
                             `Log-Normal` = "blue",
                             Exponential = "green"
                           )) + theme(legend.position = "none",
                                      legend.title = element_blank()) + ggtitle(title_clauset)
}



# This is a big plot with lots of panels


mp = grid.arrange(grobs = pl.clauset,
                  nrow = 4,
                  ncol = 4)

# fname = paste0(dirname(mypath),'/figures/keep_half_data_figure3.png')
fname = paste0('../plots/figureS3.png')
ggsave(fname, mp, width = 30, height = 20, units = 'in')

cat('\nSaved Plot at: ', fname,'\n')
# =====================================================

