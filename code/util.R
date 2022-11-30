read_dcx_mat_data <- function(filename) {
  
  var_data = readMat(filename)
  
  #str(var_data)
  
  # order is year, month, day, hour, val 
  
  year  = var_data$Dcx[,1]
  
  month = var_data$Dcx[,2]
  
  day   = var_data$Dcx[,3]
  
  hour  = var_data$Dcx[,4]
  
  val   = var_data$Dcx[,5]
  
  # create an object of class POSIXct
  
  dates = ISOdate(year, month, day, hour, tz = 'UTC')
  
  df = data.frame(datetime = dates, val = val)
  
  return(df)
}

read_kyoto_data <- function(filename) {
  
  var_data = read.table(filename, skip = 18, header = FALSE)
  
  colnames(var_data) = c('date', 'time', 'doy', 'dst')
  
  val = var_data[,'dst']
  
  var_data[,'date'] = as.Date(var_data[,'date'])
  
  year  = year( var_data[,'date'])
  
  month = month(var_data[,'date'])
  
  day   = day(  var_data[,'date'])
  
  hour = as.numeric(substr(var_data[,'time'], start = 1, stop = 2))
  
  dates = ISOdate(year, month, day, hour, tz = 'UTC')
  
  df = data.frame(datetime = dates, val = val) 
  
}

read_prov_realtime_data <- function(filename) {
  
  var_data = read.table(filename, header = FALSE)
  
  year  = var_data[,2]
  
  month = var_data[,3]
  
  day   = var_data[,4] 
  
  hour  = var_data[,5]
  
  val   = var_data[,6]
  
  dates = ISOdate(year, month, day, hour, tz = 'UTC')
  
  df = data.frame(datetime = dates, val = val)  
  
  return (df)
  
}

read_mat_data <- function(filename) {
  
  var_data = readMat(filename)
  
  data = var_data[[1]]
  
  # order is year, month, day, hour, val 
  
  year  = data[,1]
  
  month = data[,2]
  
  day   = data[,3]
  
  hour  = data[,4]
  
  val   = data[,5]
  
  # create an object of class POSIXct
  
  dates = ISOdate(year, month, day, hour, tz = 'UTC')
  
  df = data.frame(datetime = dates, val = val)
  
  return(df)
}

read_sn_txt <- function (filename) {
  
  data = read.table(filename)
  
}

calc_storms <- function(dataRaw, evntMax = NULL, evntThreshold = NULL) {
  
  if (is.null(evntThreshold)) {
    evntThreshold = -20
  }
  
  if (is.null(evntMax)) {
    evntMax = -100
  }
  
  if (evntThreshold < 0) {
    icrit = (dataRaw$evnt <= evntThreshold)
  } else {
    icrit = (dataRaw$evnt >= evntThreshold)
  }
  
  blocks = rle(icrit)
  iLoc = blocks$length
  iStorm = blocks$values
  lenBlock = length(iLoc)
  numStorm = length(iStorm[(iStorm == T)])
  
  dstStorm = replicate(numStorm, 0)
  dstStormTime <- .POSIXct(double(numStorm), tz="UTC") #.POSIXct(character(numStorm))
  iStormIndex = 1
  
  # Identify all intervals and tag with the peak value and time of that peak
  # for negative/positive thresholds use min/max functions.
  
  if (evntThreshold < 0) {
    for (i in 1:lenBlock) {
      if (iStorm[i] == T) {
        x1 = sum(iLoc[1:(i - 1)])
        x2 = x1 + iLoc[i]
        dstStorm[iStormIndex] = min(dataRaw$evnt[x1:x2])
        dstStormTime[iStormIndex] = dataRaw$evnt_time[x1]
        iStormIndex = iStormIndex + 1
      }
    }
    
  } else {
    for (i in 1:lenBlock) {
      if (iStorm[i] == T) {
        x1 = sum(iLoc[1:(i - 1)])
        x2 = x1 + iLoc[i]
        dstStorm[iStormIndex] = max(dataRaw$evnt[x1:x2])
        dstStormTime[iStormIndex] = dataRaw$evnt_time[x1]
        iStormIndex = iStormIndex + 1
      }
      
    }
    
  }
  
  dstStruct = data.frame(evnt_time = dstStormTime, evnt = dstStorm)
  
  # Only keep those values where the peak exceeded evntMax
  if (evntMax < 0) {
    dstStruct = dplyr::filter(dstStruct, evnt < evntMax)
  } else {
    dstStruct = dplyr::filter(dstStruct, evnt > evntMax)
  }
  
  return(dstStruct)
  
}


get_estimates <- function( emp_data, continous = TRUE) {
  
  emp_data <- emp_data[order(emp_data)]
  
  if ( length(which(emp_data == 0)) != 0 ) {
    emp_data <- emp_data[-1*which(emp_data == 0)]
  }
  
  if (continous == F) {     
    pl_m   = displ$new(emp_data)
  } else {
    pl_m   = conpl$new(emp_data)
  }
  
  est_pl = estimate_xmin(pl_m, distance = 'ks')
  
  KS = est_pl$gof
  xmin = est_pl$xmin
  alpha = est_pl$pars
  ntail = est_pl$ntail
  ntotal <- length(emp_data)
  
  estimates <- list( estimates = c(KS, xmin, alpha, ntail, ntotal), data_ordered = emp_data)
  return (estimates)
  
}


dst_extRemes <- function(data = NULL, threshold = 100) {
  
  blocks = rle(year(data$evnt_time))
  iEvent = blocks$length
  iYear = blocks$values
  
  # Need average number of events per year 
  # Exclude the last year from the average because it is likley to be incomplete and less accurate 
  
  events.per.year = round(mean(iEvent[1:(length(iYear) - 1)]))
  
  time.units = paste0(events.per.year, "/year")
  
  fitGP <- fevd(abs(data$evnt), threshold = threshold, type = 'GP', time.units = time.units, units = 'nT', verbose = TRUE)
  
  fitEXP <- fevd(abs(data$evnt), threshold = threshold, type = 'Exponential', time.units = time.units, units = 'nT', verbose = FALSE)
  summary(fitGP)
  return(fitGP)
  
}

clauset <- function(data, boot = 100, bootShow = 10, xcrit = NULL, xmin = NULL, xaxmin = NULL, yaxmax = NULL, continous = T, dt = NULL, tSpan = NULL) {
  
  estimates <- data[[1]]
  emp_data <- data[[2]]
  emp_data <- emp_data[order(emp_data)]
  
  if (missing(xcrit) | is.null(xcrit)) {
    xcrit <- max(emp_data, na.rm = TRUE)
  }
  
  if (missing(xmin) | is.null(xmin)) {
    xmin <- estimates[2] # could set this to a hard value of 100 to test the model fits 
  }
  
  emp_ccdf <- data.frame(value = unique(emp_data), p_x = NA)
  n_emp_data <- length(emp_data)
  for (i in 1:nrow(emp_ccdf)) {
    emp_ccdf$p_x[i] <- length(which(emp_data >= emp_ccdf$value[i]))/n_emp_data
  }
  
  if (length(which(emp_ccdf$value == 0)) != 0) {
    emp_ccdf <- emp_ccdf[-1 * which(emp_ccdf$value == 0), ]
  }
  
  
  if (missing(yaxmax) | is.null(yaxmax)) {
    yaxmax = max(emp_ccdf$p_x)
  }
  
  
  df_emp_ccdf = data.frame(x = emp_ccdf$value, y = emp_ccdf$p_x)
  
  # these are the probabilities for an event as large as the largest event in the dataset
  qLogNorm = replicate(boot, NA)
  qExp = replicate(boot, NA)
  qPowLaw = replicate(boot, NA)
  
  # these are the probabilities for an extreme event or larger
  qXcritLogNorm = replicate(boot, NA)
  qXcritExp = replicate(boot, NA)
  qXcritPowLaw = replicate(boot, NA)
  
  list_ln_ccdf = list_pl_ccdf = list_exp_ccdf = list()
  icount = 0
  
  for (i in 1:boot) {
    
    # to avoid including NaN do a while loop here
    
    jcount = 0
    while (jcount == 0) {
      y <- emp_data[sample(length(emp_data), length(emp_data), replace = TRUE)]
      y <- y[order(y)]
      y_tail <- y[which(y >= xmin)]
      y_tail <- y_tail[order(y_tail)]
      
      newYTail = c(emp_ccdf$value)
      newYTail = newYTail[which(newYTail >= xmin)]
      
      p_tail <- length(y_tail)/length(y)
      
      # compute the log-normal traces #######################
      
      if (continous == F) {
        m_ln <- dislnorm$new(y_tail)
      } else {
        m_ln <- conlnorm$new(y_tail)
      }
      m_ln$setXmin(xmin)
      m_ln$mle()
      
      ln_ccdf <- data.frame(value = newYTail, p_x = NA)
      ln_ccdf$p_x <- 1 - dist_cdf(m_ln, ln_ccdf$value)
      ln_ccdf$p_x <- ln_ccdf$p_x * p_tail
      
      # The probability of an event as large as the biggest event in this dataset
      
      if (max(ln_ccdf$value) == (max(emp_data))) {
        qLogNorm[i] = min(ln_ccdf$p_x)
      }
      # compute the probability for an xcrit event
      qXcritLogNorm[i] = p_tail * (1 - dist_cdf(m_ln, xcrit))
      if (!is.nan(qXcritLogNorm[i] )) jcount = 1
      
    }
    
    # compute the exponential traces ##############
    
    if (continous == F) {
      m_exp <- disexp$new(y_tail)
    } else {
      m_exp <- conexp$new(y_tail)
    }
    m_exp$setXmin(xmin)
    m_exp$mle()
    
    exp_ccdf <- data.frame(value = newYTail, p_x = NA)
    exp_ccdf$p_x <- 1 - dist_cdf(m_exp, exp_ccdf$value)
    exp_ccdf$p_x <- exp_ccdf$p_x * p_tail
    if (max(exp_ccdf$value) == (max(emp_data))) {
      qExp[i] = min(exp_ccdf$p_x)
    }
    
    # compute the probability for an xcrit event
    qXcritExp[i] = p_tail * (1 - dist_cdf(m_exp, xcrit))
    
    
    # compute and plot the power-law traces ##############
    
    pl_ccdf <- data.frame(value = newYTail, p_x = NA)
    
    sum <- 0
    for (value in y_tail) {
      sum <- sum + log(value/(xmin - (1/2)))
    }
    
    alpha <- 1 + length(y_tail) * ((sum)^(-1))
    
    pl_ccdf$p_x <- (pl_ccdf$value/xmin)^(-1 * alpha + 1)
    pl_ccdf$p_x <- pl_ccdf$p_x * p_tail
    
    if (max(pl_ccdf$value) == (max(emp_data))) {
      qPowLaw[i] = min(pl_ccdf$p_x)
    }
    
    # compute the probability for an xcrit event
    qXcritPowLaw[i] = p_tail * (xcrit/xmin)^(-1 * alpha + 1)
    
    if ((i%%(boot/bootShow)) == 0) {
      icount = icount + 1
      list_exp_ccdf[[icount]] = exp_ccdf
      list_pl_ccdf[[icount]]  = pl_ccdf
      list_ln_ccdf[[icount]]  = ln_ccdf
    }
    
  } # the end of the boot do-loop. 
  
  
  ##### Plot the density plot of the extreme values (maximum observed) ###########
  
  nEvents = length(emp_data)
  
  pHatLogNorm = 1. - exp(-nEvents*qLogNorm)
  pHatExp     = 1. - exp(-nEvents*qExp)
  pHatPowLaw  = 1. - exp(-nEvents*qPowLaw)
  
  #Compute probabilities that a particular event will occur over next X years
  pLogNormXyears = 1. - exp(-nEvents*dt*qXcritLogNorm/tSpan)
  pExpXyears     = 1. - exp(-nEvents*dt*qXcritExp/tSpan)
  pPowLawXyears  = 1. - exp(-nEvents*dt*qXcritPowLaw/tSpan)
  
  ci95_pl = quantile(pPowLawXyears , c(0.025, 0.975),na.rm=T)
  ci95_ln = quantile(pLogNormXyears, c(0.025, 0.975),na.rm=T)
  ci95_se = quantile(pExpXyears    , c(0.025, 0.975),na.rm=T)
  
  c2 = c(ci95_pl[1], ci95_ln[1], ci95_se[1]) * 100.0
  c3 = c(ci95_pl[2], ci95_ln[2], ci95_se[2]) * 100.0
  c1 = c(median(pPowLawXyears, na.rm = TRUE), median(pLogNormXyears, na.rm = TRUE), median(pExpXyears, na.rm = TRUE)) * 100.0
  
  c1 = signif(c1, digits = 3)
  c2 = signif(c2, digits = 3)
  c3 = signif(c3, digits = 3)
  
  pEventXyears_df <- data.frame(Model = c('Power-Law', 'Log-Normal', 'Exponential'), c1, c2, c3)
  
  pHat_df = data.frame(model = factor(rep(c("Power-Law", "Log-Normal", "Exponential"), each = boot)), prob = c(pHatPowLaw, pHatLogNorm, pHatExp))
  
  
  return(list(df_emp_ccdf = df_emp_ccdf, pl_ccdf = list_pl_ccdf, 
              ln_ccdf = list_ln_ccdf, exp_ccdf = list_exp_ccdf, 
              pHat_df = pHat_df, pEventXyears_df = pEventXyears_df, xcrit = xcrit,
              qXcritLogNorm = qXcritLogNorm, qXcritExp = qXcritExp, qXcritPowLaw = qXcritPowLaw, 
              fac = nEvents*dt/tSpan))
  
}

getrlpoints <- function(fit) {
  
  xp2 <- ppoints(fit$n, a = 0)
  ytmp <- datagrabber(fit)
  y <- c(ytmp[, 1])
  sdat <- sort(y)
  npy <- fit$npy
  u <- fit$threshold
  rlpoints.x <- -1/log(xp2)[sdat > u]/npy
  rlpoints.y <- sdat[sdat > u]
  rlpoints <- data.frame(rlpoints.x, rlpoints.y)
  
  return(rlpoints)
}

## Return a data frame of estimated returned levels with confidence intervals

getcidf <- function(fit, alpha = NULL, rperiods = NULL) {
  if (is.null(alpha)) alpha = 0.05
  if (is.null(rperiods)) rperiods = c(2, 5, 10, 20, 50, 100,200, 500)
  bds <- ci.fevd.mle(fit, return.period = rperiods, alpha = 0.5)
  c1 <- as.numeric(bds[, 1])
  c2 <- as.numeric(bds[, 2])
  c3 <- as.numeric(bds[, 3])
  ci_df <- data.frame(c1, c2, c3, rperiods)
  colnames(ci_df) <- c('2.5%CI', 'Median', '97.5%CI', 'Return-Period')
  return(ci_df)
}


compDistributions <- function(data,  continous = T) {
  
  # now compare different distributions. 
  
  if (continous == F) {
    pl_m = displ$new(data)
    ln_m = dislnorm$new(data)
    exp_m = disexp$new(data)
  } else {
    pl_m = conpl$new(data)
    ln_m = conlnorm$new(data)
    exp_m = conexp$new(data)
  }
  
  est_pl = estimate_xmin(pl_m)
  pl_m$setXmin(est_pl)
  
  # need to set xmin the same, so must get from previous estimate of PL
  
  ln_m$setXmin(pl_m$getXmin())
  est_ln = estimate_pars(ln_m)
  ln_m$setPars(est_ln$pars)
  
  exp_m$setXmin(pl_m$getXmin())
  est_exp = estimate_pars(exp_m)
  exp_m$setPars(est_exp)
  
  comp_pl_ln = compare_distributions(pl_m, ln_m)
  comp_pl_exp = compare_distributions(pl_m, exp_m)
  comp_ln_exp = compare_distributions(ln_m, exp_m)
  return(list(comp_pl_ln = comp_pl_ln, comp_pl_exp = comp_pl_exp, comp_ln_exp = comp_ln_exp))
}

make_ytit1 <- function(input.var='dst') {
  if (input.var == 'dst') {
    ytit1 = "-Dst [nT]"
  } else if (input.var == 'ae') {
    ytit1 = 'AE [nT]'
  } else if (input.var == 'al') {
    ytit1 = '-AL [nT]'
  } else if (input.var == 'au') {
    ytit1 = 'AU [nT]'
  } else if (input.var == 'ap') {
    ytit1 = 'AP [nT]'
  } else if (input.var == 'ey') {
    ytit1 = 'Ey Electric Field mV/m'
  } else if (input.var == 'bz') {
    ytit1 = '-Bz (GSM) [nT]'
  } else if (input.var == 'sw') {
    ytit1 = 'Solar Wind Speed [km/sec]'
  } else if (input.var == 'beta') {
    ytit1 = 'Plasma Beta'
  } else if (input.var == 'sunspot') {
    ytit1 = 'Sunspot Number [new version]'
  } else if (input.var == 'f107') {
    ytit1 = 'Solar Index F10.7'
  } else if (input.var == 'cme') {
    ytit1 = 'CME Speed [km/sec]'
  } else if (input.var == 'bn') {
    ytit1 = 'IMF |Bn| [nT]'
  }else {
    ytit1 = 'Magnitude'
  }
  return(ytit1)
  
}

make_ytit11 <-function(input.var='dst') {
  if (input.var == 'dst') {
    ytit11 = "-Dst > 100 [nT] "
  } else if (input.var == 'ae') {
    ytit11 = 'AE > 1000 [nT]' 
  } else if (input.var == 'al') {
    ytit11 = '-AL > 600 [nT]' 
  } else if (input.var == 'au') {
    ytit11 = 'AU > 400 [nT]' 
  } else if (input.var == 'ap') {
    ytit11 = 'AP > 100 [nT]' 
  } else if (input.var == 'ey') {
    ytit11 = 'Ey Electric Field > 8 [mV/m]' 
  } else if (input.var == 'bz') {
    ytit11 = '-Bz > 10 [nT]' 
  } else if (input.var == 'sw') {
    ytit11 = 'Solar Wind Speed > 600 [km/sec]' 
  } else if (input.var == 'beta') {
    ytit11 = 'Plasma Beta > 20' 
  } else if (input.var == 'sunspot') {
    ytit11 = 'Sunspot Number (new version) > 150' 
  } else if (input.var == 'f107') {
    ytit11 = 'Solar Index F10.7 > 100' 
  } else if (input.var == 'cme') {
    ytit11 = 'CME Speed > 700 [km/sec]'
  } else if (input.var == 'bn') {
    ytit11 = 'IMF |Bn| > 10 [nT]'
  }else {
    ytit11 = 'Magnitude > Median Magnitude'
  }
  return(ytit11)
  
}

make_ytit11 <-function(input.var='dst') {
  if (input.var == 'dst') {
    ytit11 = "-Dst > 100 [nT] "
  } else if (input.var == 'ae') {
    ytit11 = 'AE > 1000 [nT]' 
  } else if (input.var == 'al') {
    ytit11 = '-AL > 600 [nT]' 
  } else if (input.var == 'au') {
    ytit11 = 'AU > 400 [nT]' 
  } else if (input.var == 'ap') {
    ytit11 = 'AP > 100 [nT]' 
  } else if (input.var == 'ey') {
    ytit11 = 'Ey Electric Field > 8 [mV/m]' 
  } else if (input.var == 'bz') {
    ytit11 = '-Bz > 10 [nT]' 
  } else if (input.var == 'sw') {
    ytit11 = 'Solar Wind Speed > 600 [km/sec]' 
  } else if (input.var == 'beta') {
    ytit11 = 'Plasma Beta > 20' 
  } else if (input.var == 'sunspot') {
    ytit11 = 'Sunspot Number (new version) > 150' 
  } else if (input.var == 'f107') {
    ytit11 = 'Solar Index F10.7 > 100' 
  } else if (input.var == 'cme') {
    ytit11 = 'CME Speed > 700 [km/sec]'
  } else if (input.var == 'bn') {
    ytit11 = 'IMF |Bn| > 10 [nT]'
  }else {
    ytit11 = 'Magnitude > Median Magnitude'
  }
  return(ytit11)
  
}



