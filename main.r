#!/usr/bin/env Rscript

#rm(list=ls()) # Cleaning the memory

initial_year = 2000
state_code = 1 # 1=PR, 2=SC, 3=RS
n_forecast = 6 # observations (months) to forecast to validate the model
verbose = TRUE # if it should print results

# ARIMA, ßARMA, and KARMA models
if (state_code == 1) { # PR
  diff_order = 1
  arima_model_order = c(1,1,1)
  barma_model_order = list(ar=NA, ma=1:2)
  karma_model_order = list(ar=1:4, ma=1)
}
if (state_code == 2) { # SC
  diff_order = 1
  arima_model_order = c(0,1,1)
  barma_model_order = list(ar=1:2,ma=2:3)
  karma_model_order = list(ar=1:2,ma=1:4)
}
if (state_code == 3) { # RS
  diff_order = 0
  arima_model_order = c(1,1,2)
  barma_model_order = list(ar=1:3,ma=1:4)
  karma_model_order = list(ar=1:3,ma=1:3)
}

source('utils.r')

title = paste('Rate of fatal work-related accidents -', get_state_acronym(state_code))
title_twolines = paste('Rate of fatal work-related\naccidents -', get_state_acronym(state_code))

# Loads data
data = load_mortality_data(state_code, initial_year, n_forecast)
data$diff = diff(data$forecast) # First difference

catn(title)

# Uncomment the next line to plot all results into a PDF file
# pdf(paste('plots',get_state_acronym(state_code),'series.pdf', sep='_'))

# Margins order: bottom, left, top, right
original_parameters = par() # save the original margins for plotting

###########################
##   Visual inspection   ##
###########################
catn('Visual inspection')

# svg(paste('plot_original',get_state_acronym(state_code),'series.svg', sep='_'), width=10, height=2.5, pointsize=12)
# suppressMessages(library(forecast)) # Needed to ggplot2 autoplot support ts (Hyndman and Khandakar, 2008)
# ggplot2::autoplot(data$forecast, xlab = "Years", ylab = title)
par(mar=c(3.5, 4, 1, 2), mgp=c(1.8,0.6,0)) # Reduce margins
plot(data$forecast, xlab = "Years", ylab = title_twolines)
par(mar=original_parameters$mar, mgp=original_parameters$mgp) # Restore margins
# tryCatch(dev.off(), error=function(e){}) # close the SVG file if it was opened

# Decomposing time series
plot(stats::stl(data$forecast, "periodic"))
title("STL in level")

plot(stats::stl(data$diff, "periodic"))
title("STL first difference")

plot(decompose(data$forecast))
title("Decompose in level")

plot(decompose(data$diff))
title("Decompose first difference")

# ACF and PACF plots

acf_data <- forecast::Acf(data$forecast, plot=FALSE)
pacf_data <- forecast::Pacf(data$forecast, plot=FALSE)


# svg(paste('plot_acf_pacf',get_state_acronym(state_code),'series.svg', sep='_'), width=10, height=2.3, pointsize=12)
par(mfrow=c(1,2),mar=c(3, 3, 1, 1), mgp=c(1.5,0.6,0))
# We do not plot the ACF with lag=0 as it is always 1
plot(acf_data[2:length(acf_data$acf)-1], main='', ylab=paste('ACF -', get_state_acronym(state_code)))
plot(pacf_data, main='', ylab=paste('PACF -', get_state_acronym(state_code)))
par(mfrow=original_parameters$mfrow, mar=original_parameters$mar, mgp=original_parameters$mgp) # Restore margins
# tryCatch(dev.off(), error=function(e){}) # close the SVG file if it was opened


###########################
##    Unit root tests    ##
###########################
catn(''); # line break
catn('Performing unit root tests')

# Lag order is calculated using (long) Schwert rule
# Schwert, C. (1989). Tests for Unit Roots: A Monte Carlo Investigation. Journal of Business & Economic Statistics, 7(2), 147-159. doi:10.1080/07350015.1989.10509723
trunc_lag_param = trunc(12*((length(data$forecast)/100)^(1/4))) # Schwert rule

# ADF - Augmented Dickey-Fuller Test
# H0: The series have a unit root (non-stationary)
# H1: The series does not have a unit root (stationary)
tseries::adf.test(data$forecast, k = trunc_lag_param)
tseries::adf.test(data$diff, k = trunc_lag_param)

# KPSS test
# H0: The series is stationary in level
# H1: The series is not stationary in level
tseries::kpss.test(data$forecast, null="Level", lshort=FALSE)
tseries::kpss.test(data$diff, null="Level", lshort=FALSE)


###########################
##      Forecasting      ##
###########################
catn(''); # line break
catn('Forecasting phase')

### ARIMA ###
source('forecast_arima.r')
# Forecasts ARIMA and prints results
mortality_arima = forecast_mortality_arima(data, arima_model_order, n_forecast, verbose)

catn(''); # line break
catn(''); # line break

### ßARMA ###
source('forecast_barma.r')

# ßARMA sources available in:
# http://www.ufsm.br/bayer/boot-barma.zip
# https://github.com/vscher/barma
source('barma/barma.fit.r')
source('barma/barma.r')

# Forecasts ßARMA and prints results
mortality_barma = forecast_mortality_barma(data, barma_model_order, initial_year, n_forecast, diff_order, verbose)


### KARMA ###
source('forecast_karma.r')

# KARMA sources available in:
# https://github.com/fabiobayer/KARMA
source('karma/kum-mu-phi.r')
source('karma/karma.fit.r')
source('karma/karma.r')

# Forecasts KARMA and prints results
mortality_karma = forecast_mortality_karma(data, karma_model_order, initial_year, n_forecast, diff_order, verbose)


###########################
##     Forecast plots    ##
###########################

### Out-sample plot ###

# svg(paste('plot_out_sample',get_state_acronym(state_code),'series.svg', sep='_'), width=10, height=3, pointsize=12)

line_legend = c("original", "ARIMA", "ßARMA", "KARMA")
line_colors = c("black","red","blue","green")
legend_text_size = 0.8
initial_year_plot = 2016
out_sample_plot_data <- load_mortality_data(state_code, initial_year_plot, n_forecast)
n_plot_data = length(out_sample_plot_data$ts)
end_year_plot = end(out_sample_plot_data$ts)[1]
out_sample_ylim = c(0.04, 0.12)

# We take the last value before forecast and put together with the forecast data to join the lines
last_data_before_forecasts = as.data.frame(out_sample_plot_data$ts)[(n_plot_data-n_forecast):(n_plot_data-n_forecast),1]
ts_arima_forecast_plot <- ts(
  c(last_data_before_forecasts,mortality_arima$forecast$mean),
  start=c(end_year_plot, (12-n_forecast)), frequency=12)
ts_barma_forecast_plot <- ts(
  c(last_data_before_forecasts, mortality_barma$forecast),
  start=c(end_year_plot, (12-n_forecast)), frequency=12)
ts_karma_forecast_plot <- ts(
  c(last_data_before_forecasts, mortality_karma$forecast),
  start=c(end_year_plot, (12-n_forecast)), frequency=12)

tss_to_plot = ts.union(out_sample_plot_data$ts, ts_arima_forecast_plot, ts_barma_forecast_plot, ts_karma_forecast_plot)

par(mar=c(2.5, 4.2, 1, 1), mgp=c(2,0.6,0)) # Reduce margins
plot(tss_to_plot,
  plot.type = "single",
  ylab = title_twolines,
  xlab = "",
  main = "",
  lty = 1:length(line_colors), col = line_colors,
  axes = FALSE, # we'll make the axis ourselves
  ylim = out_sample_ylim
)
legend("bottomleft",
  legend = line_legend,
  lty = 1:length(line_colors), col = line_colors,
  bty = "n", # "o" to put a box surrounding legend
  cex = legend_text_size, # size of the legend text
  text.col = "black",
  horiz = F,
  inset = c(0, -0.02)
)
axis(2) # Y axis
# Add ticks for each 6-month period
for (year_tick in initial_year_plot:end_year_plot) {
  axis(1, labels=paste("jan",year_tick,sep="/"), at=year_tick)
  axis(1, labels=paste("jul",year_tick,sep="/"), at=(year_tick+0.5))
}
axis(1, labels=paste("dec",end_year_plot,sep="/"), at=(end_year_plot+0.9))
box() # box surrounding plot
par(mar=original_parameters$mar, mgp=original_parameters$mgp) # Restore margins

# tryCatch(dev.off(), error=function(e){}) # close the SVG file if it was opened


### In-sample plot ###

# svg(paste('plot_in_sample',get_state_acronym(state_code),'series.svg', sep='_'), width=10, height=3.5, pointsize=12)

tss_to_plot = ts.union(data$ts, mortality_arima$model$fitted, mortality_barma$fitted, mortality_karma$fitted)
initial_year_plot = start(data$ts)[1]
end_year_plot = end(data$ts)[1]

par(mar=c(2.5, 3.4, 1, 1), mgp=c(2,0.6,0)) # Reduce margins
plot(tss_to_plot,
  plot.type = "single",
  ylab = title,
  xlab = "",
  main = "",
  lty = 1:length(line_colors), col = line_colors,
  axes = FALSE # we'll make the axes ourselves
)
legend("bottomleft", 
  legend = line_legend,
  lty = 1:length(line_colors), col = line_colors,
  bty = "n", # "o" to put a box surrounding legend
  cex = legend_text_size, # size of the legend text
  text.col = "black",
  horiz = F,
  inset = c(0, -0.02)
)
axis(2) # Y axis
# Add ticks each two years
for (year_tick in seq(initial_year_plot, end_year_plot, by=2)) {
  axis(1, labels=year_tick, at=year_tick)
}
box() # box surrounding plot
par(mar=original_parameters$mar, mgp=original_parameters$mgp) # Restore margins

# tryCatch(dev.off(), error=function(e){}) # close the SVG file if it was opened


###########################
## Accuracy measurements ##
###########################
catn(''); # line break
catn('Accuracy measurements')

arima_accuracy = forecast::accuracy(mortality_arima$forecast$mean, data$validate)
barma_accuracy = forecast::accuracy(mortality_barma$forecast, data$validate)
karma_accuracy = forecast::accuracy(mortality_karma$forecast, data$validate)

catn('ARIMA accuracy measurements:')
print(arima_accuracy)

catn('ßARMA accuracy measurements:')
print(barma_accuracy)

catn('KARMA accuracy measurements:')
print(karma_accuracy)

pct_decimal = 1
idx_rmse = 2
idx_mae = 3
idx_mape = 5

catn('ßARMA vs ARIMA:')
catn(paste('RMSE',round(((barma_accuracy[idx_rmse]/arima_accuracy[idx_rmse])-1)*100,pct_decimal),'%'))
catn(paste('MAE',round(((barma_accuracy[idx_mae]/arima_accuracy[idx_mae])-1)*100,pct_decimal),'%'))
catn(paste('MAPE',round(((barma_accuracy[idx_mape]/arima_accuracy[idx_mape])-1)*100,pct_decimal),'%'))
catn('KARMA vs ARIMA:')
catn(paste('RMSE',round(((karma_accuracy[idx_rmse]/arima_accuracy[idx_rmse])-1)*100,pct_decimal),'%'))
catn(paste('MAE',round(((karma_accuracy[idx_mae]/arima_accuracy[idx_mae])-1)*100,pct_decimal),'%'))
catn(paste('MAPE',round(((karma_accuracy[idx_mape]/arima_accuracy[idx_mape])-1)*100,pct_decimal),'%'))
catn('ßARMA vs KARMA:')
catn(paste('RMSE',round(((barma_accuracy[idx_rmse]/karma_accuracy[idx_rmse])-1)*100,pct_decimal),'%'))
catn(paste('MAE',round(((barma_accuracy[idx_mae]/karma_accuracy[idx_mae])-1)*100,pct_decimal),'%'))
catn(paste('MAPE',round(((barma_accuracy[idx_mape]/karma_accuracy[idx_mape])-1)*100,pct_decimal),'%'))


###########################
##    Normality tests    ##
###########################
catn(''); # line break

catn('ARIMA Normality tests')
# Shapiro-Wilk test
# H0: data came from a normally distributed population
shapiro.test(residuals(mortality_arima$model))
# Lilliefors (Kolmogorov-Smirnov) Test For Normality
# H0: data came from a normally distributed population
nortest::lillie.test(residuals(mortality_arima$model))
# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(residuals(mortality_arima$model))

catn('ßARMA Normality tests')
# Shapiro-Wilk test
# H0: data came from a normally distributed population
shapiro.test(mortality_barma$resid)
# Lilliefors (Kolmogorov-Smirnov) Test For Normality
# H0: data came from a normally distributed population
nortest::lillie.test(mortality_barma$resid)
# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(mortality_barma$resid)

catn('KARMA Normality tests')
# Shapiro-Wilk test
# H0: data came from a normally distributed population
shapiro.test(mortality_karma$resid)
# Lilliefors (Kolmogorov-Smirnov) Test For Normality
# H0: data came from a normally distributed population
nortest::lillie.test(mortality_karma$resid)
# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(mortality_karma$resid)

###########################
##      Portmanteau      ##
###########################
catn(''); # line break
catn('Portmanteau tests')

# Sources available in: https://github.com/vscher/barma
source("vscher_barma/dufour.r")
source("vscher_barma/Kwan_Chest.R")

initial_df = 5
n_portmanteau_tests = 6
m_tests = seq(from=5, to=5+n_portmanteau_tests)

# ARIMA portmanteau tests
arima_fitdf = arima_model_order[1] + arima_model_order[3] # p + q
arima_ljungbox = matrix(rep(NA,n_portmanteau_tests), nrow=1, ncol=n_portmanteau_tests)
arima_dufour = matrix(rep(NA,n_portmanteau_tests), nrow=1, ncol=n_portmanteau_tests)
for (m in seq(1,n_portmanteau_tests)) {
  df = initial_df-1+m+arima_fitdf
  arima_ljungbox[,m] = Box.test(residuals(mortality_arima$model), lag=df, fitdf=arima_fitdf, type="Ljung-Box")$p.value
  arima_dufour[,m] = Dufour.test(residuals(mortality_arima$model), lag=df, fitdf=arima_fitdf)$p.value
  # vscher test is specific for BARMA #TODO
}

# ßARMA portmanteau tests
barma_fitdf = max(barma_model_order$ar) + max(barma_model_order$ma) # p + q
if (any(is.na(barma_model_order$ar))) barma_fitdf = max(barma_model_order$ma)
if (any(is.na(barma_model_order$ma))) barma_fitdf = max(barma_model_order$ar)
barma_ljungbox = matrix(rep(NA,n_portmanteau_tests), nrow=1, ncol=n_portmanteau_tests)
barma_dufour = matrix(rep(NA,n_portmanteau_tests), nrow=1, ncol=n_portmanteau_tests)
barma_vscher = matrix(rep(NA,n_portmanteau_tests), nrow=1, ncol=n_portmanteau_tests)
for (m in seq(1,n_portmanteau_tests)) {
  df = initial_df-1+m+barma_fitdf
  barma_ljungbox[,m] = Box.test(mortality_barma$resid, lag=df, fitdf=barma_fitdf, type="Ljung-Box")$p.value
  barma_dufour[,m] = Dufour.test(mortality_barma$resid, lag=df, fitdf=barma_fitdf)$p.value
  barma_vscher[,m] = Kwan.sim.chest(mortality_barma$resid, lag=df, fitdf=barma_fitdf,type="partial",test=4)$p.value
}

# KARMA portmanteau tests
karma_fitdf = max(karma_model_order$ar) + max(karma_model_order$ma) # p + q
if (any(is.na(karma_model_order$ar))) karma_fitdf = max(karma_model_order$ma)
if (any(is.na(karma_model_order$ma))) karma_fitdf = max(karma_model_order$ar)
karma_ljungbox = matrix(rep(NA,n_portmanteau_tests), nrow=1, ncol=n_portmanteau_tests)
karma_dufour = matrix(rep(NA,n_portmanteau_tests), nrow=1, ncol=n_portmanteau_tests)
for (m in seq(1,n_portmanteau_tests)) {
  df = initial_df-1+m+karma_fitdf
  karma_ljungbox[,m] = Box.test(mortality_karma$resid, lag=df, fitdf=karma_fitdf, type="Ljung-Box")$p.value
  karma_dufour[,m] = Dufour.test(mortality_karma$resid, lag=df, fitdf=karma_fitdf)$p.value
}

# Plotting portmanteau results
for (model in c('arima', 'barma','karma')) {
  for (test in c('ljungbox', 'dufour', 'vscher')) {
    results_name = paste(model, test, sep='_')
    if (exists(results_name)) {
      fitdf = get(paste(model, 'fitdf', sep='_'))
      results_data = get(results_name)

      plot((fitdf+1):(fitdf+n_portmanteau_tests), results_data, xlab = "lag",ylab = "p-value",las=1, ylim = c(0,1),lwd=1.5)
      abline(h = 0.05, lty = 2, col = "blue", lwd=1.5) # 95% confidence interval
      text((fitdf+1):(fitdf+n_portmanteau_tests), results_data, round(results_data,4), pos=3, cex=0.7)
      title(paste(get_state_acronym(state_code), results_name))
    }
  }
}

tryCatch(dev.off(), error=function(e){}) # close the PDF file if it was opened
