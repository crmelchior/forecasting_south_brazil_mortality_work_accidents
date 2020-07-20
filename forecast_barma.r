
forecast_mortality_barma <- function (data, barma_model_order, initial_year, n_forecast = 6, diff_order = 1, verbose = TRUE) {
  if (verbose) {
    catn(paste('Using ßARMA(ar=',toString(barma_model_order$ar),', ma=',toString(barma_model_order$ma),
      ') to predict Brazilian fatal work-related accident rates',sep=''))
  }

  # Apply the difference
  data_to_forecast = data$forecast
  if (diff_order > 0) {
    data_to_forecast = diff(data_to_forecast, differences=diff_order)

    # Transform the data into double-bounded (0-1)
    a = -1
    b = 1
    data_to_forecast = (data_to_forecast-a)/(b-a)
  }

  n_data = length(data_to_forecast)

  # according to Cribari-Neto, Zeileis, 2009, "Beta Regression in R"
  extremes = FALSE
  if (min(data_to_forecast) == 0 || max(data_to_forecast) == 1) {
    extremes = TRUE
    data_to_forecast = (data_to_forecast*(n_data-1)+0.5)/n_data
  }

  # best.barma(data_to_forecast,sf=c(start=c(initial_year,1),frequency=12), pmax=6, qmax=6, X=NA , X_hat=NA)

  plot_diag = 0
  if (verbose) {
    plot_diag = 1
  }
  barma_model <- barma(data_to_forecast, ar=barma_model_order$ar, ma=barma_model_order$ma, link="logit", diag=plot_diag, h=n_forecast)

  barma_forecast <- barma_model$forecast
  barma_fitted <- barma_model$fitted

  # Undo transformations and differences
  if (extremes) {
    barma_forecast <- ((barma_forecast*n_data)-0.5)/(n_data-1)
    barma_fitted <- ((barma_fitted*n_data)-0.5)/(n_data-1)
  }
  if (diff_order > 0) {
    barma_forecast <- barma_forecast*(b-a)+a
    barma_forecast <- diffinv(barma_forecast, differences=diff_order, xi=data$forecast[length(data$forecast)])
    barma_forecast <- barma_forecast[-1]
    barma_fitted <- barma_fitted*(b-a)+a
  }
  # Fitted has max(ar,ma) NAs at the beginning
  barma_nas = max(max(barma_model_order$ar),max(barma_model_order$ma))
  if (any(is.na(barma_model_order$ar))) barma_nas = max(barma_model_order$ma)
  if (any(is.na(barma_model_order$ma))) barma_nas = max(barma_model_order$ar)
  # Ignore the NAs as the beginning
  barma_fitted <- barma_fitted[(barma_nas+1):length(barma_fitted)]
  if (diff_order > 0) {
    barma_fitted <- data$forecast[(barma_nas+1):(barma_nas+length(barma_fitted))]+barma_fitted
  }

  date_start_forecast = end(data$ts)
  date_start_forecast[2] = date_start_forecast[2] - (n_forecast-1)
  # These two ts objects are used in other scripts and should not be renamed
  ts_barma_forecast <- ts(barma_forecast, start=date_start_forecast, frequency=12)
  ts_barma_fitted <- ts(barma_fitted, start=c(initial_year, barma_nas+2), frequency=12)

  return(list(
    model=barma_model,
    forecast=ts_barma_forecast,
    fitted=ts_barma_fitted,
    resid=ts(barma_model$resid1,start=initial_year,frequency=12)
  ))
}
