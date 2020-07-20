
forecast_mortality_karma <- function (data, karma_model_order, initial_year, n_forecast = 6, diff_order = 1, verbose = TRUE) {
  if (verbose) {
    catn(paste('Using KARMA(ar=',toString(karma_model_order$ar),', ma=',toString(karma_model_order$ma),
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

  # best.karma(data_to_forecast,sf=c(start=c(initial_year,1),frequency=12), pmax=6, qmax=6, X=NA , X_hat=NA)

  plot_diag = 0
  if (verbose) {
    plot_diag = 1
  }
  karma_model <- karma(data_to_forecast, ar=karma_model_order$ar, ma=karma_model_order$ma, link="logit", diag=plot_diag, h=n_forecast)

  karma_forecast <- karma_model$forecast
  karma_fitted <- karma_model$fitted

  # Undo transformations and differences
  if (extremes) {
    karma_forecast <- ((karma_forecast*n_data)-0.5)/(n_data-1)
    karma_fitted <- ((karma_fitted*n_data)-0.5)/(n_data-1)
  }
  if (diff_order > 0) {
    karma_forecast <- karma_forecast*(b-a)+a
    karma_forecast <- diffinv(karma_forecast, differences=diff_order, xi=data$forecast[length(data$forecast)])
    karma_forecast <- karma_forecast[-1]
    karma_fitted <- karma_fitted*(b-a)+a
  }
  # Fitted has max(ar,ma) NAs at the beginning
  karma_nas = max(max(karma_model_order$ar),max(karma_model_order$ma))
  if (any(is.na(karma_model_order$ar))) karma_nas = max(karma_model_order$ma)
  if (any(is.na(karma_model_order$ma))) karma_nas = max(karma_model_order$ar)
  # Ignore the NAs as the beginning
  karma_fitted <- karma_fitted[(karma_nas+1):length(karma_fitted)]
  if (diff_order > 0) {
    karma_fitted <- data$forecast[(karma_nas+1):(karma_nas+length(karma_fitted))]+karma_fitted
  }

  date_start_forecast = end(data$ts)
  date_start_forecast[2] = date_start_forecast[2] - (n_forecast-1)
  # These two ts objects are used in other scripts and should not be renamed
  ts_karma_forecast <- ts(karma_forecast, start=date_start_forecast, frequency=12)
  ts_karma_fitted <- ts(karma_fitted, start=c(initial_year, karma_nas+2), frequency=12)

  return(list(
    model=karma_model,
    forecast=ts_karma_forecast,
    fitted=ts_karma_fitted,
    resid=ts(karma_model$resid1,start=initial_year,frequency=12)
  ))
}
