
forecast_mortality_arima <- function (data, arima_model_order, n_forecast = 6, verbose = TRUE) {
  if (verbose) {
    catn(paste('Using ARIMA(',arima_model_order[1],',',arima_model_order[2],
      ',',arima_model_order[3], ') to predict Brazilian fatal work-related accident rates', sep=''))
  }

  # forecast::auto.arima(data$forecast, approximation=FALSE, allowdrift=FALSE)

  arima_model = forecast::Arima(data$forecast, order = arima_model_order)
  if (verbose) {
    print(lmtest::coeftest(arima_model))
  }

  arima_forecast = forecast::forecast(arima_model, h=n_forecast)

  if (verbose) {
    catn('ARIMA Forecasts:')
    print(arima_forecast)
  }

  return(list(
    model=arima_model,
    forecast=arima_forecast
  ))
}
