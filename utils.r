
# state_code is 1=PR, 2=SC, 3=RS
get_state_acronym <- function(state_code) { return(c('PR', 'SC', 'RS')[state_code]) }

catn <- function(x, append='\n'){ cat(x); cat(append) }

load_mortality_data <- function (state_code, initial_year, n_forecast = 6) {
    original_data = read.csv('data.csv')

    filtered_data = original_data[original_data$'state'==state_code,]
    filtered_data = filtered_data[filtered_data$'year'>=initial_year,]
    mortality_rate = filtered_data['rate']
    ts_rate = ts(mortality_rate, start=c(initial_year,1), frequency=12)

    df_rate = as.data.frame(ts_rate)
    
    n_rates = dim(df_rate)[1] # sample size (number of entries that we have on df_rate = 204)
    n = n_rates-n_forecast # remove the last k entries
    data_to_forecast = ts(df_rate[1:n,1],start=c(initial_year,1),frequency=12)
    data_to_validate = c(df_rate[(n+1):n_rates,1]) # take the k last entries
    
    return(list(
        raw=mortality_rate[['rate']],
        ts=ts_rate,
        df=df_rate,
        forecast=data_to_forecast,
        validate=data_to_validate
    ))
}
