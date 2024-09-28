#'Sequential Jump and Kink Detection Iteration
#'
#'An Iteration of the Sequential Jump Detection and the Sequential Kink
#'Detection Algorithm
#'
#'@param new_data the new data point
#'@param detector list provided by the last iteration or the
#'SeqJaKD_initialize function
#'
#'@return list with
#'detection - True if a detection was made in this iteration
#'data - The data window needed in storage for the next iteration
#'Tstat - the test statistic
#'N - vector of integers that describes the window sizes
#'change_type - vector of strings "jump" or "kink", defines which algorithm is
#'used for the corresponding window size
#'t - numerical vector that describes the thresholds to the corresponding
#'window sizes
#'slope_est - estimate of the slope from the pre change mean
#'
#'@export
#'
#'@importFrom utils tail
SeqJaKD_iterate <- function(new_data, detector) {

  #updating the test statistics
  Tstat <- mapply(update_Tstat,
                  detector$Tstat,
                  detector$N,
                  detector$change_type,
                  MoreArgs = list(new_data = new_data,
                                  old_data = detector$data,
                                  slope_est = detector$slope_est))

  #Flagging if a change was detected
  detection <- any(abs(Tstat) >= detector$t)

  #shifting the data window
  data <- c(detector$data[-1], new_data)

  return(list(detection = detection,
              data = data,
              Tstat = Tstat,
              N = detector$N,
              change_type = detector$change_type,
              t = detector$t,
              slope_est = detector$slope_est))
}


#'Update Test statistic
#'
#'Internal function that updates the test statistic in SeqJaKD_iteration
#'
#'@param new_data new data point
#'@param old_data window of the last data points
#'@param N window size
#'@param change_type "jump" or "kink", defines the algorithm used
#'@param slope_est estimate of the slope
#'
#'@return the updatet test statistic
#'
#'@keywords internal
update_Tstat <- function(new_data,
                         old_data,
                         old_Tstat,
                         N,
                         change_type,
                         slope_est) {

  #only the last N data points are needed
  old_data <- tail(old_data,N)

  #update test statisttic differently for jump or kink
  if (change_type == "jump") {

    Tstat <- old_Tstat + (new_data - old_data[1] - slope_est * N) / N

  } else if (change_type == "kink") {

    Tstat <- old_Tstat +
            (N * new_data - sum(old_data) - slope_est * N * (N + 1) / 2) / d(N)

  } else {

    stop("change_type has to be 'jump' or 'kink'")

  }

  return(Tstat)
}


#'d_N
#'
#'helperfunction for scaling of the kink test statistic
#'
#'@param N window size
#'
#'@keywords internal
d <- function(N){
  sum((1:N)*(1:N))
}


#'Sequential Jump and Kink Detection Initialization
#'
#'Initializes the Sequential Jump and Kink Detection for SeqJaKD_iterate
#'
#'@param hist_data vector of historical data for the initialization
#'@param N vector of window sizes
#'@param change_type vector of strings "jump" or "kink" defines the type of algorithm
#'for the corresponding window size
#'@param t numerical vector of thresholds for the corresponding window sizes
#'
#'@return list with
#'detection - always set as False
#'data - The data window needed in storage for the first iteration
#'Tstat - the test statistic
#'N - vector of integers that describes the window sizes
#'change_type - vector of strings "jump" or "kink", defines which algorithm is
#'used for the corresponding window size
#'t - numerical vector that describes the thresholds to the corresponding
#'window sizes
#'slope_est - estimate of the slope of the pre change mean from the
#'historical data
#'
#'@export
#'
#'@importFrom stats coef lm predict
#'@importFrom utils tail
SeqJaKD_initialize <- function(hist_data, N, change_type, t) {

  #estimating the linear function f- on the historical data
  x <- 1:length(hist_data)
  lin <- lm(hist_data ~ x)
  slope_est <- unname(coef(lin)[2])
  pred <- unname(predict(lin,data.frame(x)))

  #initializing the test statistics with the historical data
  Tstat <- mapply(initialize_Tstat,
                  N,
                  change_type,
                  MoreArgs = list(hist_data = hist_data, pred = pred))

  return(list(detection = FALSE,
              data = tail(hist_data, max(N)),
              Tstat = Tstat,
              N = N,
              change_type = change_type,
              t = t,
              slope_est = slope_est))

}



#'Intialize test statistic
#'
#'Internal function to initialize the test statistic
#'
#'@param N window size
#'@param change_type "jump" or "kink"
#'@param hist_data vector of historical data
#'@param pred prediction of the linear function on the historical data
#'
#'@return the test statistic
#'
#'@keywords internal
initialize_Tstat <- function(N, change_type, hist_data, pred) {

  #only data in the window of size N are needed
  data <- tail(hist_data, N)
  pred <- tail(pred, N)

  #Different Test statistics for jump or kink
  if (change_type == "jump") {

    Tstat <- sum(data - pred)/N

  } else if (change_type == "kink"){

    Tstat <- sum((1:N) * (data - pred))/ d(N)

  } else {

    stop("change_type has to be 'jump' or 'kink'")

  }

  return(Tstat)
}
