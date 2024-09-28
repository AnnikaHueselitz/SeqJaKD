#'threshold tuning
#'
#'Tunes the threshold for given window sizes to bound the type I error such
#'that the probability of a false alarm before tau is approximately eta.
#'
#'@param N vector of window sizes
#'@param change_type vector of "jump" and "kink" deciding the algorithm used for
#'each window size
#'@param tau integer that defines to which time the algorithm has to run with
#'a probability of eta without a false detection
#'@param eta numerical in (0,1) that is the probability with which the algorithm
#'runs until eta without a false detection
#'@param k integer that decides the number of available historical data
#'@param r integer that decides the number of maximal test statistic that is
#'drawn to decide the threshold
#'
#'@return vector of tuned thresholds
#'
#'@export
#'
#'@importFrom stats lm predict rnorm
t_tuning <- function(N, change_type, tau, eta, k, r) {

  #calculate the maximal test statistics r-times
  M <- replicate(r, max_Tstat(N, change_type, tau, k))

  #format and sort the maximal test statistics
  M <- matrix(M,ncol = r)
  M_sorted <- t(apply(M,1,sort))

  #calculate the position of the quantile
  quantile <- floor(r*eta)

  #calculate the start for checking for the quantile
  q <- r-quantile

  #Checks if the desired quantile is fulfilled
  while(M_quantile(M,M_sorted[,q]) > quantile){
    q <- q + 1
  }

  #returns the quantile of the maximal test statistics as the threshold
  return(M_sorted[,q])
}

#'maximal test statistic
#'
#'Internal function that calculates the maximal test statistic
#'
#'@param N vector of integers; window sizes
#'@param change_type vector of "jump" and "kink"; type of algorithm
#'@param tau integer; number of test statistics to calculate
#'@param k integer; amount of historical data
#'
#'@return maximal test statistic
#'
#'@keywords internal

max_Tstat <- function(N, change_type, tau, k){

  #historical data for the test statistics
  data <- rnorm(k)

  #predicts the linear function on the data
  x <- 1:k
  lin <- lm(data ~ x)
  pred <- unname(predict(lin,data.frame(x = 1:(k+tau))))

  #add data of length tau
  data <- c(data,rnorm(tau))
  data <- data - pred


  #calculate the max_Tstat for all window sizes N
  M <- mapply(max_Tstat_N,
              N,
              change_type,
              MoreArgs = list(data = data, tau = tau))

  return(M)
}


#'Maximal test statistic quantile
#'
#'Internal function returns the quantile of the threshold and
#' maximal test statistics
#'
#'@param M matrix of maximal test statistics
#'@param t possible threshold
#'
#'@return the quantile that is reached with this threshold
#'
#'@keywords internal
M_quantile <- function(M, t){

  #sum of max test statistics over the threshold
  quantile <- sum(apply(M, 2, function(x){any(x > t)}))

  return(quantile)
}

#'maximal test statistic for single window
#'
#'calculates the maximal test statistic for a single window
#'
#'@param N window size; integer
#'@param change_type type of algorithm; "jump" or "kink"
#'@param data vector of data
#'@param tau number of test statistics to calculate; integer
#'
#'@return maximal test statistic
#'
#'@keywords internal
max_Tstat_N <- function(N, change_type, data, tau){

  #initialize the test statistics as the data
  Tstat <- data

  #different test statistic for jump and kink
  if (change_type == "jump") {

    #summing the data in windows of size N
    for (i in 2:N) {

      data <- data[-1]
      Tstat <- Tstat[-length(Tstat)]
      Tstat <- Tstat + data

    }

    #scaling
    Tstat <- Tstat / N

  } else if (change_type == "kink") {

    #summing and weighting the data in windows of size N
    for(i in 2:N){

      data <- data[-1]
      Tstat <- Tstat[-length(Tstat)]
      Tstat <- Tstat + (i * data)

    }

    #scaling
    Tstat <- Tstat / d(N)

  } else {

    stop("change_type has to be 'jump' or 'kink'")

  }

  #only the last tau are test statistics
  Tstat <- tail(Tstat,tau)

  return(max(abs(Tstat)))
}
