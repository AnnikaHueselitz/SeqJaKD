#'Estimate the expected detection delay
#'
#'Estimate the expected detection delay for the Sequential Jump and Kink
#'Detection algorithm.
#'
#'@param N window sizes; vector of integers
#'@param change_type decides the type of algortihm for the corresponding
#'window sizes; vector of "jump" and "kink"
#'@param t threshold; numerical vector
#'@param k number of historical data available; integer
#'@param alpha_plus intercept at the changepoint of the second segment;
#'numerical
#'@param beta_plus slope of the second segment; numerical
#'@param alpha_minus = 0; intercept at the changepoint of the first segment;
#'numerical
#'@param beta_minus = 0; slope of the first segment; numerical
#'@param max_n = 1e5; number of data point after which the algorithm stops if no
#'change is detected; integer
#'@param m = 200; number of repetitions; integer
#'
#'@return estimation of expected detection delay
#'
#'@export
#'
#'@importFrom stats rnorm
estimate_detdel <- function(N,
                            change_type,
                            t,
                            k,
                            alpha_plus,
                            beta_plus,
                            alpha_minus = 0,
                            beta_minus = 0,
                            max_n = 1e5,
                            m = 200) {

  #detection delay over 200 repetitions
  detdel <- replicate(m, detdel(N,
                                change_type,
                                t,
                                k,
                                alpha_minus,
                                alpha_plus,
                                beta_minus,
                                beta_plus,
                                max_n))

  return(mean(detdel))
}

#'Detection delay
#'
#'Internal function that generates data and returns the detection delay on it.
#'
#'@param N window sizes; vector of integers
#'@param change_type decides the type of algortihm for the corresponding
#'window sizes; vector of "jump" and "kink"
#'@param t threshold; numerical vector
#'@param k number of historical data available; integer
#'@param alpha_plus intercept at the changepoint of the second segment;
#'numerical
#'@param beta_plus slope of the second segment; numerical
#'@param alpha_minus = 0; intercept at the changepoint of the first segment;
#'numerical
#'@param beta_minus = 0; slope of the first segment; numerical
#'@param max_n = 1e5; number of data point after which the algorithm stops if no
#'change is detected; integer
#'
#'@return detection delay
#'
#'@keywords internal
detdel <- function(N,
                   change_type,
                   t,
                   k,
                   alpha_plus,
                   beta_plus,
                   alpha_minus = 0,
                   beta_minus = 0,
                   max_n = 1e5) {

  #generate historical data for initializing
  hist_data <- rnorm(k) + seg_lin_fun(1:k,k,alpha_minus, beta_minus,alpha_plus, beta_plus)

  #initialize detector to iterate over
  detector <- SeqJaKD_initialize(hist_data, N, change_type, t)

  #initialize time
  n <- 0

  #iterate until a detection was made
  while (!detector$detection) {

    n <- n+1

    #stop if detection takes too long
    if (n > max_n) {

      stop("max. number of iterations exceeded")

    }

    #generate data at time n
    new_data <- rnorm(1)+ seg_lin_fun(k + n, k ,
                             alpha_minus, beta_minus,
                             alpha_plus, beta_plus)

    detector <- SeqJaKD_iterate(new_data, detector)

  }

  #return time of detection
  return(n)

}


#'Estimate false alarm probability
#'
#'Estimates the false alarm probability of the Sequential Jump and Kink Detector,
#'that is the probability of a false alarm befor tau
#'
#'@param N window sizes; vector of integers
#'@param change_type decides the type of algortihm for the corresponding
#'window sizes; vector of "jump" and "kink"
#'@param t threshold; numerical vector
#'@param k number of historical data available; integer
#'@param tau estimates the probaility of a false alarm before tau; integer
#'@param m number of repetitions; integer
#'
#'@return estimated false alarm probability
#'
#'@export
#'
#'@importFrom stats rnorm
estimate_faprob <- function(N, change_type,t ,k ,tau, m = 200) {

  #check for 200 repetitions if an false alarm occurs before tau
  false_alarm <- replicate(m, fa(N, change_type, t, k, tau))

  return(sum(false_alarm)/m)
}


#'False alarm
#'
#'Internal functiion that generates data and reports if a false alarm was made
#'by the Sequential Jump and Kink Detector before tau observations.
#'
#'@param N window sizes; vector of integers
#'@param change_type decides the type of algorithm for the corresponding
#'window sizes; vector of "jump" and "kink"
#'@param t threshold; numerical vector
#'@param k number of historical data available; integer
#'@param tau reports if a false alarm occurs before tau observations; integer
#'
#'@return TRUE if a false alarm occurred
#'
#'@keywords internal
fa <- function(N, change_type, t, k, tau){

  #initialize detector
  detector <- SeqJaKD_initialize(rnorm(k), N, change_type, t)

  n <- 1

  #iterate until a false detection was made
  while (!detector$detection) {

    n <- n+1

    #Return no false detection, if at time tau no false detection was made
    if (n > tau) {

      return(FALSE)

    }

    detector <- SeqJaKD_iterate(rnorm(1), detector)
  }

  #return false detection
  return(TRUE)
}

#'segmented linear function
#'
#'A function that returns a segmented linear function at times x
#'
#'@param x points to evaluate the function at; vector of integers
#'@param tau changepoint; integer
#'@param alpha_minus intercept at tau of first segment; numeric
#'@param beta_minus slope on first segment; numeric
#'@param alpha_plus intercept at tau of second segment; numeric
#'@param beta_plus slope on second segment; numeric
#'
#'@return vector that describes a segmented linear function
#'
#'@export
seg_lin_fun <- function(x,
                        tau,
                        alpha_minus,
                        beta_minus,
                        alpha_plus,
                        beta_plus) {

  #first segment
  f1 <- (x[x <= tau] - tau) * beta_minus + alpha_minus

  #second segment
  f2 <- (x[x > tau] - tau) * beta_plus + alpha_plus

  return(c(f1,f2))
}

