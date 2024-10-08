% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{detdel}
\alias{detdel}
\title{Detection delay}
\usage{
detdel(
  N,
  change_type,
  t,
  k,
  alpha_plus,
  beta_plus,
  alpha_minus = 0,
  beta_minus = 0,
  max_n = 1e+05
)
}
\arguments{
\item{N}{window sizes; vector of integers}

\item{change_type}{decides the type of algortihm for the corresponding
window sizes; vector of "jump" and "kink"}

\item{t}{threshold; numerical vector}

\item{k}{number of historical data available; integer}

\item{alpha_plus}{intercept at the changepoint of the second segment;
numerical}

\item{beta_plus}{slope of the second segment; numerical}

\item{alpha_minus}{= 0; intercept at the changepoint of the first segment;
numerical}

\item{beta_minus}{= 0; slope of the first segment; numerical}

\item{max_n}{= 1e5; number of data point after which the algorithm stops if no
change is detected; integer}
}
\value{
detection delay
}
\description{
Internal function that generates data and returns the detection delay on it.
}
\keyword{internal}
