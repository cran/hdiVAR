% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CV_VARMLE.R
\name{CV_VARMLE}
\alias{CV_VARMLE}
\title{cross-validation for transition matrix update in maximization step}
\usage{
CV_VARMLE(tol_seq, ht_seq, S0_train, S1_train, Y_test, is_echo = FALSE)
}
\arguments{
\item{tol_seq}{vector; grid of tolerance parameter in Dantzig selector for cross-validation.}

\item{ht_seq}{vector; grid of hard-thresholding levels for transition matrix estimate.
To avoid hard thresholding, set \code{ht_seq=0}.}

\item{S0_train}{a p by p matrix; average (over time points in training data) of conditional expectation of \eqn{x_t x_t^\top} on \eqn{y_1, \ldots, y_T} and parameter estimates, obtained from expectation step.}

\item{S1_train}{a p by p matrix; average (over time points in training data) of conditional expectation of \eqn{x_t x_{t+1}^\top}on \eqn{y_1, \ldots, y_T} and parameter estimates, obtained from expectation step.}

\item{Y_test}{a p by T_test matrix; observations of time series in test set.}

\item{is_echo}{logical; if true, display the information of CV-optimal (tol, ht).}
}
\value{
a list of CV-optimal parameters and test prediction error.
\tabular{ll}{
\code{tol_min}  \tab  CV-optimal tolerance parameter in Dantzig selector. \cr
\code{ht_min}  \tab  CV-optimal hard thresholding level for the output of Dantzig selector. \cr
\code{test_loss}  \tab  a matrix of prediction error in test data; columns match \code{tol_seq}, and rows match \code{ht_seq}. \cr
}
}
\description{
Tune the tolerance parameter of generalized Dantzig selector and hard thresholding
level via prediction error in test data.
}
\author{
Xiang Lyu, Jian Kang, Lexin Li
}
