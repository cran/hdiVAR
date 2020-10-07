#' statistical inference for transition matrix in high-dimensional vector autoregression with measurement error
#'
#' @description Conduct global and simultaneous testing on the transition matrix.
#'
#' @param Y observations of time series, a p by T matrix.
#' @param A_est a p by p matrix of transition matrix \eqn{A} estimate.
#' @param sig2_eta scalar; estimate of propagation error variance \eqn{\sigma_\eta^2}.
#' @param sig2_epsilon scalar;  estimate of measurement error variance \eqn{\sigma_\epsilon^2}.
#' @param global_H0 a p by p matrix of global null hypothesis for transition matrix \eqn{A}.
#' If \code{global_H0=NULL}, global testing will not be conducted, and \code{global_idx} will not be used.
#' @param global_idx a p by p boolean matrix. The TRUE/nonzero entry indicates the entry of interest
#' in global hypothesis testing. If \code{global_idx=NULL}, all p*p entries are included in global testing.
#' @param simul_H0 a p by p matrix of simultaneous null hypothesis for transition matrix \eqn{A}.
#' If \code{simul_H0=NULL}, simultaneous testing will not be conducted, and (\code{simul_idx}, \code{FDR_levels}, \code{grid_num}) will not be used.
#' @param simul_idx a p by p boolean matrix. The TRUE/nonzero entry indicates the entry of interest
#' in simultaneous hypothesis testing. If \code{simul_idx=NULL}, all p*p entries are included in simultaneous testing.
#' @param FDR_levels a vector of FDR control levels
#' @param grid_num scalar; the number of grids for cutoff search in FDR control.
#'
#'
#' @return a list of testing results and gaussian test statistic matrices.
#' \tabular{ll}{
#' \code{pvalue}  \tab  scalar; p-value of global testing. Exist if \code{global_H0} is not NULL. \cr
#' \code{global_test_stat} \tab a p by p matrix of gaussian test statistic for global null hypothesis.
#' Exist if \code{global_H0} is not NULL. \cr
#' \code{simul_test_stat}  \tab  a p by p matrix of gaussian test statistic for simultaneous null hypothesis.
#' Exist if \code{simul_H0} is not NULL. \cr
#' \code{FDR_levels}  \tab a vector of FDR control levels. The same as input argument \code{FDR_levels}. \cr
#' \code{crt} \tab a vector of critical values for rejecting entries in simultaneous hypothesis
#' under corresponding FDR control levels.  \cr
#' \code{selected} \tab a three-way tensor. The first two modes are p by p, and the third mode is for FDR control levels.
#' Nonzero elements indicate rejected entries (the first two modes) in simultanous hypothesis at correspoding FDR control levels (the third mode).
#' The entries outside of \code{simul_idx} is set at zero.
#'
#' }
#'
#'
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#'
#' @examples
#'
#' p= 3; Ti=200  # dimension and time
#' A=diag(1,p) # transition matrix
#' sig_eta=sig_epsilon=0.2 # error std
#' Y=array(0,dim=c(p,Ti)) #observation t=1, ...., Ti
#' X=array(0,dim=c(p,Ti)) #latent t=1, ...., T
#' Ti_burnin=300 # time for burn-in to stationarity
#' for (t in 1:(Ti+Ti_burnin)) {
#'   if (t==1){
#'     x1=rnorm(p)
#'   } else if (t<=Ti_burnin) { # burn in
#'     x1=A%*%x1+rnorm(p,mean=0,sd=sig_eta)
#'   } else if (t==(Ti_burnin+1)){ # time series used for learning
#'     X[,t-Ti_burnin]=x1
#'     Y[,t-Ti_burnin]=X[,t-Ti_burnin]+rnorm(p,mean=0,sd=sig_epsilon)
#'   } else {
#'     X[,t- Ti_burnin]=A%*%X[,t-1- Ti_burnin]+rnorm(p,mean=0,sd=sig_eta)
#'     Y[,t- Ti_burnin]=X[,t- Ti_burnin]+rnorm(p,mean=0,sd=sig_epsilon)
#'   }
#' }
#'
#' # null hypotheses are true
#' hdVARtest(Y,A,sig_eta^2,sig_epsilon^2,global_H0=A,global_idx=NULL,
#'          simul_H0=A,simul_idx=NULL,FDR_levels=c(0.05,0.1))
#'
#'
#' # null hypotheses are false
#' hdVARtest(Y,A,sig_eta^2,sig_epsilon^2,global_H0=matrix(0,p,p),global_idx=NULL,
#'           simul_H0=matrix(0,p,p),simul_idx=NULL,FDR_levels=c(0.05,0.1))
#'
#' @export
#'
#'



hdVARtest=function(Y, A_est, sig2_eta, sig2_epsilon, global_H0=NULL,
                   global_idx=NULL, simul_H0=NULL, simul_idx=NULL, FDR_levels=0.05, grid_num=2000){

  if (is.null(global_H0) & is.null(simul_H0)){
    stop("Must specify one null hypothesis. Either global_H0 or simul_H0.")
  }

  Ti=dim(Y)[2] # the number of time points
  p=dim(Y)[1] # dimension

  if (is.null(simul_idx)){
    simul_idx=matrix(1,p,p)
  }

  if (is.null(global_idx)){
    global_idx=matrix(1,p,p)
  }

  res=Y[,2:Ti]-A_est%*%Y[,1:(Ti-1)]-apply(Y[,2:Ti]-A_est%*%Y[,1:(Ti-1)],1,mean) # centered residual
  res_cov=res[,2:(Ti-1)] %*% t(res[,1:(Ti-2)]) /(Ti-2) # lag-1 covariance of residual
  sig2_mat=entry_var(A_est=A_est,sig2_eta=sig2_eta,sig2_epsilon=sig2_epsilon)

  result=list()

  if (!is.null(global_H0)){ # gloabl testing
    # gaussian matrix for global testing
    bias_global= (sig2_eta+sig2_epsilon)*A_est-sig2_eta*global_H0 # bias correction term
    test_unstd_global=sqrt(Ti-2)*(res_cov+bias_global) # entry-wise asymptotic normal before standarization
    test_std_global=test_unstd_global/sqrt(sig2_mat) # entry-wise standard normal

    pvalue=global_pval(test_std_global[global_idx!=0])
    result=c(result,list(pvalue=pvalue,global_test_stat=test_std_global))
  }

  if (!is.null(simul_H0)){ # simul testing
    # gaussian matrix for simultaneous testing
    bias_simul= (sig2_eta+sig2_epsilon)*A_est-sig2_eta*simul_H0 # bias correction term
    test_unstd_simul=sqrt(Ti-2)*(res_cov+bias_simul) # entry-wise asymptotic normal before standarization
    test_std_simul=test_unstd_simul/sqrt(sig2_mat) # entry-wise standard normal


    simul_result=simul_FDR(test_std_simul,simul_idx, FDR_levels,grid_num)

    result=c(result,list(simul_test_stat=test_std_simul),simul_result)
  }

  return(result)
}
