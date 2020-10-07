#' maximization step of sparse expectation-maximization algorithm for updating error standard deviations
#'
#' @description Update \eqn{\sigma_\eta,\sigma_\epsilon} based on estimate of A and
#' conditional expecation and covariance from expectation step.
#'
#' @usage Mstep(Y,A,EXtT,EXtt,EXtt1,is_MLE=FALSE)
#'
#' @param Y observations of time series, a p by T matrix.
#' @param A current estimate of transition matrix \eqn{A}. If \code{is_MLE=TRUE}, use naive MLE of transition matrix, by conditional expecation and covariance from expecation step, to update error standard deviations.
#' @param EXtT a p by T matrix of column \eqn{E[x_t | y_1,\ldots,y_T, \hat{A},\hat{\sigma}_\eta,\hat{\sigma}_\epsilon]} from expectation step.
#' @param EXtt   a p by p by T tensor of first-two-mode slice \eqn{E[x_t x_t^\top | y_1,\ldots,y_T, \hat{A},\hat{\sigma}_\eta,\hat{\sigma}_\epsilon]} from expectation step.
#' @param EXtt1   a p by p by T-1 matrix of first-two-mode slice \eqn{E[x_t x_{t+1}^\top | y_1,\ldots,y_T, \hat{A},\hat{\sigma}_\eta,\hat{\sigma}_\epsilon]} from expectation step.
#' @param is_MLE logic; if true, use naive MLE of transition matrix, by conditional expecation and covariance from expecation step, to update error variances. Otherwise, use input argument \code{A}.
#'
#' @return a list of estimates of error standard deviations.
#' \tabular{ll}{
#' \code{sig_eta}  \tab  estimate of \eqn{\sigma_\eta}. \cr
#' \code{sig_epsilon}  \tab   estimate of \eqn{\sigma_\epsilon}. \cr
#' \code{A}  \tab  naive MLE of transition matrix \eqn{A} by conditional expecation and covariance from expecation step. Exist if \code{is_MLE=TRUE}. \cr
#' }
#'
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#' @examples
#' p= 2; Ti=10  # dimension and time
#' A=diag(1,p) # transition matrix
#' sig_eta=sig_epsilon=0.2 # error std
#' Y=array(0,dim=c(p,Ti)) #observation t=1, ...., Ti
#' X=array(0,dim=c(p,Ti)) #latent t=1, ...., T
#' Ti_burnin=100 # time for burn-in to stationarity
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
#' # expectation step
#' Efit=Estep(Y,A,sig_eta,sig_epsilon,x1,diag(1,p))
#' EXtT=Efit[["EXtT"]]
#' EXtt=Efit[["EXtt"]]
#' EXtt1=Efit[["EXtt1"]]
#' # maximization step for error standard deviations
#' Mfit=Mstep(Y,A,EXtT,EXtt,EXtt1)
#'
#'
#'
#' @export
#'


Mstep=function(Y,A,EXtT,EXtt,EXtt1,is_MLE=FALSE){

  p=dim(EXtt)[1]
  Ti=dim(EXtt)[3]

  SXtt1=apply(EXtt1,c(1,2),sum)
  SXt1t1=apply(EXtt[,,2:Ti],c(1,2),sum)  # sum over t=2, ...,T for estimating eta
  SXt_1t_1=apply(EXtt[,,1:(Ti-1)],c(1,2),sum)  # sum over t=1, ...,T-1 for estimating eta

  if (is_MLE){
    SXtt=apply(EXtt[,,1:(Ti-1)],c(1,2),sum)  # sum over t=1, ...,T-1 for estimating MLE of A
    A_est=t(SXtt1)%*% solve(SXtt) # naive MLE
    sig_eta=sqrt((sum(diag(SXt1t1)) -sum(diag(A_est %*% SXtt1))) /p/(Ti-1)) # use MLE A
  } else {
    sig_eta=sqrt((sum(diag(SXt1t1)) -sum(diag(A %*% SXtt1))) /p/(Ti-1)) # use sparse A
  }


  SYtt=norm(Y,'f')^2
  SYXt= sum(diag(Y%*%t(EXtT)))
  SXtt=apply(EXtt,c(1,2),sum) # sum over t=1, ..., T for estimating epsilon
  sig_epsilon=sqrt( (SYtt-2*SYXt + sum(diag(SXtt))) /p/Ti)

  if (!is_MLE){
    return(list(sig_eta=sig_eta,sig_epsilon=sig_epsilon))
  } else {
    return(list(sig_eta=sig_eta,sig_epsilon=sig_epsilon,A=A_est))
  }

}
