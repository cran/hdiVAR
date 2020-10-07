#' kalman filtering and smoothing for vector autoregression with measurement error
#'
#'
#' @usage kalman(Y,A,sig_eta,sig_epsilon,X_init=NULL,P_init=NULL)
#'
#' @param Y observations of time series, a p by T matrix.
#' @param A current estimate of transition matrix.
#' @param sig_eta current estiamte of \eqn{\sigma_\eta}.
#' @param sig_epsilon current estiamte \eqn{\sigma_\epsilon}.
#' @param X_init inital estimate of latent \eqn{x_1} at the first iteration, a p-dimensional vector.
#' @param P_init inital covariance estimate of latent \eqn{x_1} at the first iteration, a p by p matrix.
#'
#' @return a list of conditional expectations and covariances of \eqn{x_t}'s.
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#'
#' @export
#'


kalman=function(Y,A,sig_eta,sig_epsilon,X_init=NULL,P_init=NULL){

  p=dim(Y)[1];Ti=dim(Y)[2]

  ## Kalman filter
  Kt=array(0,dim=c(p,p,Ti-1)) # t=2, ..., T
  Xt1t=array(0,dim=c(p,Ti-1)) # t=1, ..., T-1
  Xtt=array(0,dim=c(p,Ti)) # t=1, ..., T
  Pt1t=array(0,dim=c(p,p,Ti-1)) # t=1, ..., T-1
  Ptt=array(0,dim=c(p,p,Ti)) # t=1, ..., T

  Xtt[,1]=X_init; Ptt[,,1]=P_init # init t=1

  for (t in 1:(Ti-1)) {
    Xt1t[,t]=A %*% Xtt[,t]
    Pt1t[,,t]=A%*% Ptt[,,t] %*% t(A)+diag(sig_eta^2,p)
    Kt[,,t]=Pt1t[,,t] %*% solve(Pt1t[,,t]+diag(sig_epsilon^2,p))
    Xtt[,t+1]=Xt1t[,t]+Kt[,,t]%*%(Y[,t+1]-Xt1t[,t])
    Ptt[,,t+1]=Pt1t[,,t]-Kt[,,t]%*% Pt1t[,,t]
  }

  ## Kalman smoother
  XtT=array(0,dim=c(p,Ti)) # t=1, ..., T
  Lt=array(0,dim=c(p,p,Ti-1)) # t=1, ..., T-1
  PtT=array(0,dim=c(p,p,Ti)) # t=1, ..., T

  XtT[,Ti]=Xtt[,Ti] # assign XTT from filtering
  PtT[,,Ti]=Ptt[,,Ti] # assign PTT from filtering
  for (t in (Ti-1):1){
    Lt[,,t]=Ptt[,,t] %*% t(A) %*% solve(Pt1t[,,t])
    XtT[,t]=Xtt[,t] + Lt[,,t] %*% (XtT[,t+1]-Xt1t[,t])
    PtT[,,t]=Ptt[,,t]+Lt[,,t] %*% (PtT[,,t+1]-Pt1t[,,t])%*% t(Lt[,,t])
  }

  return(list(XtT=XtT,Xtt=Xtt,Xt1t=Xt1t,PtT=PtT,Lt=Lt))
}
