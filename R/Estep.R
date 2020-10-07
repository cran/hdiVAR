#' expectation step in sparse expectation-maximization algorithm
#'
#' @description Compute conditional expectation and covariance of \eqn{x_t} given \eqn{y_1,\ldots,y_T}
#' and current parameter estimates of \eqn{A, \sigma_\eta,\sigma_\epsilon}
#' via kalman filter and smoothing.
#'
#' @usage Estep(Y,A_init,sig_eta_init,sig_epsilon_init,X_init,P_init)
#'
#' @param Y observations of time series, a p by T matrix.
#' @param A_init current estimate of transition matrix \eqn{A}.
#' @param sig_eta_init  current estiamte of \eqn{\sigma_\eta}.
#' @param sig_epsilon_init current estiamte \eqn{\sigma_\epsilon}.
#' @param X_init current estimate of latent \eqn{x_1} at the first iteration, a p-dimensional vector.
#' @param P_init current covariance estimate of latent \eqn{x_1} at the first iteration, a p by p matrix.
#'
#' @return a list of conditional expectations and covariances for the sequential Maximization step.
#' \tabular{ll}{
#' \code{EXtT}  \tab  a p by T matrix of column \eqn{E[x_t | y_1,\ldots,y_T, \hat{A},\hat{\sigma}_\eta,\hat{\sigma}_\epsilon]}. \cr
#' \code{EXtt}  \tab  a p by p by T tensor of first-two-mode slice \eqn{E[x_t x_t^\top | y_1,\ldots,y_T, \hat{A},\hat{\sigma}_\eta,\hat{\sigma}_\epsilon]}. \cr
#' \code{EXtt1}  \tab  a p by p by T-1 matrix of first-two-mode slice \eqn{E[x_t x_{t+1}^\top | y_1,\ldots,y_T, \hat{A},\hat{\sigma}_\eta,\hat{\sigma}_\epsilon]}. \cr
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
#' Efit=Estep(Y,A,sig_eta,sig_epsilon,x1,diag(1,p))
#'
#'
#' @export
#'

Estep=function(Y,A_init,sig_eta_init,sig_epsilon_init,X_init,P_init){

  # kalman filter and smoothing
  Kalman_result=kalman(Y=Y,A=A_init,sig_eta=sig_eta_init,
                       sig_epsilon=sig_epsilon_init,
                       X_init=X_init,P_init=P_init)

  XtT=Kalman_result[["XtT"]]
  Xtt=Kalman_result[["Xtt"]]
  Xt1t=Kalman_result[["Xt1t"]]
  PtT=Kalman_result[["PtT"]]
  Lt=Kalman_result[["Lt"]]

  p=dim(XtT)[1]
  Ti=dim(XtT)[2]

  EXtt=array(0,dim=c(p,p,Ti)) # t= 1, ..., T
  EXtt1=array(0,dim=c(p,p,Ti-1)) # t=1, ..., T-1

  for ( t in 1:Ti) {
    EXtt[,,t]=PtT[,,t]+XtT[,t] %*% t(XtT[,t])
    if (t<Ti) {
      EXtt1[,,t]=Xtt[,t]%*% t(XtT[,t+1])+Lt[,,t] %*% (PtT[,,t+1]+ (XtT[,t+1]-Xt1t[,t])%*%t(XtT[,t+1]))
    }
  }

  return(list(EXtT=XtT,EXtt=EXtt,EXtt1=EXtt1,PtT=PtT,Xtt=Xtt))
}
