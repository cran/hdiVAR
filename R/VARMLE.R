#' generalized Dantzig selector for transition matrix update in maximization step
#'
#' @description Sparse estimation of transtion matrix in vector autoregression given conditional autocovariance matrices.
#'
#'
#' @param S0   a p by p matrix; average (over time points) of conditional expectation of \eqn{x_t x_t^\top} on \eqn{y_1, \ldots, y_T} and parameter estimates, obtained from expectation step.
#' @param S1   a p by p matrix; average (over time points) of conditional expectation of \eqn{x_t x_{t+1}^\top}on \eqn{y_1, \ldots, y_T} and parameter estimates, obtained from expectation step.
#' @param tol  tolerance parameter in Dantzig selector.
#'
#' @return Sparse estimate of transition matrix by Dantzig selector.
#'
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#' @importFrom lpSolve lp
#'
#' @export
#'


VARMLE = function(S0,S1,tol){


  p=dim(S0)[1]

  # set up constraints for linear programming
  lp_constr=matrix(rep(0,8*p^2),nrow=4*p)
  lp_constr[1:p, 1:p]=-S0;lp_constr[1:p, (p+1):(2*p)]=S0
  lp_constr[(p+1):(2*p), 1:p]=S0;lp_constr[(p+1):(2*p), (p+1):(2*p)]=-S0
  lp_constr[(2*p+1):(4*p),]=diag(1,2*p)
  lp_obj=rep(1,2*p);lp_ineq=rep(">=",4*p)

  # linear programming/Dantzig selector
  A_est_stat_tmp=sapply(1:p, FUN=function(x){
    lp_rhs=matrix(c(-tol-S1[,x],-tol+S1[,x],rep(0,2*p)),4*p,1)
    lp_out=lp('min',lp_obj,lp_constr,lp_ineq,lp_rhs)$solution
    out=matrix(lp_out[1:p]-lp_out[(1+p):(p+p)],p,1)
    return(out)
  })

  A_est_stat=(A_est_stat_tmp)
  return(t(A_est_stat))
}
