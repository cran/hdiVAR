#' entry variance estimation in residual autocovariance matrix
#'
#' @description Generate a matrix of entry variance estimate to rescale test statistic into marginally standard normal.
#'
#'
#' @noRd
#'
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'

entry_var=function(A_est,sig2_eta,sig2_epsilon){

  V1=dim(A_est)[1];V2=dim(A_est)[2]
  sig2_mat=matrix(0,V1,V2)
  for (i in 1:V1){
    for (j in 1:V2){
      sig2_mat[i,j] = (sig2_eta+sig2_epsilon)^2 + sig2_epsilon^2 * A_est[i,j]^2+
        2*sig2_epsilon^2*A_est[i,i]*A_est[j,j]+
        sig2_epsilon^2*sum(A_est[i,]^2)*sum(A_est[j,]^2)+
        (sig2_epsilon^2+sig2_epsilon*sig2_eta)*(sum(A_est[i,]^2)+sum(A_est[j,]^2))
    }
  }

  return(sig2_mat)
}
