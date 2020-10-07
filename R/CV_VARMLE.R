#' cross-validation for transition matrix update in maximization step
#'
#' @description Tune the tolerance parameter of generalized Dantzig selector and hard thresholding
#' level via prediction error in test data.
#'
#' @param tol_seq  vector; grid of tolerance parameter in Dantzig selector for cross-validation.
#' @param ht_seq  vector; grid of hard-thresholding levels for transition matrix estimate.
#' To avoid hard thresholding, set \code{ht_seq=0}.
#' @param S0_train   a p by p matrix; average (over time points in training data) of conditional expectation of \eqn{x_t x_t^\top} on \eqn{y_1, \ldots, y_T} and parameter estimates, obtained from expectation step.
#' @param S1_train   a p by p matrix; average (over time points in training data) of conditional expectation of \eqn{x_t x_{t+1}^\top}on \eqn{y_1, \ldots, y_T} and parameter estimates, obtained from expectation step.
#' @param Y_test     a p by T_test matrix; observations of time series in test set.
#' @param is_echo  logical; if true, display the information of CV-optimal (tol, ht).
#'
#' @return a list of CV-optimal parameters and test prediction error.
#' \tabular{ll}{
#' \code{tol_min}  \tab  CV-optimal tolerance parameter in Dantzig selector. \cr
#' \code{ht_min}  \tab  CV-optimal hard thresholding level for the output of Dantzig selector. \cr
#' \code{test_loss}  \tab  a matrix of prediction error in test data; columns match \code{tol_seq}, and rows match \code{ht_seq}. \cr
#' }
#'
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#' @export
#'



CV_VARMLE=function(tol_seq,ht_seq,S0_train,S1_train,Y_test,is_echo=FALSE){

  test_loss=c() # test losses along ht_seq * tol_seq

  Ti_test=dim(Y_test)[2]
  # cross validation for tol and ht
  for (i_tol in tol_seq) {
    A_est=VARMLE(S0_train,S1_train,i_tol) # fit at current tol
    pred_loss=c()
    for (i_ht in ht_seq){
      A_est_thres=A_est
      A_est_thres[abs(A_est_thres)<i_ht]=0 # thresholding
      # prediction error in test data
      pred_loss=c(pred_loss,sum((Y_test[,2:Ti_test] - A_est_thres%*%Y_test[,1:(Ti_test-1)])^2))
    }
    test_loss=cbind(test_loss,pred_loss)
  }

  rownames(test_loss)=ht_seq
  colnames(test_loss)=tol_seq

  l_tol=which(test_loss == min(test_loss), arr.ind = TRUE)[1,2]
  l_ht=which(test_loss == min(test_loss), arr.ind = TRUE)[1,1]
  tol_min=tol_seq[l_tol]
  ht_min=ht_seq[l_ht]

  if (is_echo){
    print(paste("CV-tuned (lamda,ht) is in (", l_tol,",",l_ht,")/(",length(tol_seq),
              ",",length(ht_seq),") of the parameter grid.",sep=''))
  }

  return(list(tol_min=tol_min,ht_min=ht_min,test_loss=test_loss))
}
