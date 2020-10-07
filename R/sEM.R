#' sparse expectation-maximization algorithm for high-dimensional vector autoregression with measurement error
#'
#' @description Alteranting between expectation step (by kalman filter and smoothing) and maximization step (by generalized Dantzig selector for transiton matrix)
#' to estimate transtion matrix and error variances.
#'
#'
#' @param Y observations of time series, a p by T matrix.
#
#' @param A_init  a p by p matrix as initial value of transition matrix \eqn{A} estimate.
#' @param sig2_eta_init scalar; initial value of propagation error variance \eqn{\sigma_\eta^2} estimate in latent signal process.
#' @param sig2_epsilon_init scalar; initial value of measurement error variance \eqn{\sigma_\epsilon^2} estimate in observation process.
#' @param Ti_train scalar; the number of time points in training data in cross-validation.
#' @param Ti_gap scalar; the number of time points between test data and train data in cross-validation.
#' @param tol_seq  vector; grid of tolerance parameter in Dantzig selector for cross-validation. If \code{is_cv=FALSE}, use the first element.
#' @param ht_seq  vector; grid of hard-thresholding levels for transition matrix estimate. If \code{is_cv=FALSE}, use the first element.
#' To avoid hard thresholding, set \code{ht_seq=0}.
#' @param is_cv logical; if true, run cross-validation to tune Dantzig selector tolerance parameter each sparse EM iteration.
#' @param count_vanish scalar; if the difference between updates of two consecutive
#' iterations is less that \code{thres} up to \code{count_vanish} times, the algorithm is terminated due to vanishing updates.
#' @param thres scalar; if the difference between updates of two consecutive iterations is less that \code{thres}, record one hit.
#' The algorithm is terminated due to vanishing updates if hit times accumulate up to \code{count_vanish}. If \code{thres=NULL}, the algorithm will not be terminated
#' due to vanishing updates, but too many iterations instead.
#' @param n_em scalar; the maximal allowed number of EM iterations, otherwise the algorithm is terminated due to too many iterations.
#' If \code{n_em=NULL}, the algorithm will not be terminated due to too many iterations, but vanishing updates instead.
#' @param is_echo logical; if true, display the information of CV-optimal (tol, ht) each iteration, and of algorithm termination.
#' @param is_sparse logical; if false, use standard EM algorithm.
#'
#'
#' @return a list of parameter estimates.
#' \tabular{ll}{
#' \code{A_est}  \tab  estimate of transition matrix \eqn{A}. \cr
#' \code{sig2_eta_hat}  \tab  estimate of propagation error variance \eqn{\sigma_\eta^2}. \cr
#' \code{sig2_epsilon_hat}  \tab   estimate of measurement error variance \eqn{\sigma_\epsilon^2}. \cr
#' \code{iter_err} \tab   the difference between updates of two consecutive iterations. \cr
#' \code{iter_err_ratio} \tab  the difference ratio (over the previous estimate) between updates of two consecutive iterations. \cr
#' }
#'
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#' @examples
#' p= 3; Ti=20  # dimension and time
#' A=diag(1,p) # transition matrix
#' sig_eta=sig_epsilon=0.2 # error std
#' Y=array(0,dim=c(p,Ti)) #observation t=1, ...., Ti
#' X=array(0,dim=c(p,Ti)) #latent t=1, ...., T
#' Ti_burnin=30 # time for burn-in to stationarity
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
#' sEM_fit=sEM(Y,diag(0.5,p),0.1,0.1,Ti*0.5,Ti*0.2,c(0.01,0.1))
#'
#'
#' @export
#'
#'



sEM = function(Y, A_init, sig2_eta_init, sig2_epsilon_init, Ti_train, Ti_gap, tol_seq, ht_seq=0, is_cv=TRUE, thres=1e-3, count_vanish=1, n_em=NULL, is_echo=FALSE, is_sparse=TRUE){

  if (is.null(thres)&is.null(n_em)){
    stop("Must choose one stopping criterion. Either vanishing updates (thres) or too many iterations (n_em).")
  }

  if (is.null(thres) != is.null(count_vanish)){
    stop("If choose vanishing update stopping criterion, please also specify tolerance (count_vanish).")
  }

  if (is_sparse){
    if(is.null(tol_seq)){
      stop("tol_seq must be non-empty.")
    }

    if(is.null(ht_seq)){
      stop("ht_seq must be non-empty.")
    }
  }


  # record update difference among iterations
  iter_err=c() # estiamte_new - estimate_old
  # stopping criteria
  count_stop=0 # stopping count
  count_em=1 # iteration count


  Ti=dim(Y)[2]
  p=dim(Y)[1]

  # assign initial values
  A_est=A_init
  sig_eta_hat=sqrt(sig2_eta_init)
  sig_epsilon_hat=sqrt(sig2_epsilon_init)

  X_hat=Y[,1] #E[X1|Y1] =\= Y1
  P_hat=diag(sig_epsilon_hat^2,nrow=p)

  if (is_cv & is_sparse){ # initals for CV
    X_hat_cv=Y[,Ti-Ti_train+1]
    P_hat_cv=diag(sig_epsilon_hat^2,nrow=p)
  }



  repeat { # EM iteration

    ####################################################
    ############### CV for transition matrix ###########
    ####################################################
    # choosing tuning parameters for A by EM intermediates from training data
    # time series is splitted into (test, gap, train)
    if (is_cv & is_sparse){
      Estep_train=Estep(Y=Y[,(Ti-Ti_train+1):Ti],A_init=A_est,
                        sig_eta_init=sig_eta_hat,
                        sig_epsilon_init=sig_epsilon_hat,
                        X_init=X_hat_cv,P_init=P_hat_cv)
      EXtT_train=Estep_train[["EXtT"]]
      Xtt_train=Estep_train[["Xtt"]]
      EXtt_train=Estep_train[["EXtt"]]
      EXtt1_train=Estep_train[["EXtt1"]]
      PtT_train=Estep_train[["PtT"]]

      # cross validation sparse MLE
      S0_train = apply(EXtt_train,c(1,2),mean)
      S1_train=apply(EXtt1_train,c(1,2),mean)
      sEM_train=CV_VARMLE(tol_seq,ht_seq,S0_train,S1_train,Y[,1:(Ti-Ti_train-Ti_gap)],is_echo = is_echo)
      tol_min=sEM_train$tol_min
      ht_min=sEM_train$ht_min
    } else { # if not CV, use the first element.
      tol_min=tol_seq[1]
      ht_min=ht_seq[1]
    }

    ###################################################
    ################ Expecation step ##################
    ###################################################
    # use all time points
    Estep_result=Estep(Y=Y,A_init=A_est,sig_eta_init=sig_eta_hat,
                       sig_epsilon_init=sig_epsilon_hat,
                       X_init=X_hat,P_init=P_hat)
    EXtT=Estep_result[["EXtT"]]
    Xtt=Estep_result[["Xtt"]]
    EXtt=Estep_result[["EXtt"]]
    EXtt1=Estep_result[["EXtt1"]]
    PtT=Estep_result[["PtT"]]


    ##########################################################
    #### estimate A by sparse MLE via linear programming  ####
    ##########################################################
    A_old=A_est # old estimate
    S0=apply(EXtt,c(1,2),mean)
    if (is_sparse){
      S1=apply(EXtt1,c(1,2),mean)
      A_est=VARMLE(S0,S1,tol_min)
      A_est_unthres=A_est
      A_est[abs(A_est)<ht_min]=0
    }

    ####################################################
    ############ MLE for error variance ################
    ####################################################
    sig_eta_old=sig_eta_hat # old estimate
    sig_epsilon_old=sig_epsilon_hat

    Mstep_result=Mstep(Y=Y,A=A_est,EXtT=EXtT,EXtt=EXtt,EXtt1=EXtt1,is_MLE=!is_sparse)
    # update error variance
    sig_eta_hat=Mstep_result[["sig_eta"]]  # new estimate
    sig_epsilon_hat=Mstep_result[["sig_epsilon"]]

    if (!is_sparse){ # naive MLE of A
      A_est=Mstep_result[["A"]]
    }

    #####################################################
    #### update initials given current EM estimates #####
    #####################################################
    # estimate of covariance of intial x
    cov_X_est=S0
    proj_var=cov_X_est %*% solve(cov_X_est+diag(sig_epsilon_hat^2,p))
    P_hat=cov_X_est- proj_var %*% t(cov_X_est)
    X_hat= proj_var%*%Y[,1]

    if (is_sparse){
      X_hat_cv=proj_var %*%Y[,Ti-Ti_train+1]
      P_hat_cv=P_hat
    }


    # record update difference
    iter_err_ratio=cbind(iter_err,c(sum((A_est-A_old)^2)/(sum(A_old^2)+1e-10),
                                    abs(sig_eta_hat^2-sig_eta_old^2)/sig_eta_old^2,
                                    abs(sig_epsilon_hat^2-sig_epsilon_old^2)/sig_epsilon_old^2))
    iter_err=cbind(iter_err,c(sum((A_est-A_old)^2),
                              abs(sig_eta_hat^2-sig_eta_old^2),
                              abs(sig_epsilon_hat^2-sig_epsilon_old^2)))

    rownames(iter_err)=rownames(iter_err_ratio)=c("A","sig2_eta","sig2_epsilon")
    colnames(iter_err)=colnames(iter_err_ratio)=paste("iter",1:ncol(iter_err),sep="")

    # check stopping criteria
    if (!is.null(thres)){
      if(min(iter_err[,count_em])<thres) {# negligible update
        count_stop=count_stop+1
        if (count_stop==count_vanish){
          if (is_echo){
            cat(paste("sparse EM is terminated due to vanishing updates",sep=''),fill=TRUE)
          }
          break
        }
      } else {
        count_stop=0
      }
    }


    count_em=count_em+1
    if (!is.null(n_em)){
      if (count_em==n_em) { # too many iterations
        if(is_echo){
          cat(paste("sparse EM is terminated due to too many iterations.",sep=''),fill=TRUE)
        }
        break
      }
    }
  }


  return(list(A_est=A_est,sig2_epsilon_hat=sig_epsilon_hat^2,
              sig2_eta_hat=sig_eta_hat^2,iter_err=iter_err,
              iter_err_ratio=iter_err_ratio))

}
