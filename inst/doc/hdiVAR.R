## -----------------------------------------------------------------------------
library(hdiVAR)

set.seed(123)
p=3; Ti=400  # dimension and time
A=diag(1,p) # transition matrix
sig_eta=sig_epsilon=0.2 # error std
Y=array(0,dim=c(p,Ti)) #observation t=1, ...., Ti
X=array(0,dim=c(p,Ti)) #latent t=1, ...., T
Ti_burnin=400 # time for burn-in to stationarity
for (t in 1:(Ti+Ti_burnin)) {
 if (t==1){
   x1=rnorm(p)
 } else if (t<=Ti_burnin) { # burn in
   x1=A%*%x1+rnorm(p,mean=0,sd=sig_eta)
 } else if (t==(Ti_burnin+1)){ # time series used for learning
   X[,t-Ti_burnin]=x1
   Y[,t-Ti_burnin]=X[,t-Ti_burnin]+rnorm(p,mean=0,sd=sig_epsilon)
 } else {
   X[,t- Ti_burnin]=A%*%X[,t-1- Ti_burnin]+rnorm(p,mean=0,sd=sig_eta)
   Y[,t- Ti_burnin]=X[,t- Ti_burnin]+rnorm(p,mean=0,sd=sig_epsilon)
 }
}

## -----------------------------------------------------------------------------
# cross-validation grid of tolerance parameter \tau_k in Dantzig selector.
tol_seq=c(0.0001,0.0003,0.0005) 

# cross-validation grid of hard thresholding levels in transition matrix estimate. 
# Set as zero to avoid thresholding. The output is \hat{A}_k.
ht_seq=0 


A_init=diag(0.1,p) # initial estimate of A 
# initial estimates of error variances 
sig2_eta_init=sig2_epsilon_init=0.1 

# the first half time points are training data 
Ti_train=Ti*0.5

# The latter 3/10 time points are test data (drop out train (1/2) and gap (1/5) sets).
Ti_gap=Ti*0.2 

# sparse EM algorithm 
sEM_fit=sEM(Y,A_init,sig2_eta_init,sig2_epsilon_init,Ti_train,Ti_gap,tol_seq,ht_seq,is_echo = TRUE)

# estimate of A 
sEM_fit$A_est 

# estimate of error variances 
c(sEM_fit$sig2_epsilon_hat,sEM_fit$sig2_eta_hat) 

## -----------------------------------------------------------------------------

# use sparse EM estimates to construct test. Alternative consistent estimators can also be adopted if any. 
# test the entire matrix.


# FDR control levels for simultaneous testing 
FDR_levels=c(0.05,0.1)


# if null hypotheses are true (null hypothesis is true A): 
# p-value should > 0.05, and simultaneous testing selects no entries. 
true_null=hdVARtest(Y,sEM_fit$A_est,sEM_fit$sig2_eta_hat,sEM_fit$sig2_epsilon_hat,
                    global_H0=A,global_idx=NULL,simul_H0=A,
                    simul_idx=NULL,FDR_levels=FDR_levels)

# global pvalue: 
true_null$pvalue
# selection at FDR=0.05 control level 
true_null$selected[,,FDR_levels==0.05]



# if null hypotheses are false (null hypothesis is zero matrix): 
# p-value should < 0.05, and simultaneous testing selects diagnoal entries. 
false_null=hdVARtest(Y,sEM_fit$A_est,sEM_fit$sig2_eta_hat,sEM_fit$sig2_epsilon_hat,
                     global_H0=matrix(0,p,p),global_idx=NULL,simul_H0=matrix(0,p,p),
                     simul_idx=NULL,FDR_levels=c(0.05,0.1))

# global pvalue: 
false_null$pvalue
# selection at FDR=0.05 control level 
false_null$selected[,,FDR_levels==0.05]

