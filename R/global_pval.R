#' P-value of global testing
#'
#' @param test_stat a vector of test statistics on interested entries.
#'
#'
#' @noRd
#'
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#'

global_pval=function(test_stat){
  S=length(test_stat)
  gumbel_test=max(test_stat^2)
  crt=gumbel_test-2*log(S)+ log(log(S))
  pvalue=1-exp(-exp(-crt/2)/sqrt(pi))  # pvalue
  return(pvalue)
}
