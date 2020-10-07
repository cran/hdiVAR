#' Simultaneous testing under FDR control
#'
#'
#' @param test_stat a p by p matrix of test statistics.
#' @param idx_set a p by p boolen matrix. The TRUE/nonzero entries indicate interested entries in testing.
#' @param FDR_levels a vector of FDR control levels.
#' @param grid_num scalar; the number of grids for cutoff search in FDR control.
#'
#' @noRd
#' @author Xiang Lyu, Jian Kang, Lexin Li
#'
#' @importFrom abind abind
#' @importFrom stats pnorm
#'
#'


simul_FDR=function(test_stat,idx_set, FDR_levels,grid_num){


  S=sum(idx_set!=0)
  stat_vec=test_stat[idx_set!=0]

  crt_result=NULL
  selected=NULL
  for ( i_level in FDR_levels){
    # choose thresholding value
    for ( crt_FDR in seq(0,sqrt(2*log(S)),length.out = grid_num)){
      if ((2-2*pnorm(crt_FDR))*S/max(1,sum(abs(stat_vec)>crt_FDR))<=i_level){
        break
      }
    }

    tmp= (abs(test_stat)>crt_FDR) # 0-1 matrix of model selection
    tmp[idx_set==0]=0
    selected=abind(selected,tmp,along=3)
    crt_result=c(crt_result,crt_FDR)
  }

  if (length(FDR_levels)==1){selected=selected[,,1]}

  return(list(crt=crt_result,FDR_levels=FDR_levels,selected=selected))
}
