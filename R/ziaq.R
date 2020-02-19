#' Zero-inflation adjusted quantile regression for single cell RNA sequencing data
#'
#' This function fits the zero-inflation adjusted quantile regression model for the
#' differential expression analysis in single cell RNA sequencing data
#' @param Y_matrix a matrix for expression values with row representing indiviudal genes
#' and column representing cells
#' @param colDat a dataframe including the individual cell information
#' @param formula a formula with the predictors included in \code{colDat}. The default
#' is ~condition.
#' @param group the variable name in \code{colDat} for the factor used in group comparsion. The default
#' is condition.
#' @param probs the quantile levels for the quantile regressions.
#' The default is \code{c(0.25, 0.5, 0.75)}
#' @param log_i TRUE or FALSE indicate whether to apply log transformation.
#' The default is TRUE.
#' @param parallel TRUE or FALSE indicate whether to  apply parallel computing.
#' The default is TRUE.
#' @param no.core The number of cores used in parallel computing. The default
#' is all available cores \code{detectCores()}
#' @return \item{pvalue}{The p-values of all genes for testing the signficance
#' of the specified \code{group} variable.}
#' \item{res}{The full results from function \code{ziaq_fit} for all genes}
#' @keywords ziaq_fit
#' @export
#' @import quantreg
#' @import metap
#' @import parallel
#' @import stats
#'
#' @examples
#' #Use simuluated data
#'ymatrix = matrix(round(100* runif(100*150)), ncol = 100)
#'rownames(ymatrix) = paste0('gene', 1:150)
#'
#'colDat = data.frame(condition = rep(c(1, 0), e = 50))
#'
#'res = ziaq(ymatrix, colDat, formula = ~ condition,
#'           group = 'condition', probs = c(0.25, 0.5, 0.75),
#'           log_i = TRUE, parallel = FALSE, no.core = 1)
#'
#'print(res$pvalue)


ziaq <-function (Y_matrix, colDat, formula = ~ condition,
                     group = 'condition', probs = c(0.25, 0.5, 0.75),
                     log_i = TRUE, parallel = FALSE, no.core = detectCores() ) {

  #require(parallel)
  if(parallel ){
    cl <- makeCluster(getOption("cl.cores", no.core))
    res = parApply(cl = cl, Y_matrix, 1, ziaq_fit, colDat = colDat, formula = formula,
                   group = group, probs =probs,log_i = log_i )
  }else{
    res = apply(Y_matrix, 1, ziaq_fit, colDat = colDat, formula = formula,
               group = group, probs =probs,log_i = log_i )
    }

  pval = sapply(res, function(x) return(x$pvalue))
  return(list(pvalue = pval, full_results = res))
}
