#' Zero-inflation adjusted quantile regression
#'
#' This function fits the zero-inflation adjusted quantile regression model
#' for one individual gene.
#' @param y a numeric vector including gene expression value for one gen
#' @param colDat a dataframe including the individual cell information
#' @param formula a formula with the predictors included in \code{colDat}. The default
#' is ~condition.
#' @param group the variable name in \code{colDat} for the factor used in group comparsion. he default
#' is condition.
#' @param probs the quantile levels for the quantile regressions.
#' The default is \code{c(0.25, 0.5, 0.75)}
#' @param log_i TRUE or FALSE indicate whether to apply log transformation.
#' The default is TRUE.
#' @return \item{pvalue}{The p-value for testing the signficance of the specified \code{group} variable.}
#' \item{logistic_model}{The \code{glm} object of the logistics regression for fitting binary
#' outcome (zero or not zero) vs the predictors specified in \code{formula}.}
#' \item{zero_proportion}{The proportion of zero measurements}
#' \item{quantile_levels}{The input quantile levels and zero-proportion adjusted quantile levels}
#' \item{qr}{The objects of the fitted quantile regressions}
#' @keywords ziaq_fit
#' @export
#' @import quantreg
#' @import metap
#' @import parallel
#' @import stats
#'
#' @examples
#' # simulated data
#'y = round(100* runif(100))
#'colDat = data.frame(condition = rep(c(1, 0), e = 50))
#'res = ziaq_fit(y, colDat = colDat,  formula = ~ condition,
#'               group = 'condition', probs = c(0.25, 0.5, 0.75),
#'               log_i = T )
#'print(res)


ziaq_fit<-function (y, colDat, formula = ~ condition,
                    group = 'condition', probs = c(0.25, 0.5, 0.75),
                    log_i = T) {

  set.seed(124)
  #require(quantreg)
  #require(metap)
  data = data.frame(y = y, colDat)
  out = list()

  # Step 1: logistics regression on dropout events
  y2 = I(y > 0 )
  data[, 'y_ind'] = y2

  log_p = 1
  if(length(table(y2))>1){
    formula_logit = paste0('y_ind', deparse(formula))
    log_m = glm (formula_logit, data =data, family = 'binomial')
    coef_h = summary(log_m)$coef
    log_p = coef_h[grep(group,rownames(coef_h)),  "Pr(>|z|)"  ]
    out$logistic_model = log_m
  }
  out$zero_proportion = 1 - sum(y2, na.rm = T)/length(y2)

  # Step 2: impuate y by adding perturbation around zero
  yi  = y
  yi[!y2] = NA
  m = min(yi, na.rm = T)
  if(!is.finite(m)) m = 10
  yi[!y2] = runif(sum(!y2), min = 0, max = m/5)

  if (log_i){
    yi = log(yi+1)
  }
  data[, 'y_impute'] = yi

  # Step 4: zero inflation adjusted quantile regression at various level

  probs_nonzero = probs
  probs = 1-(1-probs)*sum(y2)/length(y2)
  out$quantile_levels = data.frame(Input = probs_nonzero,
                                   Adjusted = probs)
  out$qr = list()

  formula_q = paste0('y_impute', deparse(formula))
  prob_p = NULL
  for (p in probs){
    rq_m = rq(formula_q, tau=p, data= data, method="br")
    coef_h = summary(rq_m, se = 'boot')$coef
    rq_p = coef_h[grep(group,rownames(coef_h)),  "Pr(>|t|)"  ]
    rq_p = max(rq_p, 1e-16)
    prob_p = c(prob_p, rq_p)
    out$qr[[as.character(p)]] = rq_m
  }

  rq_p_all = sumlog(c(prob_p, log_p))$p
  out$pvalue = rq_p_all

  return(out)
}


