#' Compute the spearman rank correlation between
#' the columns of UNI and the columns of DUO DECONV DATA.
#'
#' @param UNI
#' @param DUO
#'
#' @examples
#'
#' spearm <- compute_spearman_accuracy(DECONV_uni, DECONV_sub)
#'
compute_spearman_accuracy <- function(UNI, DUO) {
  x <- UNI@mat
  y <- DUO@mat[UNI@rid, UNI@cid]
  rho <- cor(t(x), t(y), method="spearman") # see ?cor
  diag(rho)
}

#
#' Compute AUC for hi and low extreme modulations from DE data
#'
#' @param UNI
#' @param DUO
#' @param th
#'
#' @examples
#'
#' auc <- compute_auc(DE_uni, DE_sub)
#'
compute_auc <- function(UNI, DUO, th=2) {
  require(pROC)
  x <- UNI@mat
  y <- DUO@mat[UNI@rid, UNI@cid]
  auc_hi <- roc(response=as.vector(x>=th), predictor=as.vector(y))$auc
  auc_lo <- roc(response=as.vector(x<=-th), predictor=as.vector(y))$auc
  return(c(auc_hi, auc_lo))
}

report <- function(DECONV.uni, DECONV.sub, DE.uni, DE.sub) {
  spearm <- compute_spearman_accuracy(DECONV.uni, DECONV.sub)
  auc <- compute_auc(DE.uni, DE.sub)
  message("Median Spearman correlation: ", round(median(spearm), 3))
  message("AUC (hi and lo): ", paste(round(auc, 3), collapse=" and "))
  message("Overall accuracy score 1e6 * COR * AUC = ", 1e6* mean(auc)*median(spearm))
}
