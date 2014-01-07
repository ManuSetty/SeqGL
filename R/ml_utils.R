
#' @import Matrix
#' @import Biostrings
#' @export


# Determine area under ROC

get.auc <- function (labels, predictions, pdf=NULL) {

  if (!all (labels %in% c(1, -1)))
    stop ('Labels have to +1/-1')
  
  ## Ranks for predictions
  sign.preds <- sign (predictions)
  ranks <- rank (predictions, ties.method='random')
  preds <- ranks[c(which (labels==1), which (labels==-1))]
  labels <- sort (labels, decreasing=TRUE)
  labels[labels == -1] <- 0

  ## Sorted order
  inds <- sort (preds, decreasing=TRUE, index.return=TRUE)$ix
  preds <- preds[inds]
  labels <- labels[inds]
  sign.preds <- sign.preds[inds]

  ## Determine true and false positive rates
  truepos <- cumsum (labels)
  falsepos <- 1:length (predictions) - truepos
  flags <- c(diff(predictions), 1) != 0
  truepos <- truepos[flags]/length (which (sign.preds == 1 & labels == 1 | sign.preds != 1 & labels == 1))
  falsepos <- falsepos[flags]/length (which (sign.preds == 1 & labels == 0 | sign.preds != 1 & labels == 0))

  ## Calculate auc
  end <- length (which (flags))
  auc <- sum((falsepos[2:end] - falsepos[1:(end-1)]) *
             (truepos[2:end] + truepos[1:(end-1)])/2)

  # Plot auc
  if (!is.null (pdf)) {
    pdf (pdf)
    plot (falsepos, truepos, type='l', xlim=c(0, 1), ylim=c(0, 1),
      xlab='False positive rate', ylab='True positive rate')
    lines (c(-10, 10), c(-10, 10), col='grey', lty=2)

    # Add auc
    legend ('bottomright', legend=sprintf ("Test auROC: %.4f", auc), bty='n', lty=1, col='white')
    dev.off ()

  }

  res <- list()
  res$auc <- auc; res$truepos <- truepos; res$falsepos <- falsepos
  invisible (res)
}


