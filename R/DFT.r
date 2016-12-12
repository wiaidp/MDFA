#' Generic function for computing the multivariate DFT: assigns a correct weight to frequency zero (1/sqrt(2)).
#'
#' @param insamp Length of in-sample portion of data (by default insamp=nrow(x))
#' @param x (Multivariate) Data: first column is the target series, the remaining columns are the explanatory series.
#' @param d For d=1 a pseudo-DFT is implemented (default is d=0)
#' @return weight_func (Multivariate) Discrete Fourier Transform
#'


spec_comp <- function(insamp, x, d) {
  # non-stationarity
  if(d == 1) {
    weight_func <- periodogram_bp(diff(x[1 : insamp, 1]), 1,
                                  insamp - 1)$fourtrans
    # explaining variables
    if(length(weight_func) > 1) {
      for(j in 2 : ncol(x)) {
        # since the data is integrated one uses the pseudo-periodogram:
        # diff(data) and d <- 1
        per <- periodogram_bp(diff(x[1 : insamp, j]), 1, insamp - 1)$fourtrans
        weight_func <- cbind(weight_func, per)
      }
    }
  } else {
    weight_func <- periodogram_bp(x[1 : insamp, 1], 0, insamp)$fourtrans
    # explaining variables
    if(length(weight_func) > 1) {
      for(j in 2 : ncol(x)) {
        per <- periodogram_bp(x[1 : insamp, j], 0, insamp)$fourtrans
        weight_func <- cbind(weight_func, per)
      }
    }
  }
  colnames(weight_func) <- colnames(x)

  return(list(weight_func = weight_func))
}



#' Periodogram and DFT function: can handle more general cases than function per()
#'
#' @param x Data
#' @param dd Integration order (default is d=0)
#' @param n.pg Length of time span considered for computing periodogram and DFT statistics
#' @return perall Periodogram
#' @return fourtrans DFT
#'

periodogram_bp <- function(x, dd, n.pg) {
  # preparations
  n.fit  <- length(x)
  xx     <- as.vector(x[((n.fit - n.pg + 1) : n.fit)])
  npg2   <- n.pg / 2
  perall <- fourtrans <- 0 * 0 : npg2

  # case without a seasonal component
  if (dd < 3) {
    for (j in 0 : npg2) {
      fourtrans[j + 1] <- xx %*% exp((1 : (2* npg2)) * 1.i * j * pi / npg2) #length(xx)
      fourtrans[j + 1] <- fourtrans[j + 1] / sqrt(pi * n.pg)
      perall[j + 1] <- abs(fourtrans[j + 1])^2
    }
  }

  # case with a seasonal component, special treatment for pi/6
  if (dd >= 3) {
    for (j in (1 : npg2)[(-npg2 / 6) * (1 : 6)]) {
      fourtrans[j + 1] <- xx %*% exp((1 : (2 * npg2)) * 1.i * j * pi / npg2)
      term2 <- abs(1 - exp(j * 1.i * pi / npg2))^2
      term3 <- abs(1 - exp(12 * j * 1.i * pi / npg2))^2
      perall[j + 1] <- abs(fourtrans[j + 1]) / (term2 * term3)
    }
    perall[(npg2 / 6) * (1 : 6) + 1] <- max(perall) * 100000
  }
# Weights wk: if length of data sample is even then DFT in frequency pi is scaled by 1/sqrt(2) (Periodogram in pi is weighted by 1/2)
  if (abs(as.integer(n.pg/2)-n.pg/2)<0.1)
  {
    fourtrans[npg2+1]<-fourtrans[npg2+1]/sqrt(2)
    perall[npg2+1]<-perall[npg2+1]/(2)
  }

  # output
  return(list(perall = perall, fourtrans = fourtrans))
}
