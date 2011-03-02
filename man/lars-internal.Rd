\name{lars-internal}
\title{Internal lars functions}
\alias{betabreaker}
\alias{backsolvet}
\alias{cv.folds}
\alias{delcol} 
\alias{downdateR}
\alias{error.bars}
\alias{nnls.lars}
\alias{plotCVLars}
\alias{print.lars}
\alias{updateR}
\description{Internal lars functions}
\usage{
betabreaker(object)
backsolvet(r,x,k = ncol(r))
cv.folds(n, folds = 10)
delcol(r, z, k = p)
downdateR (R, k = p)
error.bars(x, upper, lower, width = 0.02, ...)
nnls.lars(active, Sign, R, beta, Gram, eps = 1e-10, trace = FALSE, 
    use.Gram = TRUE)
plotCVLars(cv.lars.object, se = TRUE)
updateR(xnew, R = NULL, xold, eps = .Machine$double.eps, Gram = FALSE)
}
\author{Trevor Hastie}
\details{
  These are not to be called by the user. \code{betabreaker} figures out
  if coefficients (other than lasso) pass through zero, since the L1
  norm is discontinuous there, and this has an impact on
  predict/plot. Suggested by Yann-Ael Le Borgne. \code{backsolvet} is included
  to make the R code compatible with the Splus code, since
  \code{backsolve} in R has a \code{transpose=TRUE} option already.
}
\keyword{internal}
