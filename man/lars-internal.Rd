\name{lars-internal}
\title{Internal lars functions}
\alias{backsolvet}
\alias{cv.folds}
\alias{delcol} 
\alias{downdateR}
\alias{error.bars}
\alias{nnls.lars}
\alias{plotCVLars}
\alias{print.lars}
\alias{updateR}
\alias{.First.lib}
\description{Internal lars functions}
\usage{
backsolvet(r,x,k = ncol(r))
cv.folds(n, folds = 10)
delcol(r, z, k = p)
downdateR (R, k = p)
error.bars(x, upper, lower, width = 0.02, ...)
nnls.lars(active, Sign, R, beta, Gram, eps = 1e-10, trace = FALSE, 
    use.Gram = TRUE)
plotCVLars(cv.lars.object, se = TRUE)
updateR(xnew, R = NULL, xold, eps = .Machine$double.eps, Gram = FALSE)
.First.lib(lib, pkg)
}
\author{Trevor Hastie}
\details{
  These are not to be called by the user. \code{backsolvet} is included
  to make the R code compatible with the Splus code, since
  \code{backsolve} in R has a \code{transpose=TRUE} option already.
}
\keyword{internal}
