\name{cv.lars}
\alias{cv.lars}
\title{Computes K-fold cross-validated error curve for lars}
\description{
Computes the K-fold cross-validated mean squared prediction error for
lars, lasso, or forward stagewise.}
\usage{
cv.lars(x, y, K = 10, fraction = seq(from = 0, to = 1, length = 100), 
    trace = FALSE, plot.it = TRUE, se = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Input to lars}
  \item{y}{Input to lars}
  \item{K}{Number of folds}
  \item{fraction}{Abscissa values at which CV curve should be computed,
      as a fraction of the saturated |beta|. Default is \code{seq(from = 0, to = 1, length = 100)}}
  \item{trace}{Show computations?}
  \item{plot.it}{Plot it?}
  \item{se}{Include standard error bands?}
  \item{\dots}{Additional arguments to \code{lars}}
}
\value{
Invisibly returns a list with components (which can be plotted using \code{plotCVlars})
  \item{fraction}{As above}
  \item{cv}{The CV curve at each value of fraction}
  \item{cv.error}{The standard error of the CV curve}
}
\references{
Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}; see also
\url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.ps}.
}
\author{Trevor Hastie}
\examples{
data(diabetes)
attach(diabetes)
cv.lars(x2,y,trace=TRUE,max.steps=80)
detach(diabetes)
}
\keyword{regression}

