\name{cv.lars}
\alias{cv.lars}
\title{Computes K-fold cross-validated error curve for lars}
\description{
Computes the K-fold cross-validated mean squared prediction error for
lars, lasso, or forward stagewise.}
\usage{
cv.lars(x, y, K = 10, index, 
    trace = FALSE, plot.it = TRUE, se = TRUE,type = c("lasso", "lar", "forward.stagewise", "stepwise"),
         mode=c("fraction", "step"), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Input to lars}
  \item{y}{Input to lars}
  \item{K}{Number of folds}
  \item{index}{Abscissa values at which CV curve should be computed.
    If \code{mode="fraction"} this is the  fraction of the saturated
    |beta|. The default value in this case  is \code{index=seq(from = 0, to =
      1, length =100)}.
   If \code{mode="step"}, this is the number of steps in lars
   procedure. The default is complex in this case, and depends on
   whether \code{N>p} or not. In principal it is \code{index=1:p}.
 Users can supply their own values of index (with care).}
  \item{trace}{Show computations?}
  \item{plot.it}{Plot it?}
  \item{se}{Include standard error bands?}
  \item{type}{type of \code{lars} fit, with default \code{"lasso"}}
  \item{mode}{This refers to the index that is used for
    cross-validation. The default is \code{"fraction"} for
    \code{type="lasso"} or \code{type="forward.stagewise"}. For
    \code{type="lar"} or \code{type="stepwise"} the default is \code{"step"}}
  \item{\dots}{Additional arguments to \code{lars}}
}
\value{
Invisibly returns a list with components (which can be plotted using \code{plotCVlars})
  \item{index}{As above}
  \item{cv}{The CV curve at each value of index}
  \item{cv.error}{The standard error of the CV curve}
  \item{mode}{As above}
}
\references{
Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}; see also
\url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}.
}
\author{Trevor Hastie}
\examples{
data(diabetes)
attach(diabetes)
cv.lars(x2,y,trace=TRUE,max.steps=80)
detach(diabetes)
}
\keyword{regression}

