\name{plot.lars}
\alias{plot.lars}
\title{Plot method for lars objects}
\description{
Produce a plot of a lars fit. The default is a complete coefficient
path.
}
\usage{
plot.lars(x, xvar= c("norm", "df", "arc.length"), breaks = TRUE, plottype = c("coefficients", 
    "Cp"), omit.zeros = TRUE, eps = 1e-10, ...)
}
\arguments{
  \item{x}{lars object}
  \item{xvar}{The type of x variable against which to
plot. \code{xvar=norm} plots against the L1 norm of the coefficient
vector (default). \code{xvar=df} plots against the step number (which is
essentially degrees of freedom for LAR; not for LASSO or Forward
Stagewise). \code{xvar=arc.length} plots against the arc.length of the
fitted vector; this is useful for a LAR object, because the L1 norm of
its coefficient vector need not be monotone in the steps.}
  \item{breaks}{If \code{TRUE}, then vertical lines are drawn at each
      break point in the piecewise linear coefficient paths}
  \item{plottype}{Either \code{coefficients} (default) or \code{Cp}. The
      coefficient plot shows the path of each coefficient as a function
      of the norm fraction or Df. The Cp plot shows the Cp curve.}
  \item{omit.zeros}{When the number of variables  is much greater than
    the number of observations, many coefficients will never be nonzero;
      this logical (default \code{TRUE}) avoids plotting these zero coefficents}
  \item{eps}{Definition of zero above, default is \code{1e-10}}
  \item{\dots}{Additonal arguments for generic plot. Can be used to set xlims,
      change colors, line widths, etc}
}
\value{NULL}
\references{
Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}; see also
\url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.ps}
}
\author{Trevor Hastie}
\examples{
data(diabetes)
attach(diabetes)
object <- lars(x,y)
plot(object)
detach(diabetes)
}
\keyword{hplot}
\keyword{methods}
