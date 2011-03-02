\name{plot.lars}
\alias{plot.lars}
\title{Plot method for lars objects}
\description{
Produce a plot of a lars fit. The default is a complete coefficient
path.
}
\usage{
\method{plot}{lars}(x, xvar= c("norm", "df", "arc.length", "step"), breaks = TRUE, plottype = c("coefficients", 
    "Cp"), omit.zeros = TRUE, eps = 1e-10, ...)
}
\arguments{
  \item{x}{lars object}
  \item{xvar}{The type of x variable against which to
plot. \code{xvar=norm} (default) plots against the L1 norm of the coefficient
vector, as a fraction of the maximal L1 norm. \code{xvar=step} plots against the step number (which is
essentially degrees of freedom for LAR; not for LASSO or Forward
Stagewise). \code{xvar=arc.length} plots against the arc.length of the
fitted vector; this is useful for a LAR object, because the L1 norm of
its coefficient vector need not be monotone in the steps.
\code{xvar=df} plots against the estimated df, which is the size of the
active set at each step.}
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
  \details{The default plot uses the fraction of L1 norm as the
    xvar. For forward stagewise and LAR, coefficients can pass through
    zero during a step, which causes a change of slope of L1 norm vs
    arc-length. Since the coefficients are piecewise linear in
    arc-length between each step, this causes a change in slope of the
    coefficients.}
\value{NULL}
\references{
Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}; see also
\url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}
Yann-Ael Le Borgne (private communication) pointed out the problems in
plotting forward stagewise and LAR coefficients against L1 norm, and the
solution we have implemented.
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
