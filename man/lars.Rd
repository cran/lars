\name{lars}
\alias{lars}
\title{
Fits Least Angle Regression, Lasso and Infinitesimal Forward Stagewise regression models
}
\description{
These are all variants of Lasso, and provide the entire sequence of
coefficients and fits, starting from zero, to the least squares fit.
}
\usage{
lars(x, y, type = c("lasso", "lar", "forward.stagewise", "stepwise"), 
    trace = FALSE, normalize = TRUE, intercept = TRUE, Gram, eps = .Machine$double.eps, max.steps, use.Gram = TRUE)
}
\arguments{
\item{x}{
matrix of predictors
}
\item{y}{
response 
}
\item{type}{
One of  "lasso", "lar", "forward.stagewise" or "stepwise". The names can
be abbreviated to any unique substring. Default is "lasso".
}
\item{trace}{
If TRUE, lars prints out its progress
}
\item{normalize}{
If TRUE, each variable is standardized to have unit L2 norm, otherwise
it is left alone. Default is TRUE.
}
\item{intercept}{
if TRUE, an intercept is included in the model (and not penalized),
otherwise no intercept is included. Default is TRUE.
}
\item{Gram}{
The X'X matrix; useful for repeated runs (bootstrap) where a large X'X
stays the same.
}
\item{eps}{
An effective zero
}
\item{max.steps}{
  Limit the number of steps taken; the default is \code{8 * min(m,
    n-intercept)}, with m the number of variables, and n the number of samples.
  For \code{type="lar"} or \code{type="stepwise"}, the maximum number of steps is
  \code{min(m,n-intercept)}. For \code{type="lasso"} and especially
  \code{type="forward.stagewise"}, there can be many more terms, because
  although no more than \code{min(m,n-intercept)} variables can be active during
  any step, variables are frequently droppped and added as the algorithm
  proceeds. Although the default usually guarantees that the algorithm
  has proceeded to the saturated fit, users should check.
}
\item{use.Gram}{
When the number m of variables is very large, i.e. larger than N, then
you may not want LARS to precompute the Gram matrix. Default is use.Gram=TRUE
}}
\value{
A "lars" object is returned, for which print, plot, predict, coef and summary
methods exist.
}
\details{
LARS is described in detail in Efron, Hastie, Johnstone and Tibshirani
(2002). With the "lasso" option, it computes the complete lasso
solution simultaneously for ALL values of the shrinkage parameter in
the same computational cost as a least squares fit. A "stepwise" option
has recently been added to LARS.
}
\references{
Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}; see also \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}.
Hastie, Tibshirani and Friedman (2002) Elements of Statistical
Learning, Springer, NY.
}
\author{Brad Efron and Trevor Hastie}
\seealso{
print, plot, summary and predict methods for lars, and  cv.lars
}
\examples{
data(diabetes)
par(mfrow=c(2,2))
attach(diabetes)
object <- lars(x,y)
plot(object)
object2 <- lars(x,y,type="lar")
plot(object2)
object3 <- lars(x,y,type="for") # Can use abbreviations
plot(object3)
detach(diabetes)
}
\keyword{regression}
% Converted by Sd2Rd version 1.21.


