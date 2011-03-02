\name{summary.lars}
\alias{summary.lars}
\title{Summary method for lars objects}
\description{
Produce an anova-type summary for a lars object.
}
\usage{
\method{summary}{lars}(object, sigma2=NULL, ...)
}
\arguments{
  \item{object}{lars object}
  \item{sigma2}{optional variance measure (for p>n)}
  \item{\dots}{Additional arguments for summary generic}
}
\details{
  An anova summary is produced, with Df, RSS and Cp for each step.
  Df is tricky for some models, such as forward stagewise and stepwise,
  and is not likely to be accurate. When p>n, the user is responsible
  for supplying sigma2. 
}
\value{
  An anova object is returned, with rownames the step number, and with components:
  \item{Df }{Estimated degree of freedom}
  \item{Rss}{The Residual sum of Squares}
  \item{Cp}{The Cp statistic}
}
\references{
Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}; see also \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}.
Hastie, Tibshirani and Friedman (2002) Elements of Statistical
Learning, Springer, NY.
}
\author{Brad Efron and Trevor Hastie}
\seealso{
lars, and  print, plot,and predict methods for lars, and  cv.lars
}
\examples{
data(diabetes)
attach(diabetes)
object <- lars(x,y)
summary(object)
detach(diabetes)
}
\keyword{regression}
