\name{predict.lars}
\alias{predict.lars}
\alias{coef.lars}
\title{
Make predictions or extract coefficients from a fitted lars model
}
\description{
While lars() produces the entire path of solutions, predict.lars
allows one to extract a prediction at a particular point along the path.
}
\usage{
predict.lars(object, newx, s, type = c("fit", "coefficients"), mode = c("step", 
    "fraction", "norm"), ...)
coef.lars(object, ...)
}
\arguments{
\item{object}{
A fitted lars object
}
\item{newx}{
If type="fit", then newx should be the x values at which the fit is
required. If type="coefficients", then newx can be omitted.
}
\item{s}{
a value, or vector of values, indexing the path. Its values depends on the mode= argument. By
default (mode="step"), s should take on values between 0 and p (e.g.,
a step of 1.3 means .3 of the way between step 1 and 2.)
}
\item{type}{
If type="fit", predict returns the fitted values. If
type="coefficients", predict returns the coefficients.
Abbreviations allowed.
}
\item{mode}{
Mode="step" means the s= argument indexes the lars step number, and
the coefficients will be returned corresponding to the values
corresponding to step s. If mode="fraction", then s should be a number
between 0 and 1, and it refers to the ratio of the L1 norm of the
coefficient vector, relative to the norm at the full LS solution.
Mode="norm" means s refers to the L1 norm of the coefficient vector.
Abbreviations allowed.
}
\item{\dots}{Any arguments for \code{predict.lars} should work for \code{coef.lars}}
}
\value{
Either a vector/matrix of fitted values, or a vector/matrix of coefficients.
}
\details{
LARS is described in detail in Efron, Hastie, Johnstone and Tibshirani
(2002). With the "lasso" option, it computes the complete lasso
solution simultaneously for ALL values of the shrinkage parameter in
the same computational cost as a least squares fit.
}
\references{
Efron, Hastie, Johnstone and Tibshirani (2002)
Hastie, Tibshirani and Friedman (2002) Elements of Statistical
Learning, Springer, NY.
}
\author{Trevor Hastie}
\seealso{
print, plot,  lars, cv.lars
}
\examples{
data(diabetes)
attach(diabetes)
object <- lars(x,y,type="lasso")
### make predictions at the values in x, at each of the
### steps produced in object
fits <- predict.lars(object, x, type="fit")
### extract the coefficient vector with L1 norm=4.1
coef4.1 <- coef(object, s=4.1, mode="norm") # or
coef4.1 <- predict(object, s=4.1, type="coef", mode="norm")
detach(diabetes)
}
\keyword{regression}
\keyword{methods}
% Converted by Sd2Rd version 1.21.
