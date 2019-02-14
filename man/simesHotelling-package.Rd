\name{simesHotelling-package}
\alias{simesHotelling-package}
\alias{simesHotelling}
\docType{package}
\title{
\packageTitle{simesHotelling}
}
\description{
\packageDescription{simesHotelling}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{simesHotelling}
\packageIndices{simesHotelling}
}
\author{
\packageAuthor{simesHotelling}

Maintainer: \packageMaintainer{simesHotelling}
}
\references{

}
\keyword{ Global null testing, Simes, Hotelling, Multivariate}
\seealso{

}
\examples{

library(MASS)
X <- mvrnorm(50, rep(0, 100), diag(100)) 
### One-sample global null test 
SHTest(X = X,
       samp.size = 30,
       iterations = 100,
       delta = rep(0, 100))
\dontrun{


Y <- mvrnorm(50, rep(0, 100), diag(100)) 
### Two-sample global null test 
SHTest(X = X,
       Y = Y, 
       samp.size = 30,
       iterations = 100,
       delta = rep(0, 100))
       
### Non-equal covariance two-sample global null test        
### Generating data 
X <- mvrnorm(50, rep(0, 100), diag(100)) 
Y <- mvrnorm(79, rep(0, 100), 2 * diag(100))


SHTest(X, 
       Y, 
       samp.size = 40, 
       iterations = 100, 
       equal.cov = FALSE,
       delta = rep(0, 100))
       
### Non parametric global null test 
L <- PrepareNonPara(X, Y,
               permutation.num = 100, 
               iter.num = 100, 
               samp.size = 50) 
SHnonPara(L, method = 'SH')
SHnonPara(L, method = 'Thulin')
} 
}
