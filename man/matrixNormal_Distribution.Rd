% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixNormal_Distribution.R
\name{matrixNormal_Distribution}
\alias{matrixNormal_Distribution}
\alias{dmatnorm}
\alias{pmatnorm}
\alias{rmatnorm}
\title{The Matrix Normal Distribution}
\usage{
dmatnorm(A, M, U, V, tol = .Machine$double.eps^0.5, use.log = TRUE)

pmatnorm(Lower = -Inf, Upper = Inf, M, U, V,
  tol = .Machine$double.eps^0.5, algorithm = mvtnorm::GenzBretz(), ...)

rmatnorm(M, U, V, tol = .Machine$double.eps^0.5, method = "chol")
}
\arguments{
\item{A}{The numeric n x p matrix that follows the matrix-normal.}

\item{M}{The mean n x p matrix that is numeric and real.}

\item{U}{The individual scale n x n real positive-definite matrix (rows)}

\item{V}{The parameter scale p x p  real positive-semidefinite matrix (columns)}

\item{tol}{A numeric tolerance level used to check if a matrix is symmetric; that is if the difference between the matrix and its transpose is between -tol and tol.}

\item{use.log}{Logical; if TRUE, densities d are given as log(d).}

\item{Lower}{The n x p matrix of lower limits for CDF}

\item{Upper}{The n x p matrix of upper limits for CDF}

\item{algorithm}{ an object of class \code{\link{GenzBretz}},
                    \code{\link{Miwa}} or \code{\link{TVPACK}}
                    specifying both the algorithm to be used as well as 
                    the associated hyper parameters.}

\item{...}{ additional parameters (currently given to \code{GenzBretz} for 
              backward compatibility issues). }

\item{method}{string specifying the matrix decomposition used to
   determine the matrix root of \code{sigma}.  Possible methods are
   eigenvalue decomposition (\code{"eigen"}, default),
   singular value decomposition (\code{"svd"}), and
   Cholesky decomposition (\code{"chol"}).  The
   Cholesky is typically fastest, not by much though.}
}
\description{
The density (dmatnorm), cumulative distribution function (CDF, pmatnorm), and generation of 1 random number from the matrix normal (rmatnorm) is produced:
\deqn{ A ~ MatNorm_{n,p}(M, U, V) }
}
\details{
Ideally, both scale matrices are positive-definite. However, they may not appear to be symmetric; you may want to increase the tolerance.

These functions rely heavily on this following property of matrix normal distribution. Let function `koch()` refer to the Kronecker product of a matrix. For a n x p matrix A, if  \eqn{ A \sim MatNorm(M, U, V)}, then \deqn{ vec(A) \sim MVN_{np} (M, Sigma = koch(U,V) ) }. Thus, we can find the probability that Lower < A < Upper by finding the CDF of vec(A), which is given in \code{\link[mvtnorm]{pmvnorm}} function in \pkg{mvtnorm}. Also, we can simulate 1 random matrix A from a matrix normal by sampling vec(A) from \pkg{mvtnorm}::\code{\link[mvtnorm]{rmvnorm}}() function. This matrix A takes the rownames from U and the colnames from V.
}
\examples{
#Data Used
A <- datasets::CO2[1:10, 4:5]
M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11) )
V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)
V  #Right covariance matrix (2 x 2), say the covariance between parameters.
U <- I(10) #Block of left-covariance matrix ( 84 x 84), say the covariance between subjects.

#PDF
dmatnorm(A, M, U, V )
dmatnorm(A, M, U, V, use.log = FALSE)

#Generating Probability Lower and Upper Bounds (They're matrices )
Lower <- matrix( rep(-1, 20), ncol = 2)
Upper <- matrix( rep(3, 20), ncol = 2)
Lower; Upper
#The probablity that a randomly chosen matrix A is between Lower and Upper
pmatnorm( Lower, Upper, M, U, V)
#CDF
pmatnorm( Lower = -Inf, Upper, M, U, V)
#entire domain = 1
pmatnorm( Lower = -Inf, Upper = Inf, M, U, V)

#Random generation
set.seed(123)
M <- cbind(rnorm(3, 435, 296), rnorm(3, 27, 11) )
U <- diag(1, 3)
V <- matrix(c(10, 5 ,5, 3), nrow = 2)
rmatnorm(M, U, V)
\dontrun{  #M has a different sample size than U; will return an error.
M <- cbind(rnorm(4, 435, 296), rnorm(4, 27, 11) )
rmatnorm(M, U, V)
}
}
\references{
Iranmanesh, Anis, M. Arashi, and S. M. M. Tabatabaey  On Conditional Applications of Matrix Variate Normal Distribution. \emph{Iranian Journal of Mathematical Sciences and Informatics} 5, no. 2. (November 1, 2010): 33-43. < https://doi.org/10.7508/ijmsi.2010.02.004 >
}
\seealso{
Other distribution: \code{\link{vec}}
}
\concept{distribution}
\keyword{distribution}
\keyword{matrixNormal}
