#'Generating Special Matrices
#'
#'@family matrix
#'@keywords matrix identity ones
#'
#'@name special.matrix
#'
#'@description Creates Identity Matrix I and Matrix of Ones J.
#'
#'@details
#'Create an Identity Matrix where the number of columns is n.  This is a diagonal matrix with all equal to one (1). An identity matrix is usually written as I. To make an identity matrix with r rows and columns, use \code{I}.
#'
#'A J matrix is a general matrix of any number of rows and columns, but in which all elements in the matrix are equal to one (1). \code{J} will make a n x m J matrix, given the number or rows, n, and number of columns, m. Names of rows and columns (dimnames) are included.
#'
#'@param n number of rows in I or J.
#'@param m number of columns in J. Default: same as number of rows.
#'
#'@examples
#' #To create an identity matrix of order 12
#' I(2)
#' #To make a matrix of 6 rows and 10 columns of all ones
#' J(6,10)
#' #To make a matrix of unity, dimensions 6 x 6 .
#' J(6)

#'@rdname  special.matrix
#'@export I
I <- function(n){
  #Identity Matrix where number of columns is n.
  D <- diag(1, n)
  dimnames(D) <- list(1:n, 1:n)
  return(D)
}

#'@rdname  special.matrix
#'@export J
J<- function(n, m=n){
  #A matrix of ones with n rows and m columns.
  matrix(c(1),nrow=n,ncol=m, dimnames = list(1:n, 1:m))
}



