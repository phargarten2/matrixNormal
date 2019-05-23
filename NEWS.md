# matrixNormal  0.0.1 
 * Documentation updates:  ADDED vignette, removed package R file.  Included the use and uniqueness of matrixNormal distribution. Also included an example that uses the package. Added datasets package to be imported in Documentation. 
 * is_symmetric()
 * rmatnorm(): Removed argument pre0.9_9994 that was passed from rmvnorm() function in mvtnorm library. This argument was introduced in mvtnorm library to fix a bug in version 0.9-9993, but matrixNormal uses the version at least 1.0.8. This argument simply is not needed. The origin of this argument was to fix bug in the mvtnorm library, but matrixNormal uses a version of mvtnorm library greater than update. This argument is just not needed, and if pre0.9_9994 is set to TRUE, nothing will happen. 
 * pmatnorm(), dmatnorm(), rmatnorm(): Added error if parameters of the matrix Normal Distribution -- M, U, or V -- contain any missing values. This makes the issue clearer to the user, but the function couldn't run anyway)
 * rmatnorm() now returns a matrix with rownames from U and the colnames from V. 

 * Minor Changes in Examples
      ** Changed order of examples in rmatnorm() for reproducibility
      ** Made clearer in using the dataset package in pmatnorm() example
      ** Removed installation of matrixcalc package in is.symmetric.matrix examples.


# matrixNormal 0.0.0.9000
* This is a new submission. 
* Added a `NEWS.md` file to track changes to the package.
* First Release of the Package
* Successfully passed windows check. 
 
