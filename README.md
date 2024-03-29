# matrixNormal
Computes densities, probabilities, and random deviates of the Matrix Normal

--------------------------------------- Note: The following is a copy from the vignette ----------------------------------
---
title: "Introduction to Matrix Normal Package"
author: "Paul M. Hargarten"
date: "`r Sys.Date()`"
keywords: 
geometry: margin=1in
preamble: >
  \usepackage{indentfirst}
  \usepackage{amsmath}
output: 
   rmarkdown::html_vignette: 
       toc: true
vignette: >
  %\VignetteIndexEntry{Introduction to Matrix Normal Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
link-citations: true
citation_package: biblatex
bibliography: mnvignette.bib
csl: multidisciplinary-digital-publishing-institute.csl 
# csl: /Users/Shared/Zotero/styles/multidisciplinary-digital-publishing-institute.csl
---

<!---
Formatting Issues: matrix Normal Vignette
3/20/19
  The $ ... $ is not working
  The itemized list is not working 
  Bold all matrices, somehow
  Adjust the font -- Headers are too big. 
  
Then possibly send to Mom for edits? Does it make it sense? 
--->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,  # If TRUE, all output would be in the code chunk.
  results = "markup",
  comment = NA,
  prompt = TRUE,
  strip.white = TRUE,
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 90), # options for tidy to remove blank lines [blank = FALSE] and set the approximate line width to be 80.
  fig.show = "asis",
  fig.height = 4.5,  # inches
  fig.width = 4.5   # inches
)
library("formatR")
```
 
The `matrixNormal` package in the Comprehensive R Archive Network (CRAN) [@rcoreteamLanguageEnvironmentStatistical2018] consists of two different types of functions:  distribution and matrix ones. First, one can compute densities, probabilities and quantiles of the Matrix Normal Distribution. Second, one can perform useful but simple matrix operations like creating identity matrices and matrices of ones, calculating the trace of a matrix, and implementing the matrix operator `vec`().

## Installation
You can install the released version of `matrixNormal` from [CRAN](https://CRAN.R-project.org)[@rcoreteamLanguageEnvironmentStatistical2018] with `install.packages("matrixNormal")` and can load the package by:
```{r}
library(matrixNormal)
```

## Distribution functions

  The matrix normal distribution is a generalization of the multivariate normal distribution to matrix-valued random values. The parameters consist of a *n x p* mean matrix **M** and two positive-definite covariance matrices, one for rows **U** and another for columns **V**. Suppose that any *n x p* matrix $A \sim MatNorm_{n,p}(M, U, V)$. The mean of **A** is **M**, and the variance of `vec(A)` is the Kronecker product of **V** and **U**.

  The matrix normal distribution is a conjugate prior of the coefficients used in multivariate regression. Suppose  there are *p* predictors for *k* dependent variables. In univariate regression, only one dependent variable exists, so the conjugate prior for the *k* regression coefficients is the multivariate normal distribution. For multivariate regression, the extension of the conjugate prior of the coefficient matrix is the matrix normal distribution. This package also contains the PDF, CDF, and random number generation for the matrix normal distribution (*matnorm). See matrixNormal_Distribution file for more information.

  For instance, the USArrests dataset in the dataset package examines statistics in arrests per 100,000 residents for assault, murder, and rape in each of 50 states since 1973. Suppose that a researcher wants to determine whether the outcomes of assault and murder rates are associated with urban population and rape. The researcher then decides to use a Bayesian multiple linear regression and makes the assumption that the states are independent. However, the covariance between assault and murder is nonzero and needs to be taken into account. In fact, it has a correlation of 0.411 as given below. 
  
```{r}

library(datasets)
data(USArrests)
X <- cbind(USArrests$Assault, USArrests$Murder)
Y <- cbind(USArrests$UrbanPop, USArrests$Rape)
cor(Y)
```

We can assume that the outcome **Y** follows a matrix normal distribution with mean matrix *M*, which is the matrix product of the coefficient matrix, $\Psi$, times **X**. The covariance across states is assumed to be independent, and the covariance across the predictors be $\Sigma$. Or succinctly, with $k=2$ outcomes, 
$$ Y \sim MatNorm_{nx2}( M = X\Psi, U = I(n), V = \Sigma_2) $$

of coefficient matrix times **X**, and covariance across the predictors $\Phi$. The overall covariance matrix of **Y**, $\Sigma$, is a block diagonal matrix of $\Phi$. For instance, suppose that Y has the following distribution, and if we know the parameters, we can calculate its density using `dmatnorm()`. 

```{r}
# Y is n = 50 x p = 2 that follows a matrix normal with mean matrix M, which is product of
M <- (100 * toeplitz(50:1))[, 1:2]
dim(M)
head(M)
U <- I(50) # Covariance across states: Assumed to be independent
U[1:5, 1:5]
V <- cov(X) # Covariance across predictors
V

# Find the density if Y has the density with these arguments.
matrixNormal::dmatnorm(Y, M, U, V)
```

The coefficient matrix, $Psi$, for urban population and rape has dimensions *2 x 3* for the k = 2 outcomes and the p = 3 predictors. A semi-conjugate prior can be constructed on $\Psi$ to be a Matrix Normal distribution with mean $\Psi_0$ as a matrix of ones, covariance across the predictors as $X'X$, and with no covariance across states (due to independence). The conjugate prior for the covariance between the outcomes $\Sigma$ is the inverse-Wishart distribution with mean covariance $\Sigma_0$ and degrees of freedom $\nu$.  After defining the parameters, one random matrix is generated from this prior. 

```{r}
# Generate a random matrix from this prior.
# The prior mean of regression matrix
J(2, 3)
# The prior variance between rape and population
t(X) %*% X
# The prior variance between regression parameters
I(3)

# Random draw for prior would have these values
A <- matrixNormal::rmatnorm(M = J(2, 3), U = t(X) %*% X, V = I(3))
A

# Predicted Counts for y can be given as:
ceiling(rowSums(X %*% A))
```

However, the predicted counts do not have much meaning because the prior is uninformed from the data. We should use the posterior distribution to predict **Y**. Iranmanesh 2010  [@iranmaneshConditionalApplicationsMatrix2010] shows the posterior distributions to be
$$ \Psi | \Sigma, Y, X \sim MatNorm_{kxp}( \frac{\hat{\Psi} +\Psi_0}{2}, (2X^TX)^{-1}, \Sigma)$$
  and
$$ \Sigma | Y, X \sim W^{-1}_p(S^*, 2\nu+n)$$
  where:  
* $\hat{\Psi}$ is the maximum likelihood estimate (MLE) for the coefficient $\Psi$ matrix.
* *S* is the MLE for the covariance matrix $\Sigma$: $(Y-X\cdot\hat{\Psi})^T*(Y-X\cdot\hat{\Psi})$
* $S^{*}$ is the data-adjusted matrix of inverse Wishart, 
   $$ S^{*} = 2\Sigma_0 + S + (\hat{\Psi}-\Psi_0)'(2X^TX)^{-1}(\hat{\Psi}-\Psi_0) $$
* *c* is the dimension of $\Sigma$
* *n* is the sample size, the number of rows in **Y** and **X**. 

Additional details can be found in Iranmanesh et al. 2010 [@iranmaneshConditionalApplicationsMatrix2010] and Pocuca et al. 2019 [@Pocuca2019AssessingAV]. At any rate, the `matrixNormal` package can be used in Bayesian Multivariate Linear Regression. 

## Matrix functions
In the `matrixNormal` package, useful but simple matrix operations have been coded:

| Filename  | Description |
|---------- | ------------|
| is.symmetric.matrix   | Is a matrix square? symmetric? positive definite? or positive semi-definite? A tolerance is included here. |
|  Special_matrices | Creates the identity matrix *I* and matrix of 1's *J*.| 
|  tr               | Calculates the usual trace of a matrix| 
|  vec              | Stacks a matrix using matrix operator vec() and has option to keep names.| 
|  vech             | Stacks elements of a numeric symmetric matrix A in lower triangular only (using half vectorization, vech()).| 


For example, these functions can be applied to the following matrices.
```{r}
# Make a 3 x 3 Identity matrix
I(3)
# Make a 3 x 4 J matrix
J(3, 4)
# Make a 3 x 3 J matrix
J(3, 3)

# Calculate the trace of a J matrix
tr(J(3, 3))  # Should be 3

# Stack a matrix (used in distribution functions)
A <- matrix(c(1:4), nrow = 2, dimnames = list(NULL, c("A", "B")))
A
vec(A)

# Test if matrix is symmetric (used in distribution function)
is.symmetric.matrix(A)
```

## Conclusion
  Although other packages on CRAN have some matrixNormal functionality, this package provides a general approach in randomly sampling a matrix normal random variate. The `MBSP::matrix.normal`  , `matrixsampling::rmatrixnormal`, and `LaplacesDemon::rmatrixnorm` functions also randomly sample from the matrix normal distribution [@baiMBSPMultivariateBayesian2018; @laurentMatrixsamplingSimulationsMatrix2018; @statisticatLaplacesDemonCompleteEnvironment2018]. The function in [MBSP] (https://CRAN.R-project.org/package=MBSP/) uses Cholesky decomposition of individual matrices **U** and **V** [@baiMBSPMultivariateBayesian2018]. Similarly, the function in [matrixsampling] (https://CRAN.R-project.org/package=matrixsampling/) uses the Spectral decomposition of the individual matrices **U** and **V** [@laurentMatrixsamplingSimulationsMatrix2018]. Comparatively, the new function, `rmatnorm()`, in `matrixNormal` package is flexible in the decomposition of the covariance matrix, which is Kronecker product of **U** and **V**.  While you can simulate many samples using the `rmatrixnormal()` function from `matrixsampling` package, the new `rmatnorm`() function only generates one variate, but can generate many samples by placing `rmatnorm`() function inside a for-loop.  The [`LaplacesDemon`]  (https://cran.r-project.org/package=LaplacesDemon) package also has a density function, `dmatrixnorm`(), which calculates the log determinant of Cholesky decomposition of a positive definite matrix [@statisticatLaplacesDemonCompleteEnvironment2018]. However, the random matrix **A** that follows a matrix normal distribution does not need to be positive definite [@iranmaneshConditionalApplicationsMatrix2010]. There is no such restriction on the random matrix in `matrixNormal` package.
  
  In conclusion, the `matrixNormal` package collects all forms of the Matrix Normal Distribution in one place: calculating the PDF and CDF of the Matrix Normal distribution and simulating a random variate from this distribution. The package allows the users to be flexible in finding the random variate. Its main application in using this package is the Bayesian multivariate regression. 

## Computational Details 

This vignette is successfully processed using the following. 
```{r echo=FALSE}
# sessioninfo::session_info() # makes a mess! Instead

cat(" -- Session info ---------------------------------------------------")
sessioninfo::platform_info()
cat("--  Packages -------------------------------------------------------")
tmp.df <- sessioninfo::package_info(
  pkgs =c("LaplacesDemon", "MBSP", "matrixsampling", "matrixcalc"), dependencies = FALSE
  )
print(tmp.df)
```

 
## References
