## ----setup, include = FALSE----------------------------------------------
library(simesHotelling)
library(MASS)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(MASS)
X <- mvrnorm(50, rep(0, 100), diag(100)) 

## ------------------------------------------------------------------------
SHTest(X = X,
       samp.size = 30,
       iterations = 100,
       delta = rep(0, 100))

## ------------------------------------------------------------------------
SHTest(X = X,
       samp.size = 30,
       iterations = 100,
       delta = rep(0, 100))$`P-value`

set.seed(999)
SHTest(X = X,
       samp.size = 30,
       iterations = 100,
       delta = rep(0, 100))$`P-value`
set.seed(999)
SHTest(X = X,
       samp.size = 30,
       iterations = 100,
       delta = rep(0, 100))$`P-value`


## ------------------------------------------------------------------------
Y <- mvrnorm(50, rep(0, 100), diag(100)) 

SHTest(X = X,
       Y = Y, 
       samp.size = 30,
       iterations = 100,
       delta = rep(0, 100))

## ------------------------------------------------------------------------
ar.mat <- 0.8^abs(outer(1:100, 1:100, "-"))
### Generating data 
X <- mvrnorm(50, rep(0, 100), diag(100)) 
Y <- mvrnorm(79, rep(0, 100), ar.mat)


SHTest(X, 
       Y, 
       samp.size = 40, 
       iterations = 100, 
       equal.cov = FALSE,
       delta = rep(0, 100))

## ------------------------------------------------------------------------
X <- mvrnorm(50, rep(0, 200), diag(200)) 
Y <- mvrnorm(100, rep(0, 200), diag(200))

## ------------------------------------------------------------------------
L <- PrepareNonPara(X, Y,
               permutation.num = 100, 
               iter.num = 100, 
               samp.size = 50) 

## ------------------------------------------------------------------------
SHnonPara(L, method = 'SH')
SHnonPara(L, method = 'Thulin')

