library(simesHotelling)
library(testthat)
library(MASS)

set.seed(999)
X <- mvrnorm(100, rep(0, 100), diag(100))
Y <- mvrnorm(100, rep(0, 100), diag(100))

test_that("Sample sizes", {  
  expect_warning(SHTest(X, Y), 
                 'Sample is larger then p, reverting to standard Hotelling test')
  expect_warning(SHTest(X, Y, samp.size = 200), 
                 'Sample is larger then p, reverting to standard Hotelling test')
  
})





test_that("Preparing permutation test", {  
  set.seed(999)
  X <- mvrnorm(100, rep(0, 100), diag(100))
  Y <- mvrnorm(100, rep(50, 99), diag(99))
  expect_error(PrepareNonPara(X, Y, permutation.num = 100, iter.num = 300, samp.size = 25), 
               'Samples are not of same dimesnions')
  Y <- mvrnorm(100, rep(50, 100), diag(100))
  expect_error(PrepareNonPara(X, Y, permutation.num = 100, iter.num = 300, samp.size = 500))  
  
})


test_that("Permutation test", {  
  set.seed(999)
  X <- mvrnorm(100, rep(0, 100), diag(100))
  Y <- mvrnorm(100, rep(50, 100), diag(100))
  L <- PrepareNonPara(X, Y, permutation.num = 100, iter.num = 300, samp.size = 25)
  expect_equal(SHnonPara(L, 'SH')$`Permutation P-value`, 1 / 100) 
  expect_equal(SHnonPara(L, 'Thullin')$`Permutation P-value`, 1 / 100)
  expect_warning(SHnonPara(L, 'TalTol'), 
                 'Method is not recognized, acceptable methods are SH and Thullin')
  
  
})
