### Added for useRcpp 
#' @useDynLib simesHotelling, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
usethis::use_package('dplyr')
#' @importFrom dplyr "%>%"
NULL

## Simes Hoteling Test For equal covariance matrices and non equal covariance matrices
## Accepts two matrices returns Simes P-value
#' Conduct Simes Hotelling test
#' @param X A matrix of nx observations on p dimensions  
#' @param Y A matrix of ny observations on p dimensions
#' @param samp.size Number of dimensions (m) sampled at each iteration 
#' @param iterations Number of iterations (B) conducted  
#' @param equal.cov Logical, indicating if the equal.cov version of the test should be conducted 
#' @param delta Vector of expected differences between the mean vector of X and Y under the null 
#' @return A list including the dimensions used for analysis, the simes P-value, and the parameters used to conduct the test, including samp.size, iteration and the equal.cov flag. 
#' @export
SHTest <- function(X,
                   Y = NULL, 
                   samp.size  = NULL,
                   iterations = nrow(X) * log(nrow(X)),
                   equal.cov  = TRUE,
                   delta      = rep(0, ncol(X))) {
  ## Parameters of X 
  p.x        <- ncol(X)
  n.x        <- nrow(X)
  cov.mat.x  <- cov(X)
  mean.vec.x <- apply(X, 2, mean)
  two.sample <- !is.null(Y) 
  ### Checks 
  if (two.sample) {  ### Two samples checking 
    p.y <- ncol(Y)
    n.y <- nrow(Y)
    if (p.x != p.y) { 
      stop('Samples are not of same dimesnions')
    }
    samp.size <- ifelse(is.null(samp.size), floor((n.x + n.y) / 2), samp.size)
  }
  if (!two.sample) { 
    samp.size <- ifelse(is.null(samp.size), floor((n.x) / 2), samp.size)
  }
  ### General checks 
  if (!is.numeric(samp.size) 
      | !is.numeric(iterations) 
      | !is.logical(equal.cov) 
      | !is.numeric(delta)) { 
    stop('One of the parameters is not initalized correctly')  
  }
  if (samp.size >= p.x) { 
    warning('Sample is larger then p, reverting to standard Hotelling test') 
    samp.size <- p.x
    iteration <- 1 
  }
  if (!two.sample) {   ### One sample SH T2 test 
    print(paste('Conducting one SH test, sample size',
                samp.size, ', iteration', iterations))
    samp.size   <- ifelse(is.null(samp.size), floor(p.x / 2), samp.size)
    test.vec    <- mean.vec.x 
    cov.test    <- cov.mat.x / n.x 
    df.2        <- n.x - samp.size - 1
    t.to.f      <- (n.x - samp.size) / ((n.x - 1) * samp.size)
  }
  ### Two sample SH test (equal covariance assumption)
  if (equal.cov & two.sample) {
    print(paste('Conducting two sample equal covariance SH test, sample size',
                samp.size, ', iteration', iterations))
    mean.vec.y  <- apply(Y, 2, mean)
    cov.test    <- ((n.x - 1) * cov(X) + (n.y - 1) * cov(Y)) / (n.x + n.y - 2)
    test.vec    <- mean.vec.x - mean.vec.y
    df.2        <- n.x + n.y - samp.size - 1
    t.to.f      <- (n.x * n.y * (n.x + n.y - samp.size - 1)) / ((n.x + n.y - 2) * samp.size * (n.x + n.y)) 
  }
  if (equal.cov) { 
    t.statistic <-  RandomHoteling(test.vec,
                                   center = delta,
                                   cov = cov.test,
                                   n = iterations,
                                   samp_size = samp.size) 
    f.vec       <- t.to.f * t.statistic 
    p.vec       <- 1 - pf(t.to.f * t.statistic, samp.size, df.2) 
    return(list('P-value' = simesTest(p.vec),
                'F.vec'   = f.vec,
                'p.vec'   = p.vec))
  }
  ### Two sample SH test (no equal covariance assumption)
  if (!equal.cov & two.sample) {
    print(paste('Conducting non-equal covariance SH test, sample size',
                samp.size, ', iteration', iterations))
    SHNoEqual <- SimesHotellingNonEqual(X, 
                                        Y, 
                                        samp.size, 
                                        iterations, 
                                        cent = delta)
    return(SHNoEqual)
  } 
  else { 
    stop('Could not find test configurations')  
  }
}


## Simes Hoteling wrapper
## Input: X - matrix of multivariate normal, Y - matrix of multivariate norm we want to compare
##        B1 - number of repeatition
## Output: P- value of simes hotelling test
SimesHotellingNonEqual <- function(X , Y = NULL, samp.size, iterations, cent = rep(0, ncol(X))) {
  p        <- ncol(X)
  samp.mat <- t(replicate(iterations, sample(p, samp.size)))
  hot.test <- HotelingNonEqual(X, Y, samp.mat, cent)
  F.vec    <- (hot.test[ ,2] - samp.size  + 1) * hot.test[ ,1] / (hot.test[ ,2] * samp.size) 
  p.vec    <- 1 - pf(F.vec,
                     samp.size,
                     hot.test[,2] - samp.size - 1)
  return(list('P-value' = simesTest(p.vec),
              'F.vec'   = F.vec,
              'p.vec'   = p.vec,
              'df'      = hot.test[ ,2]))
}

## Simes test 
#' Conduct Simes test
#' @param p.vec A vector of p-values  
#' @return The simes p-value 
#' @export
simesTest <- function(p.vec)
{
  if (all(p.vec == 0)) {
    return(0)
  }
  n   <- length(p.vec)
  sort.vec<- sort(p.vec)
  return('Min Value' =  min(sort.vec * n / rep(1:n)))
}



# Non-parametric version of the Simes Hoteling test based on permutation tests. 
#' Conducts non parametric test, based on sampling of dimensions 
#' @param permute.list a list output from PrepareNonPara
#' @param method the type of test to conduct SH or Thullin 
#' @return A list containing the original statistic defined by method, permutation statistic vector, 
#' the p-value of the permutaiton test, method used in the test and number of permutations 
#' @export
SHnonPara <- function(permute.list, method = 'SH') {
  n.x <- nrow(permute.list$Original$X.org)
  n.y <- nrow(permute.list$Original$Y.org)
  p   <- ncol(permute.list$samp.mat)
  K   <- length(permute.list$Permutation)
  org.stat     <- SHnonParaIteration(permute.list$Original, permute.list$samp.mat)
  permute.vec  <- sapply(permute.list$Permutation, SHnonParaIteration, permute.list$samp.mat)
  if (method == 'SH') { 
    org.result      <- calcF(org.stat, 
                       n.x, 
                       n.y,
                       p) %>% simesTest()
    permute.result  <- calcF(permute.vec, 
                       n.x, 
                       n.y,
                       p) %>% apply(2, simesTest)
    permute.pval     <- mean(org.result > permute.result) + 1 / K
  }
  if (method == 'Thullin') { 
    org.result       <- org.stat %>% mean()
    permute.result   <- permute.vec %>% apply(2, mean)
    permute.pval     <- mean(org.result <= permute.result) + 1 / K
  }
  if (!(method %in% c('SH', 'Thullin'))) {
    warning('Method is not recognized, acceptable methods are SH and Thullin')
    org.result     <- org.stat 
    permute.result <- permute.vec
    permute.pval   <- NA
    method         <- 'Not recognized'
  }
  res.list    <- list('Original.statistic'     = org.result, 
                      'Permutation'            = permute.result, 
                      'Permutation P-value'    = permute.pval,
                      'Method'                 = method, 
                      'Number of permutations' = K)
  return(res.list) 
}

### Transform malhanobis distance to F statistic 

calcF <- function(stat, n.x, n.y, p) { 
 step1 <- n.x * n.y * stat / (n.x + n.y)
 step2 <- (n.x + n.y - p - 1) * step1 / ((n.x + n.y - 2) * p)
 return(1 - pf(step2, p, n.x + n.y - p - 1))  
}

### Permutes between X and Y 
PermuteTal <- function(X, Y) {
  new.mat <- rbind(X, Y)[sample(nrow(X) + nrow(Y)),]
  x.new <- new.mat[1:nrow(X),]
  y.new <- new.mat[(nrow(X) + 1):(nrow(Y) + nrow(X)),]
  return(list('X.perm' = x.new, 'Y.perm' = y.new))
}

### Prepares the non-parametric tests 
#' Prepares a list for non-parametric testing 
#' @param X A matrix of nx observations on p dimensions  
#' @param Y A matrix of ny observations on p dimensions
#' @param permutation.mat Number of permutations to conduct 
#' @param iter.num Number of samples of dimensions to conduct 
#' @param samp.size Number of dimensions to sample at each sampling 
#' @return A list which consists of Original X and Y matrices, Permutation a list length of permutation.num
#' consists of shuffled X and Y and samp.mat a matrix of iter.num on samp.size dimensions 
#' @export
PrepareNonPara <- function(X, Y, permutation.num, iter.num, samp.size) { 
  if (ncol(X) != ncol(Y)) { 
    stop('Samples are not of same dimesnions')
  }
  permute.mat <- replicate(permutation.num, PermuteTal(X, Y), simplify = FALSE)
  samp.mat    <- replicate(iter.num, sample(ncol(X), samp.size, replace = FALSE))
  return(list('Original'    = list('X.org' = X, 'Y.org' = Y),
              'Permutation' = permute.mat, 
              'samp.mat'    = t(samp.mat)))
}


## Conducts an single iteration of the SH non parametric test 
## Grouping of dimensions is according to samp.mat 

SHnonParaIteration <- function(list.mat, samp.mat) { 
  X   <- list.mat[[1]] 
  Y   <- list.mat[[2]]
  n.x <- nrow(X)
  p.x <- ncol(X)
  n.y <- nrow(Y) 
  mean.vec.x  <- colMeans(X)
  mean.vec.y  <- colMeans(Y)
  cov.mat.xy  <- ((n.x - 1) * cov(X) + (n.y - 1) * cov(Y)) / (n.x + n.y - 2)
  samp.size   <- ncol(samp.mat)
  t.statistic <- IndHoteling(mean.vec.x - mean.vec.y,
                             center = rep(0, p.x),
                             cov = cov.mat.xy,
                             samp_mat = samp.mat) 
  return(t.statistic) 
}
