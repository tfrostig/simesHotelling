#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>


using namespace Rcpp ;

// [[Rcpp::export]]
arma::mat PowerMat(arma::mat tempMat, double power) { 
  int n = tempMat.n_rows; 
  int p = tempMat.n_cols;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < p; i++)  {
      tempMat(i,j) = pow(tempMat(i,j), power);
    }
  }
  return(tempMat);
}

// [[Rcpp::export]]
NumericVector colMeanRcpp(NumericMatrix X) { 
  int p = X.ncol();
  NumericVector meanVec(p);
  for (int i = 0; i < p; i++) { 
    meanVec(i) = mean(X(_,i)); 
  }
  return(meanVec); 
}


// [[Rcpp::export]] 
arma::mat CovFind(NumericMatrix X) { 
  int p = X.ncol();
  int n = X.nrow(); 
  NumericVector meanVec = colMeanRcpp(X);  
  NumericMatrix centerMat(n, p);  
  for (int i = 0; i < p; i++) {
    centerMat(_,i) = X(_,i) - meanVec(i);  
  }
  arma::mat tempMat = as<arma::mat>(centerMat); 
  arma::mat covMat = tempMat.t() * tempMat / (n - 1);
  return(covMat);
}  


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double MahalanobisArma(arma::vec x, arma::vec center, arma::mat cov) {
  arma::vec temp_vector = x - center;
  arma::vec sol_vector  = temp_vector.t() * inv(cov) * temp_vector;
  return(sol_vector(0));
}


// [[Rcpp::export]]
NumericVector RandomHoteling(arma::vec x, arma::vec center, arma::mat cov, const int n, int samp_size){
  NumericVector store_vec(n);
  IntegerVector indices_vector = seq_len(x.size());
  arma::uvec temp_samp(samp_size);
  for (int i=0; i < n; ++i) {
    IntegerVector temp_samp_2 = RcppArmadillo::sample(indices_vector, samp_size, false) - 1;
    Rcpp::wrap(temp_samp_2);
    temp_samp = as<arma::uvec>(Rcpp::wrap(temp_samp_2));
    store_vec[i] = MahalanobisArma(x.elem(temp_samp),
                                   center.elem(temp_samp),
                                   cov.submat(temp_samp, temp_samp)) ;
  }
  return(store_vec);
}

// [[Rcpp::export]] 
NumericVector IndHoteling(arma::vec x, arma::vec center, arma::mat cov, IntegerMatrix samp_mat){
  int samp_size = samp_mat.ncol();
  int samp_time = samp_mat.nrow();
  NumericVector store_vec(samp_time);
  arma::uvec temp_samp(samp_size);
  for (int i=0; i < samp_time; ++i) {
    IntegerVector temp_samp_2 = samp_mat(i,_) - 1;
    Rcpp::wrap(temp_samp_2);
    temp_samp = as<arma::uvec>(Rcpp::wrap(temp_samp_2));
    store_vec[i] = MahalanobisArma(x.elem(temp_samp),
                                   center.elem(temp_samp),
                                   cov.submat(temp_samp, temp_samp)) ;
  }
  return(store_vec);
}

// [[Rcpp::export]] 
arma::uvec CombineVectors(arma::uvec vec_1, int num) { 
  int n = vec_1.size();
  arma::uvec comb_vec(n + 1);
  for (int i = 0; i < n; i++){
    comb_vec(i) = vec_1(i) ; 
  }
  comb_vec(n) = num;
  return(comb_vec);
}
  
// [[Rcpp::export]] 
NumericVector HotelingFind(arma::uvec selected, arma::uvec not_select, arma::vec x, arma::vec center, arma::mat cov) {
  int n = not_select.size();  
  NumericVector store_vec(n);
  arma::uvec temp_vec(selected.size() + 1);
  for (int i=0; i < n; i++) {
    temp_vec = CombineVectors(selected, not_select(i)) - 1;
    store_vec[i] = MahalanobisArma(x.elem(temp_vec),
                                    center.elem(temp_vec),
                                    cov.submat(temp_vec, temp_vec)) ;
  }
  return(store_vec); 
}

// [[Rcpp::export]]
NumericVector CreateSol(arma::vec x, arma::vec center, arma::mat cov, IntegerVector group_vector) {
  IntegerVector uniq_group = unique(group_vector);
  int n = group_vector.size();
  IntegerVector indices_vector = seq_len(n) - 1;
  int group_num = max(group_vector);
  NumericVector sol_vector(group_num);
  for (int i=1; i < group_num + 1; i++) {
    IntegerVector temp_samp_2 = indices_vector[group_vector == i];
    if (temp_samp_2.isNULL()) {
      sol_vector[i - 1] = NA_REAL; 
    }
    else {
      Rcpp::wrap(temp_samp_2);
      arma::uvec temp_samp = as<arma::uvec>(Rcpp::wrap(temp_samp_2));
      sol_vector[i - 1] =  MahalanobisArma(x.elem(temp_samp),
                                           center.elem(temp_samp),
                                           cov.submat(temp_samp, temp_samp));
    }
  } 
  return(sol_vector);
}


// [[Rcpp::export]]
NumericMatrix pMatMaker(NumericMatrix statMat) {
  int p = statMat.ncol();
  int n = statMat.nrow();
  NumericMatrix Pmat(n, p);
  NumericVector permVec(p);  
  
  for (int i=0; i < n; i++) {
    for (int j=0; j < p; j++ ) {
      permVec = statMat(_,j); 
      Pmat(i,j) = 1 - mean(statMat(i,j) > permVec);
    }
  }
  return(Pmat);
}



// [[Rcpp::export]]
double FreedomDegreeFind(arma::mat covX, arma::mat covY, double nX, double nY) {
  arma::mat sXtild   = covX / nX;
  arma::mat sYtild   = covY / nY;
  arma::mat sComb    = sXtild + sYtild;
  arma::mat sCombInv = inv(sComb);
  double nominator = covX.n_cols + pow(covX.n_cols, 2);
  double denom1 = (1 / nX) * (arma::trace(PowerMat(sXtild * sCombInv, 2)) + pow(arma::trace(sXtild * sCombInv), 2));
  double denom2 = (1 / nY) * (arma::trace(PowerMat(sYtild * sCombInv, 2)) + pow(arma::trace(sYtild * sCombInv), 2));
  return(nominator / (denom1 + denom2));
}



// [[Rcpp::export]]
double HotelingCalc(NumericMatrix X, NumericMatrix Y, arma::vec center) {
  int nx = X.nrow(); 
  int ny = X.nrow();
  arma::mat sX = CovFind(X);  
  arma::mat sY = CovFind(Y);  
  arma::mat combS = sX / nx + sY / ny; 
  NumericVector tempVec = colMeanRcpp(X) - colMeanRcpp(Y); 
  double tStar = MahalanobisArma(tempVec, center, combS);
  return(tStar);
}

// [[Rcpp::export]]
NumericMatrix HotelingNonEqual(NumericMatrix X, NumericMatrix Y, IntegerMatrix samp_mat, arma::vec center ){
  int samp_size = samp_mat.ncol();
  int samp_time = samp_mat.nrow();
  int nx = X.nrow(); 
  int ny = Y.nrow();
  arma::mat covX = CovFind(X); 
  arma::mat covY = CovFind(Y); 
  arma::mat cov = covX / nx + covY / ny;
  arma::uvec temp_samp(samp_size);
  arma::vec mean_diff = colMeanRcpp(X) - colMeanRcpp(Y); 
  NumericMatrix store_mat(samp_time, 2);
  for (int i=0; i < samp_time; ++i) {
    IntegerVector temp_samp_2 = samp_mat(i,_) - 1;
    Rcpp::wrap(temp_samp_2);
    temp_samp = as<arma::uvec>(Rcpp::wrap(temp_samp_2));
    store_mat(i,0) = MahalanobisArma(mean_diff.elem(temp_samp), 
                                   center.elem(temp_samp),
                                   cov.submat(temp_samp, temp_samp));
    store_mat(i,1) = FreedomDegreeFind(covX.submat(temp_samp, temp_samp), 
                                       covY.submat(temp_samp, temp_samp),
                                       nx, 
                                       ny);
  } 
  return(store_mat);
}
