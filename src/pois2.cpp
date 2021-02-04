#include <algorithm>
#include <math.h>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//' @rdname pois2
//' @export
// [[Rcpp::export]]
NumericVector dpois2(IntegerVector& x1,
                      IntegerVector& x2,
                       NumericVector& lambda1,
                       NumericVector& lambda2,
                       NumericVector& lambda3) {


  // if any of the inputs are empty, return empty vector.
  if (x1.length() == 0 || x2.length() == 0 || lambda1.length() == 0 ||
      lambda2.length() == 0 || lambda3.length() == 0 ){
    Rcpp::NumericVector empty_vec (0);
    return empty_vec;
  }

  const int N = std::max(x1.length(),
                         std::max(x2.length(),
                                  std::max(lambda1.length(),
                                           std::max(lambda2.length(),
                                                    lambda3.length()))));

  Rcpp::NumericVector res (N);


  for (int ii = 0; ii != N; ++ii){

    int x1_idx = ii % x1.length();
    int x2_idx = ii % x2.length();
    int l1_idx = ii % lambda1.length();
    int l2_idx = ii % lambda2.length();
    int l3_idx = ii % lambda3.length();

    if (Rcpp::IntegerVector::is_na(x1[x1_idx]) ||
        Rcpp::IntegerVector::is_na(x2[x2_idx]) ||
        Rcpp::NumericVector::is_na(lambda1[l1_idx]) ||
        Rcpp::NumericVector::is_na(lambda2[l2_idx]) ||
          Rcpp::NumericVector::is_na(lambda3[l3_idx])){
      res[ii] = NA_REAL;
      continue;
    } else {

      double lsum = lambda1[l1_idx] + lambda2[l2_idx] + lambda3[l3_idx];

      int minx = std::min(x1[x1_idx], x2[x2_idx]);

      double ss = 0;

      for (int jj = 0; jj != minx+1; ++jj){
        double enuml = std::pow(lambda1[l1_idx], x1[x1_idx]-jj) *
          std::pow(lambda2[l2_idx], x2[x2_idx]-jj) *
          std::pow(lambda3[l3_idx], jj);

        double log_denom = lgamma(x1[x1_idx]-jj + 1) +
          lgamma(x2[x2_idx]-jj + 1 ) + lgamma(jj + 1);

        ss += std::exp(std::log(enuml) - log_denom);

      }

      res[ii] = std::exp(-lsum) * ss;

    }



  }

  return(res);

}




