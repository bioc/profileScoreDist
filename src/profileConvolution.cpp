#include "Rcpp.h"
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector profileConvolution (NumericVector x, NumericVector y,
                                  int LS, int US) {
  int i, j, m;
  m = x.size();
  NumericVector z(m);

  for(i = 0; i < m; i++) {
    z[i] = 0;
    for(j = 0; j < m; j++) {
      if ((i - j >= -LS) && (i - j <= US))
        z[i] += x[i - j + LS] * y[j];
    }
  }

  return z;
}
