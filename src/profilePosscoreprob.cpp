#include "Rcpp.h"
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector profilePosscoreprob (IntegerVector Scores, IntegerVector Srow,
                                   NumericVector prow) {
  int i, j, m, n;
  m = Scores.length();
  n = Srow.length();
  NumericVector fout(m);

  for (i = 0; i < m; i++) {
    fout[i] = 0;
    for (j = 0; j < n; j++) {
      if (Scores[i] == Srow[j])
        fout[i] += prow[j];
    }
  }

  return fout;
}
