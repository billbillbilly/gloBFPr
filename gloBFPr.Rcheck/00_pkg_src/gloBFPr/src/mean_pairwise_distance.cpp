#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double mean_pairwise_distance(Rcpp::NumericMatrix coords) {
  int n = coords.nrow();
  int dim = coords.ncol();
  double total_dist = 0.0;
  long count = 0;

  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double d = 0.0;
      for (int k = 0; k < dim; ++k) {
        double diff = coords(i, k) - coords(j, k);
        d += diff * diff;
      }
      total_dist += std::sqrt(d);
      count += 1;
    }
  }
  return total_dist / count;
}
