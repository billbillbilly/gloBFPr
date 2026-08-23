#include <Rcpp.h>
using namespace Rcpp;

// For each building centroid, find the distance to the nearest green pixel
// within `radius` metres. All coordinates must be in the same projected CRS
// (metres). Returns NA for buildings with no green pixel within radius.
//
// [[Rcpp::export]]
NumericVector nearest_green_cpp(NumericMatrix centroids,
                                NumericMatrix green_xy,
                                double radius) {
  int n = centroids.nrow();
  int m = green_xy.nrow();
  NumericVector result(n, NA_REAL);
  double r2 = radius * radius;

  for (int i = 0; i < n; i++) {
    double cx = centroids(i, 0);
    double cy = centroids(i, 1);
    double min_d2 = R_PosInf;

    for (int j = 0; j < m; j++) {
      double dx = green_xy(j, 0) - cx;
      double dy = green_xy(j, 1) - cy;
      double d2 = dx * dx + dy * dy;
      if (d2 <= r2 && d2 < min_d2) {
        min_d2 = d2;
      }
    }

    if (min_d2 < R_PosInf) {
      result[i] = sqrt(min_d2);
    }
  }

  return result;
}
