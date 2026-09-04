#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector allocate_population_cpp(NumericVector pop_total,
                                      IntegerVector cell_group,
                                      NumericVector volume,
                                      NumericVector area) {
  R_xlen_t n = pop_total.size();
  if (cell_group.size() != n || volume.size() != n || area.size() != n) {
    stop("Population, cell group, volume, and area vectors must have equal length.");
  }

  NumericVector out(n, NA_REAL);
  if (n == 0) return out;

  int max_group = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (cell_group[i] != NA_INTEGER && cell_group[i] > max_group) {
      max_group = cell_group[i];
    }
  }
  if (max_group <= 0) return out;

  NumericVector volume_sum(max_group);
  NumericVector area_sum(max_group);
  IntegerVector count(max_group);
  NumericVector group_pop(max_group, NA_REAL);

  for (R_xlen_t i = 0; i < n; ++i) {
    int group = cell_group[i];
    if (group == NA_INTEGER || group <= 0) continue;
    int index = group - 1;
    if (NumericVector::is_na(pop_total[i])) continue;

    if (NumericVector::is_na(group_pop[index])) {
      group_pop[index] = pop_total[i];
    }
    count[index] += 1;

    if (!NumericVector::is_na(volume[i]) && volume[i] > 0) {
      volume_sum[index] += volume[i];
    }
    if (!NumericVector::is_na(area[i]) && area[i] > 0) {
      area_sum[index] += area[i];
    }
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    int group = cell_group[i];
    if (group == NA_INTEGER || group <= 0) continue;
    int index = group - 1;
    if (NumericVector::is_na(group_pop[index]) || count[index] == 0) continue;

    double denominator = volume_sum[index];
    double weight = (!NumericVector::is_na(volume[i]) && volume[i] > 0) ? volume[i] : 0.0;

    if (denominator <= 0) {
      denominator = area_sum[index];
      weight = (!NumericVector::is_na(area[i]) && area[i] > 0) ? area[i] : 0.0;
    }

    if (denominator <= 0) {
      out[i] = group_pop[index] / static_cast<double>(count[index]);
    } else {
      out[i] = group_pop[index] * weight / denominator;
    }
  }

  return out;
}
