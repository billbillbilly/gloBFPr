#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <limits>
using namespace Rcpp;

namespace {

const double EPS = 1e-12;

bool finite2(double x, double y) {
  return std::isfinite(x) && std::isfinite(y);
}

bool point_on_segment(double px, double py, double ax, double ay, double bx, double by) {
  double cross = (px - ax) * (by - ay) - (py - ay) * (bx - ax);
  double len = std::hypot(bx - ax, by - ay);
  if (len <= EPS) return std::hypot(px - ax, py - ay) <= EPS;
  if (std::fabs(cross) > EPS * std::max(1.0, len)) return false;
  double dot = (px - ax) * (px - bx) + (py - ay) * (py - by);
  return dot <= EPS;
}

bool point_in_ring(const NumericMatrix& ring, double px, double py) {
  int n = ring.nrow();
  if (n < 3 || !finite2(px, py)) return false;
  bool inside = false;
  for (int i = 0, j = n - 1; i < n; j = i++) {
    double xi = ring(i, 0), yi = ring(i, 1);
    double xj = ring(j, 0), yj = ring(j, 1);
    if (point_on_segment(px, py, xi, yi, xj, yj)) return true;
    bool crosses = ((yi > py) != (yj > py));
    if (crosses) {
      double x_intersect = (xj - xi) * (py - yi) / (yj - yi) + xi;
      if (px < x_intersect) inside = !inside;
    }
  }
  return inside;
}

bool point_in_quad(double px, double py,
                   double x1, double y1, double x2, double y2,
                   double x3, double y3, double x4, double y4) {
  double xs[4] = {x1, x2, x3, x4};
  double ys[4] = {y1, y2, y3, y4};
  bool has_pos = false;
  bool has_neg = false;
  for (int i = 0; i < 4; ++i) {
    int j = (i + 1) % 4;
    if (point_on_segment(px, py, xs[i], ys[i], xs[j], ys[j])) return true;
    double cross = (xs[j] - xs[i]) * (py - ys[i]) -
      (ys[j] - ys[i]) * (px - xs[i]);
    if (cross > EPS) has_pos = true;
    if (cross < -EPS) has_neg = true;
    if (has_pos && has_neg) return false;
  }
  return true;
}

double distance_to_segment(double px, double py, double ax, double ay, double bx, double by) {
  double vx = bx - ax;
  double vy = by - ay;
  double len2 = vx * vx + vy * vy;
  if (len2 <= EPS) return std::hypot(px - ax, py - ay);
  double t = ((px - ax) * vx + (py - ay) * vy) / len2;
  t = std::max(0.0, std::min(1.0, t));
  double cx = ax + t * vx;
  double cy = ay + t * vy;
  return std::hypot(px - cx, py - cy);
}

}

// [[Rcpp::export]]
NumericVector building_shadow_height_cpp(NumericMatrix xy,
                                         List rings,
                                         IntegerVector ring_building,
                                         IntegerVector ring_part,
                                         LogicalVector ring_is_hole,
                                         NumericVector heights,
                                         double dx_unit,
                                         double dy_unit,
                                         double tan_elev) {
  int n_points = xy.nrow();
  int n_rings = rings.size();
  int n_buildings = heights.size();
  NumericVector out(n_points, NA_REAL);

  if (!std::isfinite(dx_unit) || !std::isfinite(dy_unit) || !std::isfinite(tan_elev)) {
    return out;
  }

  for (int p = 0; p < n_points; ++p) {
    double px = xy(p, 0);
    double py = xy(p, 1);
    if (!finite2(px, py)) continue;

    double best = NA_REAL;
    for (int b = 1; b <= n_buildings; ++b) {
      double building_height = heights[b - 1];
      double dx = dx_unit * building_height;
      double dy = dy_unit * building_height;
      bool in_source = false;
      bool in_shifted = false;
      bool in_side = false;
      double min_dist = std::numeric_limits<double>::infinity();

      for (int r = 0; r < n_rings; ++r) {
        if (ring_building[r] != b) continue;
        NumericMatrix ring = rings[r];
        bool is_hole = ring_is_hole[r];

        bool inside_ring = point_in_ring(ring, px, py);
        bool inside_shifted_ring = point_in_ring(ring, px - dx, py - dy);
        if (!is_hole) {
          in_source = in_source || inside_ring;
          in_shifted = in_shifted || inside_shifted_ring;
        } else {
          if (inside_ring) in_source = false;
          if (inside_shifted_ring) in_shifted = false;
        }

        int nr = ring.nrow();
        for (int i = 0; i < nr - 1; ++i) {
          double x1 = ring(i, 0), y1 = ring(i, 1);
          double x2 = ring(i + 1, 0), y2 = ring(i + 1, 1);
          if (point_in_quad(px, py, x1, y1, x2, y2, x2 + dx, y2 + dy, x1 + dx, y1 + dy)) {
            in_side = true;
          }
          double d = distance_to_segment(px, py, x1, y1, x2, y2);
          if (d < min_dist) min_dist = d;
        }
      }

      bool shadowed = (in_shifted || in_side) && !in_source;
      if (!shadowed || !std::isfinite(min_dist)) continue;
      double h = building_height - min_dist * tan_elev;
      if (!std::isfinite(h) || h < 0) continue;
      if (NumericVector::is_na(best) || h > best) best = h;
    }
    out[p] = best;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector canopy_shadow_height_cpp(NumericMatrix xy,
                                       NumericMatrix canopy_xy,
                                       NumericVector canopy_height,
                                       NumericVector canopy_ground,
                                       NumericVector target_ground,
                                       double dx_unit,
                                       double dy_unit,
                                       double tan_elev,
                                       double half_cell) {
  int n_points = xy.nrow();
  int n_canopy = canopy_xy.nrow();
  NumericVector out(n_points, NA_REAL);

  if (!std::isfinite(dx_unit) || !std::isfinite(dy_unit) ||
      !std::isfinite(tan_elev) || !std::isfinite(half_cell)) {
    return out;
  }

  for (int i = 0; i < n_canopy; ++i) {
    double h_canopy = canopy_height[i];
    if (!std::isfinite(h_canopy) || h_canopy <= 0) continue;
    double sx = canopy_xy(i, 0);
    double sy = canopy_xy(i, 1);
    double vx = dx_unit * h_canopy;
    double vy = dy_unit * h_canopy;
    double len2 = vx * vx + vy * vy;
    if (!std::isfinite(len2) || len2 <= EPS) continue;
    double source_top = canopy_ground[i] + h_canopy;

    for (int p = 0; p < n_points; ++p) {
      double px = xy(p, 0);
      double py = xy(p, 1);
      if (!finite2(px, py)) continue;
      double relx = px - sx;
      double rely = py - sy;
      double t = (relx * vx + rely * vy) / len2;
      if (t < 0 || t > 1) continue;
      double cx = sx + t * vx;
      double cy = sy + t * vy;
      double cross_dist = std::hypot(px - cx, py - cy);
      if (cross_dist > half_cell) continue;
      double along_dist = std::hypot(relx, rely);
      double target_z = target_ground[p];
      if (!std::isfinite(target_z)) target_z = 0;
      double h = source_top - along_dist * tan_elev - target_z;
      if (!std::isfinite(h) || h < 0) continue;
      if (NumericVector::is_na(out[p]) || h > out[p]) out[p] = h;
    }
  }
  return out;
}
