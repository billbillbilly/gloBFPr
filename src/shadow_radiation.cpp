#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

namespace {

const double EPS = 1e-12;

// Flattened ring geometry, extracted once from the R List so that the hot
// loops never touch Rcpp proxy objects or R's memory manager.
struct FlatRings {
  std::vector<double> x, y;      // all vertices, concatenated
  std::vector<int>    start, count;   // slice per ring
  std::vector<double> xmin, xmax, ymin, ymax;   // per-ring bbox
  std::vector<double> top;       // per-ring obstacle top elevation
  std::vector<char>   usable;    // finite top and >= 3 vertices
  int n = 0;
};

FlatRings flatten_rings(const List& rings,
                        const IntegerVector& ring_building,
                        const NumericVector& building_height,
                        const NumericVector& building_ground) {
  FlatRings fr;
  fr.n = rings.size();
  fr.start.resize(fr.n, 0);
  fr.count.resize(fr.n, 0);
  fr.xmin.resize(fr.n, 0.0);
  fr.xmax.resize(fr.n, 0.0);
  fr.ymin.resize(fr.n, 0.0);
  fr.ymax.resize(fr.n, 0.0);
  fr.top.resize(fr.n, 0.0);
  fr.usable.resize(fr.n, 0);

  std::size_t total = 0;
  for (int r = 0; r < fr.n; ++r) {
    NumericMatrix ring = rings[r];
    total += static_cast<std::size_t>(ring.nrow());
  }
  fr.x.reserve(total);
  fr.y.reserve(total);

  for (int r = 0; r < fr.n; ++r) {
    NumericMatrix ring = rings[r];
    int nr = ring.nrow();
    fr.start[r] = static_cast<int>(fr.x.size());
    fr.count[r] = nr;

    double bx0 =  std::numeric_limits<double>::infinity();
    double bx1 = -std::numeric_limits<double>::infinity();
    double by0 =  std::numeric_limits<double>::infinity();
    double by1 = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < nr; ++i) {
      double vx = ring(i, 0), vy = ring(i, 1);
      fr.x.push_back(vx);
      fr.y.push_back(vy);
      if (std::isfinite(vx) && std::isfinite(vy)) {
        bx0 = std::min(bx0, vx); bx1 = std::max(bx1, vx);
        by0 = std::min(by0, vy); by1 = std::max(by1, vy);
      }
    }
    fr.xmin[r] = bx0; fr.xmax[r] = bx1;
    fr.ymin[r] = by0; fr.ymax[r] = by1;

    int b = ring_building[r] - 1;
    double top = std::numeric_limits<double>::quiet_NaN();
    if (b >= 0 && b < building_height.size()) {
      top = building_ground[b] + building_height[b];
    }
    fr.top[r] = top;
    fr.usable[r] = (nr >= 3 && std::isfinite(top) && std::isfinite(bx0)) ? 1 : 0;
  }
  return fr;
}

// Shortest distance from a point to an axis-aligned bounding box (0 if inside).
inline double dist_to_bbox(double px, double py,
                           double x0, double x1, double y0, double y1) {
  double dx = (px < x0) ? (x0 - px) : ((px > x1) ? (px - x1) : 0.0);
  double dy = (py < y0) ? (y0 - py) : ((py > y1) ? (py - y1) : 0.0);
  return std::sqrt(dx * dx + dy * dy);
}

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

// [[Rcpp::export]]
NumericVector svf_cpp(NumericMatrix xy,
                      NumericVector azimuth,
                      List rings,
                      IntegerVector ring_building,
                      NumericVector building_height,
                      NumericVector building_ground,
                      NumericMatrix canopy_xy,
                      NumericVector canopy_height,
                      NumericVector canopy_ground,
                      NumericVector observer_ground,
                      NumericVector observer_height,
                      double canopy_cell_size,
                      double max_distance) {
  int n_points = xy.nrow();
  int n_az = azimuth.size();
  int n_canopy = canopy_xy.nrow();
  NumericVector out(n_points, NA_REAL);
  double half_cell = canopy_cell_size / 2.0;

  if (n_az == 0) return out;

  // --- Hoisted setup: done once, not once per (point, azimuth, ring) ---
  const FlatRings fr = flatten_rings(rings, ring_building,
                                     building_height, building_ground);
  std::vector<double> ux_all(n_az), uy_all(n_az);
  for (int a = 0; a < n_az; ++a) {
    double az = azimuth[a] * M_PI / 180.0;
    ux_all[a] = std::sin(az);
    uy_all[a] = std::cos(az);
  }
  const bool has_max_dist = std::isfinite(max_distance);
  // Raw pointers keep the parallel region free of R API access.
  const double* p_xy      = REAL(xy);
  const double* p_obs_grd = REAL(observer_ground);
  const double* p_obs_h   = REAL(observer_height);
  const int     n_obs_h   = observer_height.size();
  const double* p_can_xy  = n_canopy > 0 ? REAL(canopy_xy) : nullptr;
  const double* p_can_h   = n_canopy > 0 ? REAL(canopy_height) : nullptr;
  const double* p_can_grd = n_canopy > 0 ? REAL(canopy_ground) : nullptr;
  double* p_out           = REAL(out);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
  for (int p = 0; p < n_points; ++p) {
    double px = p_xy[p];
    double py = p_xy[p + n_points];
    if (!finite2(px, py)) continue;
    double obs_ground = p_obs_grd[p];
    if (!std::isfinite(obs_ground)) obs_ground = 0;
    double obs_h = (n_obs_h == 1) ? p_obs_h[0] : p_obs_h[p];
    if (!std::isfinite(obs_h)) obs_h = 0;
    double obs_z = obs_ground + obs_h;

    // --- Per-point ring culling (computed once, reused by every azimuth) ---
    // A ring can only matter if it rises above the observer AND lies within
    // max_distance.  In dense scenes this removes the vast majority of rings.
    // The bbox distance is cached because the per-azimuth horizon test needs it.
    std::vector<int>    cand;
    std::vector<double> cand_dmin;
    cand.reserve(64);
    cand_dmin.reserve(64);
    for (int r = 0; r < fr.n; ++r) {
      if (!fr.usable[r]) continue;
      if (fr.top[r] <= obs_z) continue;
      double d = dist_to_bbox(px, py, fr.xmin[r], fr.xmax[r],
                              fr.ymin[r], fr.ymax[r]);
      if (has_max_dist && d > max_distance) continue;
      cand.push_back(r);
      cand_dmin.push_back(std::max(d, EPS));
    }

    double svf_sum = 0.0;

    for (int a = 0; a < n_az; ++a) {
      double ux = ux_all[a];
      double uy = uy_all[a];
      // Track tan(beta) rather than beta: this keeps atan2 out of the hot loop
      // entirely, and cos^2(atan(m)) == 1/(1+m^2) recovers the SVF term exactly.
      double max_tan = 0.0;

      for (std::size_t ci = 0; ci < cand.size(); ++ci) {
        int r = cand[ci];
        double dh = fr.top[r] - obs_z;
        // Upper bound on this ring's tan(beta) is dh / nearest-bbox-distance.
        // If that cannot beat the horizon so far, the ring cannot contribute.
        if (dh / cand_dmin[ci] <= max_tan) continue;
        int off = fr.start[r];
        int nr  = fr.count[r];
        const double* rx = &fr.x[off];
        const double* ry = &fr.y[off];
        for (int i = 0; i < nr - 1; ++i) {
          double x1 = rx[i],     y1 = ry[i];
          double x2 = rx[i + 1], y2 = ry[i + 1];
          double sx = x2 - x1;
          double sy = y2 - y1;
          double denom = ux * sy - uy * sx;
          if (std::fabs(denom) <= EPS) continue;
          double qx = x1 - px;
          double qy = y1 - py;
          double t = (qx * sy - qy * sx) / denom;
          if (t <= EPS) continue;
          if (has_max_dist && t > max_distance) continue;
          double u = (qx * uy - qy * ux) / denom;
          if (u < -EPS || u > 1.0 + EPS) continue;
          double tn = dh / t;
          if (tn > max_tan) max_tan = tn;
        }
      }

      for (int i = 0; i < n_canopy; ++i) {
        double top = p_can_grd[i] + p_can_h[i];
        if (!std::isfinite(top) || top <= obs_z) continue;
        double cx = p_can_xy[i];
        double cy = p_can_xy[i + n_canopy];
        // Cheap radial reject before the slab test
        if (has_max_dist) {
          double ddx = cx - px, ddy = cy - py;
          if (ddx * ddx + ddy * ddy >
              (max_distance + half_cell) * (max_distance + half_cell)) continue;
        }
        // Skip if this cell cannot beat the horizon already established
        double dh_c = top - obs_z;
        {
          double d = dist_to_bbox(px, py, cx - half_cell, cx + half_cell,
                                  cy - half_cell, cy + half_cell);
          if (dh_c / std::max(d, EPS) <= max_tan) continue;
        }
        double xmin = cx - half_cell;
        double xmax = cx + half_cell;
        double ymin = cy - half_cell;
        double ymax = cy + half_cell;
        double tmin = -std::numeric_limits<double>::infinity();
        double tmax = std::numeric_limits<double>::infinity();

        if (std::fabs(ux) <= EPS) {
          if (px < xmin || px > xmax) continue;
        } else {
          double tx1 = (xmin - px) / ux;
          double tx2 = (xmax - px) / ux;
          tmin = std::max(tmin, std::min(tx1, tx2));
          tmax = std::min(tmax, std::max(tx1, tx2));
        }
        if (std::fabs(uy) <= EPS) {
          if (py < ymin || py > ymax) continue;
        } else {
          double ty1 = (ymin - py) / uy;
          double ty2 = (ymax - py) / uy;
          tmin = std::max(tmin, std::min(ty1, ty2));
          tmax = std::min(tmax, std::max(ty1, ty2));
        }
        if (tmax <= EPS || tmin > tmax) continue;
        double t = std::max(tmin, EPS);
        if (has_max_dist && t > max_distance) continue;
        double tn = dh_c / t;
        if (std::isfinite(tn) && tn > max_tan) max_tan = tn;
      }

      // cos^2(atan(max_tan)) — identical to the previous cos(beta)^2
      svf_sum += 1.0 / (1.0 + max_tan * max_tan);
    }
    p_out[p] = svf_sum / n_az;
  }
  return out;
}
