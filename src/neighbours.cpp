// Per-cell summaries of the neighbours of each cell type within a radius, for one image: the distance to the
// nearest cell of each type (getDistances) or the number of cells of each type (getAbundances). A grid of
// bin side r replaces spatstat's closepairs() and tidyr's pivot_wider(), with the same results.
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// x, y: coordinates; type: 0-based cell-type codes (n_types of them); r: the radius (pairs with d <= r, the
// cell itself excluded, as closepairs()). mode 0: the nearest distance, r when a type has no cell within r;
// mode 1: the number of cells at 0 < d <= r. As distanceCalculator(): a cell with no neighbour at all is NA,
// and so is a type with no pair in the whole image (it is not a column of that image's pivot_wider()).
// [[Rcpp::export(.neighbourSummary)]]
NumericMatrix neighbourSummary(NumericVector x, NumericVector y, IntegerVector type, int n_types, double r, int mode) {
  const int n = x.size();
  NumericMatrix out(n, n_types);
  if (n == 0) return out;
  const double r2 = r * r;
  const double xmin = *std::min_element(x.begin(), x.end()), ymin = *std::min_element(y.begin(), y.end());
  const double xmax = *std::max_element(x.begin(), x.end()), ymax = *std::max_element(y.begin(), y.end());
  // bins of side >= r (every neighbour within r is in the 3 x 3 bins around a cell), at most about n of them
  const double area = std::max((xmax - xmin) * (ymax - ymin), 1e-300);
  const double side = std::max(r > 0 ? r : 1e-12, std::sqrt(area / n));
  const long long nbx = std::max(1LL, static_cast<long long>(std::floor((xmax - xmin) / side)) + 1);
  const long long nby = std::max(1LL, static_cast<long long>(std::floor((ymax - ymin) / side)) + 1);
  // cells sorted by bin (counting sort)
  std::vector<long long> bin(n);
  std::vector<int> start(nbx * nby + 1, 0), order(n);
  for (int i = 0; i < n; ++i) {
    long long bx = static_cast<long long>((x[i] - xmin) / side), by = static_cast<long long>((y[i] - ymin) / side);
    bin[i] = std::min(by, nby - 1) * nbx + std::min(bx, nbx - 1);
    ++start[bin[i] + 1];
  }
  for (long long b = 0; b < nbx * nby; ++b) start[b + 1] += start[b];
  std::vector<int> fill(start.begin(), start.end() - 1);
  for (int i = 0; i < n; ++i) order[fill[bin[i]]++] = i;

  std::vector<double> val(n_types);
  std::vector<char> present(n_types, 0);
  std::vector<char> any(n, 0);
  for (int i = 0; i < n; ++i) {
    std::fill(val.begin(), val.end(), mode == 0 ? r : 0.0);
    const long long bx = bin[i] % nbx, by = bin[i] / nbx;
    for (long long yy = std::max(0LL, by - 1); yy <= std::min(nby - 1, by + 1); ++yy)
      for (long long xx = std::max(0LL, bx - 1); xx <= std::min(nbx - 1, bx + 1); ++xx) {
        const long long b = yy * nbx + xx;
        for (int k = start[b]; k < start[b + 1]; ++k) {
          const int j = order[k];
          if (j == i) continue;
          const double dx = x[j] - x[i], dy = y[j] - y[i], d2 = dx * dx + dy * dy;
          if (d2 > r2) continue;
          const int t = type[j];
          any[i] = 1; present[t] = 1;
          if (mode == 0) { const double d = std::sqrt(d2); if (d < val[t]) val[t] = d; }
          else if (d2 > 0) val[t] += 1.0;
        }
      }
    for (int t = 0; t < n_types; ++t) out(i, t) = val[t];
  }
  for (int i = 0; i < n; ++i)
    for (int t = 0; t < n_types; ++t)
      if (!any[i] || !present[t]) out(i, t) = NA_REAL;
  return out;
}
