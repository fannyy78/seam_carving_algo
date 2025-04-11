#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP seam_carving_dp(SEXP im, int n_seams) {
  NumericVector im_vec = as<NumericVector>(im);
  IntegerVector dims = im_vec.attr("dim");

  int w = dims[0]; // largeur (x)
  int h = dims[1]; // hauteur (y)
  int c = dims[3]; // canaux

  // Charger l'image dans tableau 3D : image[x][y][ch]
  std::vector<std::vector<std::vector<double>>> image(w, std::vector<std::vector<double>>(h, std::vector<double>(c)));
  int idx = 0;
  for (int ch = 0; ch < c; ++ch)
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        image[x][y][ch] = im_vec[idx++];

  for (int seam_step = 0; seam_step < n_seams; ++seam_step) {
    // 1. Énergie avec contraste renforcé
    std::vector<std::vector<double>> energy(w, std::vector<double>(h, 0.0));
    for (int x = 1; x < w - 1; ++x) {
      for (int y = 1; y < h - 1; ++y) {
        double e = 0.0;
        for (int ch = 0; ch < c; ++ch) {
          double dx = image[x + 1][y][ch] - image[x - 1][y][ch];
          double dy = image[x][y + 1][ch] - image[x][y - 1][ch];
          e += dx * dx + dy * dy;
        }
        energy[x][y] = std::sqrt(e) * 10.0; // accentuation visuelle
      }
    }

    // Optionnel : barrière sur les bords
    for (int y = 0; y < h; ++y) {
      energy[0][y] += 1e6;
      energy[w - 1][y] += 1e6;
    }

    // 2. Matrice cumulative
    std::vector<std::vector<double>> cumulative(w, std::vector<double>(h));
    for (int x = 0; x < w; ++x)
      cumulative[x][0] = energy[x][0];

    for (int y = 1; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        double left = (x > 0) ? cumulative[x - 1][y - 1] : 1e9;
        double up   = cumulative[x][y - 1];
        double right = (x < w - 1) ? cumulative[x + 1][y - 1] : 1e9;
        cumulative[x][y] = energy[x][y] + std::min({left, up, right});
      }
    }

    // 3. Retrouver la seam minimale
    std::vector<int> seam(h);
    int min_col = 0;
    for (int x = 1; x < w; ++x) {
      if (cumulative[x][h - 1] < cumulative[min_col][h - 1]) {
        min_col = x;
      }
    }
    seam[h - 1] = min_col;

    for (int y = h - 2; y >= 0; --y) {
      int prev = seam[y + 1];
      int best_x = prev;
      double best_val = cumulative[prev][y];

      if (prev > 0 && cumulative[prev - 1][y] < best_val) {
        best_x = prev - 1;
        best_val = cumulative[prev - 1][y];
      }
      if (prev < w - 1 && cumulative[prev + 1][y] < best_val) {
        best_x = prev + 1;
      }
      seam[y] = best_x;
    }

    // 4. Supprimer la seam
    for (int y = 0; y < h; ++y) {
      int col = seam[y];
      for (int x = col; x < w - 1; ++x) {
        for (int ch = 0; ch < c; ++ch) {
          image[x][y][ch] = image[x + 1][y][ch];
        }
      }
    }

    w--; // Réduction de la largeur
  }

  // 5. Reconvertir vers R (image finale)
  NumericVector out(w * h * c);
  int out_idx = 0;
  for (int ch = 0; ch < c; ++ch)
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        out[out_idx++] = image[x][y][ch];

  IntegerVector new_dims = IntegerVector::create(w, h, 1, c);
  out.attr("dim") = new_dims;
  out.attr("class") = "cimg";
  return out;
}
