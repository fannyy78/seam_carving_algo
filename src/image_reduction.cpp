#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector greedy_carving_cpp_raw(NumericVector im, IntegerVector dim, int n_seams) {
  int w = dim[0]; // largeur (x)
  int h = dim[1]; // hauteur (y)
  int c = dim[3]; // canaux

  std::vector<std::vector<std::vector<double>>> image(w, std::vector<std::vector<double>>(h, std::vector<double>(c)));

  int idx = 0;
  for (int ch = 0; ch < c; ++ch)
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        image[x][y][ch] = im[idx++];

  for (int seam_step = 0; seam_step < n_seams; ++seam_step) {
    std::vector<std::vector<double>> energy(h, std::vector<double>(w));
    for (int y = 0; y < h; ++y) {
      for (int x = 1; x < w - 1; ++x) {
        double e = 0;
        for (int ch = 0; ch < c; ++ch) {
          double diff = image[x + 1][y][ch] - image[x - 1][y][ch];
          e += diff * diff;
        }
        energy[y][x] = sqrt(e);
      }
    }

    std::vector<int> seam(h);
    seam[0] = w / 2;
    for (int y = 1; y < h; ++y) {
      int prev = seam[y - 1];
      double left = (prev > 0) ? energy[y][prev - 1] : 1e9;
      double mid  = energy[y][prev];
      double right = (prev < w - 1) ? energy[y][prev + 1] : 1e9;

      if (left < mid && left < right)
        seam[y] = prev - 1;
      else if (right < mid && right < left)
        seam[y] = prev + 1;
      else
        seam[y] = prev;
    }

    for (int y = 0; y < h; ++y) {
      int col = seam[y];
      for (int x = col; x < w - 1; ++x) {
        for (int ch = 0; ch < c; ++ch) {
          image[x][y][ch] = image[x + 1][y][ch];
        }
      }
    }

    w--; // nouvelle largeur
  }

  NumericVector out(w * h * c);
  int idx_out = 0;
  for (int ch = 0; ch < c; ++ch)
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        out[idx_out++] = image[x][y][ch];

  return out;
}

// [[Rcpp::export]]
SEXP greedy_carving_cpp(SEXP im, int n_seams) {
  NumericVector im_vec = as<NumericVector>(im);
  IntegerVector dims = im_vec.attr("dim");

  NumericVector result = greedy_carving_cpp_raw(im_vec, dims, n_seams);

  int new_x = dims[0] - n_seams;  // largeur = x
  int y = dims[1];                // hauteur = y
  int cc = dims[3];               // canaux

  IntegerVector new_dims = IntegerVector::create(new_x, y, 1, cc);
  result.attr("dim") = new_dims;
  result.attr("class") = "cimg";
  return result;
}
