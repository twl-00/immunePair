#include <Rcpp.h>
using namespace Rcpp;

double chisq_2x2_cpp(int a, int b, int c, int d) {
  int row1 = a + b;
  int row2 = c + d;
  int col1 = a + c;
  int col2 = b + d;
  int n = row1 + row2;

  if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) {
    return NA_REAL;
  }

  double ad = static_cast<double>(a) * static_cast<double>(d);
  double bc = static_cast<double>(b) * static_cast<double>(c);
  double num = (ad - bc) * (ad - bc) * static_cast<double>(n);
  double den = static_cast<double>(row1) * static_cast<double>(row2) *
    static_cast<double>(col1) * static_cast<double>(col2);

  if (den == 0.0) {
    return NA_REAL;
  }

  double chi2_stat = num / den;
  return R::pchisq(chi2_stat, 1.0, false, false);
}

// [[Rcpp::export]]
List pairwise_chisq_cpp(NumericMatrix mat,
                        IntegerVector resp,
                        double min_prop = 0.05,
                        double max_prop = 0.95,
                        double chisq_cutoff = 0.01,
                        double diff_cutoff = 0.0,
                        double min_valid_prop = 0.5) {
  int n_gene = mat.nrow();
  int n_samp = mat.ncol();

  if (resp.size() != n_samp) {
    stop("resp length must equal the number of matrix columns.");
  }

  std::vector<int> g1_idx;
  std::vector<int> g2_idx;
  std::vector<int> nr0_vec, nr1_vec, r0_vec, r1_vec;
  std::vector<double> chisq_p_vec;

  g1_idx.reserve(100000);
  g2_idx.reserve(100000);
  nr0_vec.reserve(100000);
  nr1_vec.reserve(100000);
  r0_vec.reserve(100000);
  r1_vec.reserve(100000);
  chisq_p_vec.reserve(100000);

  for (int i = 0; i < n_gene - 1; ++i) {
    NumericVector gi = mat(i, _);

    for (int j = i + 1; j < n_gene; ++j) {
      NumericVector gj = mat(j, _);

      int nr0 = 0, nr1 = 0, r0 = 0, r1 = 0;
      int one_count = 0;
      int valid_n = 0;

      for (int s = 0; s < n_samp; ++s) {
        double v1 = gi[s];
        double v2 = gj[s];

        if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) {
          continue;
        }

        int resp_s = resp[s];
        double diff = v1 - v2;
        int cmp;

        if (diff > diff_cutoff) {
          cmp = 1;
        } else if (diff < -diff_cutoff) {
          cmp = 0;
        } else {
          continue;
        }

        if (cmp == 1) {
          one_count++;
        }
        valid_n++;

        if (resp_s == 0) {
          if (cmp == 0) {
            nr0++;
          } else {
            nr1++;
          }
        } else if (resp_s == 1) {
          if (cmp == 0) {
            r0++;
          } else {
            r1++;
          }
        }
      }

      if (valid_n == 0) {
        continue;
      }

      double valid_prop = static_cast<double>(valid_n) / static_cast<double>(n_samp);
      if (valid_prop < min_valid_prop) {
        continue;
      }

      double prop1 = static_cast<double>(one_count) / static_cast<double>(valid_n);
      if (prop1 <= min_prop || prop1 >= max_prop) {
        continue;
      }

      double p = chisq_2x2_cpp(nr0, nr1, r0, r1);
      if (R_IsNA(p) || p >= chisq_cutoff) {
        continue;
      }

      g1_idx.push_back(i + 1);
      g2_idx.push_back(j + 1);
      nr0_vec.push_back(nr0);
      nr1_vec.push_back(nr1);
      r0_vec.push_back(r0);
      r1_vec.push_back(r1);
      chisq_p_vec.push_back(p);
    }
  }

  return List::create(
    _["g1"] = g1_idx,
    _["g2"] = g2_idx,
    _["nr0"] = nr0_vec,
    _["nr1"] = nr1_vec,
    _["r0"] = r0_vec,
    _["r1"] = r1_vec,
    _["chisq_p"] = chisq_p_vec
  );
}
