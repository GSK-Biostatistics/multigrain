#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calc_local_weights(NumericVector w, NumericMatrix G) {
  // w: initial local weights, length m
  // G: initial m x m transition matrix.
  int m = w.size();

  // Create names automatically: "H1", "H2", ..., "Hm".
  CharacterVector hyp_names(m);
  for (int i = 0; i < m; i++) {
    hyp_names[i] = "H" + std::to_string(i + 1);
  }

  // Number of intersection hypotheses is 2^m - 1.
  int n_intersections = std::pow(2, m) - 1;

  // ---- Prepare initial graph ----
  List initial_graph = List::create(
    Named("w") = clone(w),
    Named("G") = clone(G)
  );
  {
    NumericVector tmp = initial_graph["w"];
    tmp.attr("names") = hyp_names;
    initial_graph["w"] = tmp;
  }

  // ---- Compute the "parents" vector ----
  std::vector<int> parents;
  for (int i = 1; i <= m; i++) {
    int times = std::pow(2, i - 1);
    for (int j = 1; j <= times; j++) {
      parents.push_back(j);
    }
  }
  if (!parents.empty()) {
    parents.pop_back();
  }

  // ---- Compute the "delete" vector ----
  std::vector<int> rev_seq;
  for (int i = m; i >= 1; i--) {
    rev_seq.push_back(i);
  }
  std::vector<int> delete_vec;
  for (int i = 0; i < m; i++) {
    int times = std::pow(2, i);
    for (int j = 0; j < times; j++) {
      delete_vec.push_back(rev_seq[i]);
    }
  }
  if (!delete_vec.empty()) {
    delete_vec.pop_back();
  }

  // ---- Prepare storage for graphs ----
  int n_graphs = n_intersections;  // 2^m - 1
  std::vector<List> graphs(n_graphs);
  graphs[0] = initial_graph;

  // ---- Initialize matrix_weights ----
  NumericMatrix matrix_weights(n_intersections, m);
  for (int j = 0; j < m; j++) {
    matrix_weights(0, j) = w[j];
  }

  // ---- Loop over all intersections (except the initial one) ----
  int n_iter = parents.size(); // should be 2^m - 2
  for (int i = 0; i < n_iter; i++) {
    int parent_idx = parents[i] - 1; // convert to 0-indexed.
    List parent_graph = graphs[parent_idx];
    NumericVector parent_w = parent_graph["w"];
    NumericMatrix parent_G = parent_graph["G"];
    CharacterVector parent_names = parent_w.attr("names");

    // Determine deletion index: find index in parent's weights whose name
    // equals hyp_names[delete_vec[i]-1]
    std::string del_name = as<std::string>(hyp_names[delete_vec[i] - 1]);
    int del_index = -1;
    for (int j = 0; j < parent_w.size(); j++) {
      if (as<std::string>(parent_names[j]) == del_name) {
        del_index = j;
        break;
      }
    }
    if (del_index == -1) {
      stop("Deletion index not found.");
    }

    NumericVector init_parent_w = clone(parent_w);
    NumericMatrix init_parent_G = clone(parent_G);

    // Local copies for updating.
    NumericVector cur_w = clone(parent_w);
    NumericMatrix cur_G = clone(parent_G);

    // Build vector of indices to keep.
    std::vector<int> keep;
    for (int j = 0; j < cur_w.size(); j++) {
      if (j != del_index) {
        keep.push_back(j);
      }
    }

    // Update remaining hypotheses.
    for (size_t idx = 0; idx < keep.size(); idx++) {
      int h = keep[idx];
      cur_w[h] = init_parent_w[h] + init_parent_w[del_index] *
        init_parent_G(del_index, h);
      double denom = 1 - init_parent_G(h, del_index) *
        init_parent_G(del_index, h);
      for (size_t k_idx = 0; k_idx < keep.size(); k_idx++) {
        int k = keep[k_idx];
        if (h == k || denom <= 0) {
          cur_G(h, k) = 0;
        } else {
          cur_G(h, k) = (init_parent_G(h, k) + init_parent_G(h, del_index) *
            init_parent_G(del_index, k)) / denom;
        }
      }
    }

    // Create a new graph by removing the deleted hypothesis.
    int new_size = cur_w.size() - 1;
    NumericVector new_w(new_size);
    NumericMatrix new_G(new_size, new_size);
    CharacterVector new_names(new_size);
    int pos = 0;
    for (int j = 0; j < cur_w.size(); j++) {
      if (j == del_index) continue;
      new_w[pos] = cur_w[j];
      new_names[pos] = parent_names[j];
      pos++;
    }
    int r_cur = cur_G.nrow();
    int c_cur = cur_G.ncol();
    pos = 0;
    for (int row = 0; row < r_cur; row++) {
      if (row == del_index) continue;
      int pos2 = 0;
      for (int col = 0; col < c_cur; col++) {
        if (col == del_index) continue;
        new_G(pos, pos2) = cur_G(row, col);
        pos2++;
      }
      pos++;
    }

    // Build the new graph.
    List new_graph = List::create(
      Named("w") = new_w,
      Named("G") = new_G
    );
    {
      NumericVector tmp = new_graph["w"];
      tmp.attr("names") = new_names;
      new_graph["w"] = tmp;
    }
    new_graph.attr("class") = "initial_graph";
    graphs[i + 1] = new_graph;

    // Update matrix_weights.
    NumericVector new_weights(m, NA_REAL);
    for (int j = 0; j < m; j++) {
      std::string name = as<std::string>(hyp_names[j]);
      bool found = false;
      for (int k = 0; k < new_w.size(); k++) {
        if (as<std::string>(new_names[k]) == name) {
          new_weights[j] = new_w[k];
          found = true;
          break;
        }
      }
      if (!found) {
        new_weights[j] = NA_REAL;
      }
    }
    for (int j = 0; j < m; j++) {
      matrix_weights(i + 1, j) = new_weights[j];
    }
  }

  // Create matrix_intersections.
  int n_rows = matrix_weights.nrow();
  int n_cols = matrix_weights.ncol();
  LogicalMatrix matrix_intersections(n_rows, n_cols);
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      if (NumericVector::is_na(matrix_weights(i, j))) {
        matrix_intersections(i, j) = false;
        matrix_weights(i, j) = 0; // replace NA with 0
      } else {
        matrix_intersections(i, j) = true;
      }
    }
  }

  // Combine matrix_intersections and matrix_weights columnwise.
  NumericMatrix result(n_rows, 2 * n_cols);
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      result(i, j) = matrix_intersections(i, j) ? 1 : 0;
      result(i, j + n_cols) = matrix_weights(i, j);
    }
  }

  return result;
}
