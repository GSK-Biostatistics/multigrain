// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>


// Helper function to generate a random number using a scaled Cauchy distribution
double randcauchy(double min, double max) {
  double mi = 0.0, t = 1.0, band = 10.0;
  double limit_inf = mi - (band * 0.5);
  double limit_sup = mi + (band * 0.5);

  double cauchy_mit;
  do {
    double na_unif = R::runif(0, 1);
    cauchy_mit = t * tan((na_unif - 0.5) * M_PI) + mi;
  } while ( (cauchy_mit < limit_inf) || (cauchy_mit > limit_sup) );

  if (cauchy_mit < 0) {
    cauchy_mit = -cauchy_mit;
  } else {
    cauchy_mit = cauchy_mit + (band * 0.5);
  }

  double valor = cauchy_mit / band;
  return min + (max - min) * valor;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix esch_population_rcpp(int popSize,
                                         int nBits,
                                         Rcpp::NumericVector lower,
                                         Rcpp::NumericVector upper) {
  Rcpp::RNGScope scope;
  Rcpp::NumericMatrix population(popSize, nBits);
  for(int i = 0; i < popSize; ++i) {
    for(int j = 0; j < nBits; ++j) {
      population(i, j) = randcauchy(lower[j], upper[j]);
    }
  }
  return population;
}

// Deprecated: included only for benchmarking purposes
// [[Rcpp::export]]
Rcpp::NumericVector esch_mutation_rcpp(Rcpp::NumericVector parent,
                                       Rcpp::NumericVector lower,
                                       Rcpp::NumericVector upper) {
  Rcpp::RNGScope scope;
  int n_params = parent.size();
  int mutation_point = R::runif(0, n_params);

  Rcpp::NumericVector mutated_parent = Rcpp::clone(parent);
  mutated_parent[mutation_point] = randcauchy(lower[mutation_point], upper[mutation_point]);
  return mutated_parent;
}
