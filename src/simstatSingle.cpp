#include <Rcpp.h>
#include <random>
using namespace Rcpp;
using namespace std;

// Helper function that runs matrix-vector multiplication.
// Note: Prioritized over armadillo for easier class transformation
NumericVector mvMultiplication(NumericMatrix A, NumericVector B){
  int Arow = A.nrow();
  int Blen = B.length();

  NumericVector result(Arow);

  for (int i = 0; i < Arow; ++i){
    for (int j = 0; j < Blen; ++j){
      result(i) += A(i,j) * B(j);
    }
  }

  return result;
}

// [[Rcpp::export]]
double simstatSingle(int m, int K, NumericVector as, NumericVector bs, NumericMatrix L){
  // initialize max
  double Zmax = 0;

  // generate all random normals that we'll need
  NumericVector norms = rnorm(m*(K-1));

  // initialize vectors that contain indices for retrieving normals
  NumericVector starts(m);
  NumericVector ends(m);

  // generate indices for retrieving normals
  for (int j = 0; j < m; ++j){
    starts(j) = j * (K-1) + 1;
  }
  for (int k = 1; k <= m; ++k){
    ends(k-1) = k * (K-1);
  }

  // # get stats at first marker
  NumericVector Zlast = mvMultiplication(L, norms[Range(starts(0)-1, ends(0)-1)]);

  // update Zmax here
  if (max(abs(Zlast)) > Zmax){
    Zmax = max(abs(Zlast));
  }

  // loop through remaining markers
  for (int p = 2; p <= m; ++p){
    // generate stats at marker i
    NumericVector Znew = as(p-1) * Zlast +
      bs(p-1) * mvMultiplication(L, norms[Range(starts(p-1)-1, ends(p-1)-1)]);
    // check if we have a new max
    if (max(abs(Znew)) > Zmax){
      Zmax = max(abs(Znew));
    }
    // update Zlast for next marker
    Zlast = Znew;
  }

  return Zmax;
}
