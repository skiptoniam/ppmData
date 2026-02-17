#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List matrices_to_polygons(List matrices) {
  List polygons(matrices.size());

  for (int i = 0; i < matrices.size(); ++i) {
    NumericMatrix coords = matrices[i];
    List ring(1);
    ring[0] = coords;
    polygons[i] = ring;
  }

  return polygons;
}
