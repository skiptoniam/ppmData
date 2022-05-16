#include <Rcpp.h>
using namespace Rcpp;
#include"deltri.hpp"

//[[Rcpp::export]]
List deltri_cpp(std::vector<double> coords) {

  deltri::deltri_cpp del(coords);

  return List::create(_["coords"] = del.coords,
                      _["triangles"] = del.triangles,
                      _["halfedges"] = del.halfedges,
                      _["convexhull"] = del.hull,
                      _["convexhull.area"] = del.get_hull_area());

}


