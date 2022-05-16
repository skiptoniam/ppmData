#include <Rcpp.h>
using namespace Rcpp;
#include "deltri.h"
#include "dirtess3.h"

// define a polygon class that we'll use to catch the
constexpr std::size_t INVALID_INDEX = (std::numeric_limits<std::size_t>::max)();

typedef std::vector<std::vector<double>> polygon;

size_t next_half_edge(size_t edge) {
  return (edge % 3 == 2) ? edge - 2 : edge + 1;
}

size_t triangle_of_edge(size_t edge) {
  return edge / 3;
}

std::vector<size_t> edges_of_triangle(size_t tri){
  return std::vector<size_t>{3* tri, 3* tri + 1, 3 * tri + 2};
}

size_t cell_to_vertex(size_t tri, size_t lv){
  size_t v = 3* tri + lv;
  return v;
}

size_t find_vertex(size_t t, size_t v) {
  for(size_t lv=0; lv<3; ++lv) {
    if(cell_to_vertex(t,lv) == v) {
      return lv;
    }
  }
  size_t lv = INVALID_INDEX;
  return lv;
}

//finds 3 points for a given triangle id
std::vector<size_t> vertices_of_triangle(const deltri::deltri_cpp& del, const size_t tri_id){

  std::vector<size_t> edges = edges_of_triangle(tri_id);
  std::vector<std::size_t> find_tri_for_edge(3);

  for (int i = 0; i < 3; i++){
    find_tri_for_edge[i] = del.triangles[edges[i]];
  }

  return find_tri_for_edge;
}

deltri::site circumcenter_pts(const double ax, const double ay,
                                     const double bx, const double by,
                                     const double cx,const double cy) {

  const double dx = bx - ax;
  const double dy = by - ay;
  const double ex = cx - ax;
  const double ey = cy - ay;

  const double bl = dx * dx + dy * dy;
  const double cl = ex * ex + ey * ey;
  const double d = dx * ey - dy * ex;

  const double x = ax + (ey * bl - dy * cl) * 0.5 / d;
  const double y = ay + (dx * cl - ex * bl) * 0.5 / d;

  deltri::site res(x,y);
  return res;
}

deltri::site circumcenter_tri_id(const size_t tri_id, const deltri::deltri_cpp& del){

  // std::vector<double> coords = del.coords;
  std::vector<size_t> tri_point_ids = vertices_of_triangle(del, tri_id);

  NumericMatrix tri_coords(3,2);

  for (int ii=0; ii < 3; ii++){

    tri_coords(ii,0) = del.coords.at(tri_point_ids[ii]*2);
    tri_coords(ii,1) = del.coords.at(tri_point_ids[ii]*2 + 1);

  }

  deltri::site result = circumcenter_pts(tri_coords(0,0),tri_coords(0,1),
                                                tri_coords(1,0),tri_coords(1,1),
                                                tri_coords(2,0),tri_coords(2,1));

  return result;
}

// Compute exterior cell rays.
//<3 i love you <3


double dirtess_poly_area(NumericVector x, NumericVector y){
  // Initialize area
  double area = 0.0;

  // Calculate value of shoelace formula
  int n = x.size();
  int j = n - 1;
  for (int i = 0; i < n; i++){
    area += (x[j] + x[i]) * (y[j] - y[i]);
    j = i;  // j is previous vertex to i
  }

  // Return absolute value
  return std::abs(area / 2.0);
}


// get the edges for each cell of the Dirichlet polygon.
void dirtess_poly_i(size_t i, std::vector<size_t> inedges,
                    dirtess::polygon& poly_i,
                    const deltri::deltri_cpp& del){

  // polygon poly_i;
  size_t e0 = inedges[i];
  poly_i.resize(0);
  if (e0 == INVALID_INDEX) return; // coincident point
  deltri::site tmp_pt;
  size_t tri, e, e_next;

  e = e0;
  do {
    tri = triangle_of_edge(e);
    tmp_pt = circumcenter_tri_id(tri,del);
    poly_i.push(tmp_pt);
    // poly_i.push_back(tmp_pt.at(1));
    e_next = next_half_edge(e);
    e = del.halfedges[e_next];
  } while (e != e0 && e != INVALID_INDEX);

}

List dirtess_poly_list(size_t e, dirtess::polygon& poly_i){

  double tmp_x, tmp_y;
  NumericVector px;
  NumericVector py;
  bool onborder = false;
  int np = poly_i.size();

  // pass polygon by reference
  // dirtess_poly_i(e, del, poly_i);

  for(int i = 0; i < np; i++){
    tmp_x = poly_i.x(i);
    tmp_y = poly_i.y(i);
    px.push_back(tmp_x);
    py.push_back(tmp_y);
  }


  NumericMatrix polygon_coords_i(px.size(),2);
  polygon_coords_i(_,0) = px;
  polygon_coords_i(_,1) = py;
  double area_i = dirtess_poly_area(px,py);

  return List::create(_["poly.id"] = e,
                      _["poly.coords"] = polygon_coords_i,
                      _["area"] = area_i);

}


List dirtess_polygons(const deltri::deltri_cpp& del){//, std::vector<double> bbox){

  // How many polygons - should number of points (one per coordina)
  int np = del.coords.size()/2;
  Rcpp::List polylist(np);
  Rcpp::List polylist_clip(np);

  std::set<size_t> seen;
  dirtess::polygon poly_i;
  std::vector<size_t> inedges(del.coords.size()/2,INVALID_INDEX);

  // Compute an index from each point to an (arbitrary) incoming halfedge
  // Used to give the first neighbor of each point; for this reason,
  // on the hull we give priority to exterior halfedges
  size_t p;
  for (size_t e = 0;  e < del.halfedges.size(); e++) {
    p = del.triangles[e % 3 == 2 ? e - 2 : e + 1];
    if (del.halfedges[e] == INVALID_INDEX || inedges[p] == INVALID_INDEX){
      inedges[p] = e;
    }
  }


  // polyclip not working properly so lets' ignore for now.
  for(size_t i = 0; i<del.coords.size()/2; i++){
    dirtess_poly_i(i, inedges, poly_i, del);
    polylist[i] = dirtess_poly_list(i, poly_i);//not_clipped; look at testDelvor/dev/dirtess2.cpp for ideas.
  }

  return Rcpp::List::create(_["poly"]=polylist);

}


//[[Rcpp::export]]
List dirtess_cpp(std::vector<double> coords){//, std::vector<double> bbox) { //, Nullable<NumericVector> bbox = R_NilValue)

  // compute the delaunay triangulation and the use to calculate dirichlet tessellation
  deltri::deltri_cpp del(coords);

  // Compute the Dirichlet tessellation polygons - a bit slow with lists
  List dirtess_polys = dirtess_polygons(del);//##, bbox);

  return List::create(_["polygons"] = dirtess_polys);
}



