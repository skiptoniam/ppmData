#include <Rcpp.h>
using namespace Rcpp;

// Helper function to check if a point is inside the bounding box
bool is_inside(double x, double y, double xmin, double xmax, double ymin, double ymax) {
  return (x >= xmin && x <= xmax && y >= ymin && y <= ymax);
}

// Helper function to find intersection point of a line segment with bounding box
bool find_intersection(double x1, double y1, double x2, double y2,
                       double xmin, double xmax, double ymin, double ymax,
                       double& x_int, double& y_int) {
  // Calculate line parameters
  double dx = x2 - x1;
  double dy = y2 - y1;

  if (dx != 0) {
    double slope = dy / dx;

    // Left boundary
    if (x1 < xmin && x2 >= xmin) {
      y_int = y1 + slope * (xmin - x1);
      if (y_int >= ymin && y_int <= ymax) {
        x_int = xmin;
        return true;
      }
    }
    // Right boundary
    if (x1 > xmax && x2 <= xmax) {
      y_int = y1 + slope * (xmax - x1);
      if (y_int >= ymin && y_int <= ymax) {
        x_int = xmax;
        return true;
      }
    }
  }

  if (dy != 0) {
    double inv_slope = dx / dy;

    // Bottom boundary
    if (y1 < ymin && y2 >= ymin) {
      x_int = x1 + inv_slope * (ymin - y1);
      if (x_int >= xmin && x_int <= xmax) {
        y_int = ymin;
        return true;
      }
    }
    // Top boundary
    if (y1 > ymax && y2 <= ymax) {
      x_int = x1 + inv_slope * (ymax - y1);
      if (x_int >= xmin && x_int <= xmax) {
        y_int = ymax;
        return true;
      }
    }
  }

  return false;
}

// [[Rcpp::export]]
List clip_polygon_to_bbox(NumericMatrix coords, double xmin, double xmax, double ymin, double ymax) {
  std::vector<double> x_clipped;
  std::vector<double> y_clipped;
  int n = coords.nrow();

  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n; // Next vertex index to form an edge

    double x1 = coords(i, 0);
    double y1 = coords(i, 1);
    double x2 = coords(j, 0);
    double y2 = coords(j, 1);

    bool inside1 = is_inside(x1, y1, xmin, xmax, ymin, ymax);
    bool inside2 = is_inside(x2, y2, xmin, xmax, ymin, ymax);

    // Case 1: First vertex inside
    if (inside1) {
      x_clipped.push_back(x1);
      y_clipped.push_back(y1);
    }

    // Case 2: Edge intersects the bounding box
    double x_int, y_int;
    if (find_intersection(x1, y1, x2, y2, xmin, xmax, ymin, ymax, x_int, y_int)) {
      x_clipped.push_back(x_int);
      y_clipped.push_back(y_int);
    }
  }

  // Output as a matrix
  int n_clipped = x_clipped.size();
  NumericMatrix clipped_coords(n_clipped, 2);
  for (int k = 0; k < n_clipped; ++k) {
    clipped_coords(k, 0) = x_clipped[k];
    clipped_coords(k, 1) = y_clipped[k];
  }

  return List::create(_["clipped_polygon"] = clipped_coords);
}
