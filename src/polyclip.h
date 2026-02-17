#ifndef POLYCLIP_H
#define POLYCLIP_H

#include <vector>
#include <cmath>

namespace polyclip {

struct Point {
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x_, double y_) : x(x_), y(y_) {}
};

// Signed area via shoelace formula (positive for CCW, negative for CW)
double signed_area(const std::vector<Point>& poly);

// Absolute area of polygon
double polygon_area(const std::vector<Point>& poly);

// Winding number point-in-polygon test (non-zero = inside)
int winding_number(const Point& p, const std::vector<Point>& poly);

// Sutherland-Hodgman: clip subject polygon against convex clip polygon
// Both polygons should be CCW. The clip polygon MUST be convex.
std::vector<Point> sutherland_hodgman(const std::vector<Point>& subject,
                                       const std::vector<Point>& clip);

// Ensure polygon vertices are in counter-clockwise order
void ensure_ccw(std::vector<Point>& poly);

} // namespace polyclip

#endif
