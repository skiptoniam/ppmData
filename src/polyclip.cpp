#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "polyclip.h"
#include "deltri3.h"

using namespace Rcpp;

// ============================================================
// polyclip namespace: geometry primitives
// ============================================================

namespace polyclip {

// 2D cross product of vectors (b - a) and (c - a)
static inline double cross2d(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// Line-line intersection: segment (p1,p2) with line through (a,b)
static Point line_intersection(const Point& p1, const Point& p2,
                               const Point& a, const Point& b) {
    double A1 = p2.y - p1.y;
    double B1 = p1.x - p2.x;
    double C1 = A1 * p1.x + B1 * p1.y;

    double A2 = b.y - a.y;
    double B2 = a.x - b.x;
    double C2 = A2 * a.x + B2 * a.y;

    double det = A1 * B2 - A2 * B1;
    if (std::abs(det) < 1e-15) {
        // Parallel lines — return midpoint as fallback
        return Point((p1.x + p2.x) * 0.5, (p1.y + p2.y) * 0.5);
    }
    return Point((C1 * B2 - C2 * B1) / det,
                 (A1 * C2 - A2 * C1) / det);
}

// Is point on the inside (left side) of directed edge a -> b?
static inline bool is_inside_edge(const Point& p, const Point& a, const Point& b) {
    return cross2d(a, b, p) >= 0.0;
}

// --- public API ---

double signed_area(const std::vector<Point>& poly) {
    double area = 0.0;
    size_t n = poly.size();
    if (n < 3) return 0.0;
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        area += (poly[i].x * poly[j].y) - (poly[j].x * poly[i].y);
    }
    return area * 0.5;
}

double polygon_area(const std::vector<Point>& poly) {
    return std::abs(signed_area(poly));
}

int winding_number(const Point& p, const std::vector<Point>& poly) {
    int wn = 0;
    size_t n = poly.size();
    for (size_t i = 0; i < n; i++) {
        const Point& v1 = poly[i];
        const Point& v2 = poly[(i + 1) % n];
        if (v1.y <= p.y) {
            if (v2.y > p.y) {
                if (cross2d(v1, v2, p) > 0.0)
                    wn++;
            }
        } else {
            if (v2.y <= p.y) {
                if (cross2d(v1, v2, p) < 0.0)
                    wn--;
            }
        }
    }
    return wn;
}

std::vector<Point> sutherland_hodgman(const std::vector<Point>& subject,
                                       const std::vector<Point>& clip) {
    if (subject.empty() || clip.empty()) return {};

    std::vector<Point> output = subject;
    size_t clip_size = clip.size();

    for (size_t i = 0; i < clip_size; i++) {
        if (output.empty()) return {};

        std::vector<Point> input;
        input.swap(output);

        const Point& a = clip[i];
        const Point& b = clip[(i + 1) % clip_size];

        size_t input_size = input.size();
        for (size_t j = 0; j < input_size; j++) {
            const Point& p = input[j];
            const Point& q = input[(j + 1) % input_size];

            bool p_in = is_inside_edge(p, a, b);
            bool q_in = is_inside_edge(q, a, b);

            if (q_in) {
                if (!p_in) {
                    output.push_back(line_intersection(p, q, a, b));
                }
                output.push_back(q);
            } else if (p_in) {
                output.push_back(line_intersection(p, q, a, b));
            }
        }
    }
    return output;
}

void ensure_ccw(std::vector<Point>& poly) {
    if (signed_area(poly) < 0.0) {
        std::reverse(poly.begin(), poly.end());
    }
}

} // namespace polyclip


// ============================================================
// Voronoi cell extraction (reuses deltri triangulation)
// ============================================================

static std::vector<polyclip::Point> get_voronoi_cell(
        size_t i,
        const std::vector<size_t>& inedges,
        const deltri::deltri_cpp& del) {

    std::vector<polyclip::Point> cell;
    size_t e0 = inedges[i];
    if (e0 == deltri::INVALID_INDEX) return cell;

    size_t e = e0;
    do {
        size_t tri = e / 3;

        // triangle vertex indices
        size_t v0 = del.triangles[3 * tri];
        size_t v1 = del.triangles[3 * tri + 1];
        size_t v2 = del.triangles[3 * tri + 2];

        // circumcenter
        double ax = del.coords[v0 * 2],     ay = del.coords[v0 * 2 + 1];
        double bx = del.coords[v1 * 2],     by = del.coords[v1 * 2 + 1];
        double cx = del.coords[v2 * 2],     cy = del.coords[v2 * 2 + 1];

        double dx = bx - ax, dy = by - ay;
        double ex = cx - ax, ey = cy - ay;
        double bl = dx * dx + dy * dy;
        double cl = ex * ex + ey * ey;
        double d  = dx * ey - dy * ex;

        double ccx = ax + (ey * bl - dy * cl) * 0.5 / d;
        double ccy = ay + (dx * cl - ex * bl) * 0.5 / d;

        cell.push_back(polyclip::Point(ccx, ccy));

        // advance to next halfedge around the point
        size_t e_next = (e % 3 == 2) ? e - 2 : e + 1;
        e = del.halfedges[e_next];
    } while (e != e0 && e != deltri::INVALID_INDEX);

    return cell;
}

// Build the inedges lookup (one incoming halfedge per point)
static std::vector<size_t> build_inedges(const deltri::deltri_cpp& del) {
    size_t np = del.coords.size() / 2;
    std::vector<size_t> inedges(np, deltri::INVALID_INDEX);

    for (size_t e = 0; e < del.halfedges.size(); e++) {
        size_t p = del.triangles[e % 3 == 2 ? e - 2 : e + 1];
        if (del.halfedges[e] == deltri::INVALID_INDEX ||
            inedges[p] == deltri::INVALID_INDEX) {
            inedges[p] = e;
        }
    }
    return inedges;
}

// A single polygon part: outer ring + holes
struct ClipPart {
    std::vector<polyclip::Point> outer_ring;
    std::vector<std::vector<polyclip::Point>> holes;
};

// Build clip polygon vectors from R inputs (multiple parts)
static std::vector<ClipPart> build_clip_data_multi(
        Rcpp::List parts_outer_x, Rcpp::List parts_outer_y,
        Rcpp::List parts_hole_x_list, Rcpp::List parts_hole_y_list) {

    int nparts = parts_outer_x.size();
    std::vector<ClipPart> parts(nparts);

    for (int p = 0; p < nparts; p++) {
        NumericVector clip_x = parts_outer_x[p];
        NumericVector clip_y = parts_outer_y[p];

        int n_outer = clip_x.size();
        parts[p].outer_ring.resize(n_outer);
        for (int i = 0; i < n_outer; i++) {
            parts[p].outer_ring[i] = polyclip::Point(clip_x[i], clip_y[i]);
        }
        polyclip::ensure_ccw(parts[p].outer_ring);

        Rcpp::List hole_x_list = parts_hole_x_list[p];
        Rcpp::List hole_y_list = parts_hole_y_list[p];
        int nholes = hole_x_list.size();
        parts[p].holes.resize(nholes);
        for (int h = 0; h < nholes; h++) {
            NumericVector hx = hole_x_list[h];
            NumericVector hy = hole_y_list[h];
            int nh = hx.size();
            parts[p].holes[h].resize(nh);
            for (int i = 0; i < nh; i++) {
                parts[p].holes[h][i] = polyclip::Point(hx[i], hy[i]);
            }
            polyclip::ensure_ccw(parts[p].holes[h]);
        }
    }

    return parts;
}

// Check if all cell vertices are inside a single polygon part
static bool cell_fully_inside_part(
        const std::vector<polyclip::Point>& cell,
        const ClipPart& part) {

    for (const auto& pt : cell) {
        if (polyclip::winding_number(pt, part.outer_ring) == 0)
            return false;
    }
    for (size_t h = 0; h < part.holes.size(); h++) {
        for (const auto& pt : cell) {
            if (polyclip::winding_number(pt, part.holes[h]) != 0)
                return false;
        }
    }
    return true;
}

// Check if all cell vertices are inside ANY polygon part
static bool cell_fully_inside(
        const std::vector<polyclip::Point>& cell,
        const std::vector<ClipPart>& parts) {

    for (const auto& part : parts) {
        if (cell_fully_inside_part(cell, part))
            return true;
    }
    return false;
}

// Compute clipped area for one Voronoi cell against a single part
static double clip_cell_area_part(
        const std::vector<polyclip::Point>& cell,
        const ClipPart& part) {

    std::vector<polyclip::Point> clipped =
        polyclip::sutherland_hodgman(part.outer_ring, cell);
    double area = polyclip::polygon_area(clipped);

    for (size_t h = 0; h < part.holes.size(); h++) {
        std::vector<polyclip::Point> hole_clipped =
            polyclip::sutherland_hodgman(part.holes[h], cell);
        area -= polyclip::polygon_area(hole_clipped);
    }

    return std::max(0.0, area);
}

// Compute clipped area summed over all polygon parts
static double clip_cell_area(
        const std::vector<polyclip::Point>& cell,
        const std::vector<ClipPart>& parts) {

    double total = 0.0;
    for (const auto& part : parts) {
        total += clip_cell_area_part(cell, part);
    }
    return total;
}


// ============================================================
// Rcpp exports
// ============================================================

// [[Rcpp::export]]
NumericVector dirtess_clip_areas_cpp(
        std::vector<double> coords,
        int ncoords,
        Rcpp::List parts_outer_x, Rcpp::List parts_outer_y,
        Rcpp::List parts_hole_x_list, Rcpp::List parts_hole_y_list) {

    // 1. Delaunay triangulation
    deltri::deltri_cpp del(coords);

    // 2. Inedges lookup
    std::vector<size_t> inedges = build_inedges(del);

    // 3. Build clip geometry (multiple parts)
    std::vector<ClipPart> parts = build_clip_data_multi(
        parts_outer_x, parts_outer_y,
        parts_hole_x_list, parts_hole_y_list);

    // 4. Compute clipped area for each real point
    NumericVector areas(ncoords);

    for (int i = 0; i < ncoords; i++) {
        std::vector<polyclip::Point> cell = get_voronoi_cell(i, inedges, del);

        if (cell.empty()) {
            areas[i] = 0.0;
            continue;
        }
        polyclip::ensure_ccw(cell);

        // Fast path: cell fully inside study region
        if (cell_fully_inside(cell, parts)) {
            areas[i] = polyclip::polygon_area(cell);
            continue;
        }

        areas[i] = clip_cell_area(cell, parts);
    }

    return areas;
}


// [[Rcpp::export]]
Rcpp::List dirtess_clip_cpp(
        std::vector<double> coords,
        int ncoords,
        Rcpp::List parts_outer_x, Rcpp::List parts_outer_y,
        Rcpp::List parts_hole_x_list, Rcpp::List parts_hole_y_list) {

    // 1. Delaunay triangulation
    deltri::deltri_cpp del(coords);

    // 2. Inedges lookup
    std::vector<size_t> inedges = build_inedges(del);

    // 3. Build clip geometry (multiple parts)
    std::vector<ClipPart> parts = build_clip_data_multi(
        parts_outer_x, parts_outer_y,
        parts_hole_x_list, parts_hole_y_list);

    int nparts = parts.size();

    // 4. For each real point: clip against all parts and collect results
    NumericVector areas(ncoords);
    Rcpp::List poly_list(ncoords);

    for (int i = 0; i < ncoords; i++) {
        std::vector<polyclip::Point> cell = get_voronoi_cell(i, inedges, del);

        if (cell.empty()) {
            areas[i] = 0.0;
            NumericMatrix empty_mat(0, 2);
            poly_list[i] = Rcpp::List::create(
                Rcpp::_["outer"] = empty_mat,
                Rcpp::_["holes"] = Rcpp::List(0));
            continue;
        }
        polyclip::ensure_ccw(cell);

        // Collect clipped polygon fragments from all parts
        double total_area = 0.0;
        std::vector<Rcpp::NumericMatrix> all_outers;
        std::vector<Rcpp::NumericMatrix> all_holes;

        for (int p = 0; p < nparts; p++) {
            bool fi = cell_fully_inside_part(cell, parts[p]);

            // --- outer ring for this part ---
            std::vector<polyclip::Point> clipped_outer;
            if (fi) {
                clipped_outer = cell;
            } else {
                clipped_outer = polyclip::sutherland_hodgman(
                    parts[p].outer_ring, cell);
            }

            if (clipped_outer.empty()) continue;

            double area_outer = polyclip::polygon_area(clipped_outer);

            // Convert outer to closed Nx2 matrix
            int no = clipped_outer.size();
            NumericMatrix outer_mat(no + 1, 2);
            for (int k = 0; k < no; k++) {
                outer_mat(k, 0) = clipped_outer[k].x;
                outer_mat(k, 1) = clipped_outer[k].y;
            }
            outer_mat(no, 0) = clipped_outer[0].x;
            outer_mat(no, 1) = clipped_outer[0].y;
            all_outers.push_back(outer_mat);

            // --- holes for this part ---
            double area_holes = 0.0;
            if (!fi) {
                int nholes = parts[p].holes.size();
                for (int h = 0; h < nholes; h++) {
                    std::vector<polyclip::Point> hclip =
                        polyclip::sutherland_hodgman(parts[p].holes[h], cell);
                    if (hclip.empty()) continue;

                    double ha = polyclip::polygon_area(hclip);
                    area_holes += ha;

                    int nh = hclip.size();
                    NumericMatrix hmat(nh + 1, 2);
                    for (int k = 0; k < nh; k++) {
                        hmat(k, 0) = hclip[k].x;
                        hmat(k, 1) = hclip[k].y;
                    }
                    hmat(nh, 0) = hclip[0].x;
                    hmat(nh, 1) = hclip[0].y;
                    all_holes.push_back(hmat);
                }
            }

            total_area += std::max(0.0, area_outer - area_holes);
        }

        areas[i] = total_area;

        // Use the first (largest) outer fragment as the primary polygon
        // and collect all holes
        if (all_outers.empty()) {
            NumericMatrix empty_mat(0, 2);
            poly_list[i] = Rcpp::List::create(
                Rcpp::_["outer"] = empty_mat,
                Rcpp::_["holes"] = Rcpp::List(0));
        } else {
            poly_list[i] = Rcpp::List::create(
                Rcpp::_["outer"] = all_outers[0],
                Rcpp::_["holes"] = Rcpp::wrap(all_holes));
        }
    }

    return Rcpp::List::create(
        Rcpp::_["areas"] = areas,
        Rcpp::_["polygons"] = poly_list);
}
