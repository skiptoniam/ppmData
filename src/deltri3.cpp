// Adapted from the delaunator-cpp code from https://github.com/delfrrr/delaunator-cpp

#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "deltri3.h"
// using namespace Rcpp;

namespace deltri {

 size_t fast_mod(const size_t i, const size_t c) {
    return i >= c ? i % c : i;
}

 double sum(const std::vector<double>& x) {
    double sum = x[0];
    double err = 0.0;

    for (size_t i = 1; i < x.size(); i++) {
        const double k = x[i];
        const double m = sum + k;
        err += std::fabs(sum) >= std::fabs(k) ? sum - m + k : k - m + sum;
        sum = m;
    }
    return sum + err;
}

 double dist(
    const double ax,
    const double ay,
    const double bx,
    const double by) {
    const double dx = ax - bx;
    const double dy = ay - by;
    return dx * dx + dy * dy;
}

 double circumradius(const site& p1, const site& p2, const site& p3){
    site d = site::vector(p1, p2);
    site e = site::vector(p1, p3);

    const double bl = d.magnitude2();
    const double cl = e.magnitude2();
    const double det = site::determinant(d, e);

    site radius((e.y() * bl - d.y() * cl) * 0.5 / det,
                 (d.x() * cl - e.x() * bl) * 0.5 / det);

    if ((bl > 0.0 || bl < 0.0) &&
        (cl > 0.0 || cl < 0.0) &&
        (det > 0.0 || det < 0.0))
        return radius.magnitude2();
    return (std::numeric_limits<double>::max)();
}

 double circumradius(
    const double ax,
    const double ay,
    const double bx,
    const double by,
    const double cx,
    const double cy) {
    const double dx = bx - ax;
    const double dy = by - ay;
    const double ex = cx - ax;
    const double ey = cy - ay;

    const double bl = dx * dx + dy * dy;
    const double cl = ex * ex + ey * ey;
    const double d = dx * ey - dy * ex;

    const double x = (ey * bl - dy * cl) * 0.5 / d;
    const double y = (dx * cl - ex * bl) * 0.5 / d;

    if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) && (d > 0.0 || d < 0.0)) {
        return x * x + y * y;
    } else {
        return (std::numeric_limits<double>::max)();
    }
}

 bool clockwise(const site& p0, const site& p1, const site& p2){
    site v0 = site::vector(p0, p1);
    site v1 = site::vector(p0, p2);
    double det = site::determinant(v0, v1);
    double dist = v0.magnitude2() + v1.magnitude2();
    double dist2 = site::dist2(v0, v1);
    if (det == 0)
    {
        return false;
    }
    double reldet = std::abs(dist / det);
    if (reldet > 1e14)
        return false;
    return det < 0;
}

 bool clockwise(double px, double py, double qx, double qy,
    double rx, double ry){
    site p0(px, py);
    site p1(qx, qy);
    site p2(rx, ry);
    return clockwise(p0, p1, p2);
}

 bool counterclockwise(const site& p0, const site& p1, const site& p2){
    site v0 = site::vector(p0, p1);
    site v1 = site::vector(p0, p2);
    double det = site::determinant(v0, v1);
    double dist = v0.magnitude2() + v1.magnitude2();
    double dist2 = site::dist2(v0, v1);
    if (det == 0)
        return false;
    double reldet = std::abs(dist / det);
    if (reldet > 1e14)
        return false;
    return det > 0;
}

 bool counterclockwise(double px, double py, double qx, double qy,
    double rx, double ry){
    site p0(px, py);
    site p1(qx, qy);
    site p2(rx, ry);
    return counterclockwise(p0, p1, p2);
}


 site circumcenter(
    const double ax,
    const double ay,
    const double bx,
    const double by,
    const double cx,
    const double cy) {
    const double dx = bx - ax;
    const double dy = by - ay;
    const double ex = cx - ax;
    const double ey = cy - ay;

    const double bl = dx * dx + dy * dy;
    const double cl = ex * ex + ey * ey;
    const double d = dx * ey - dy * ex;

    const double x = ax + (ey * bl - dy * cl) * 0.5 / d;
    const double y = ay + (dx * cl - ex * bl) * 0.5 / d;

    return site(x, y);
}

 bool in_circle(
    const double ax,
    const double ay,
    const double bx,
    const double by,
    const double cx,
    const double cy,
    const double px,
    const double py) {
    const double dx = ax - px;
    const double dy = ay - py;
    const double ex = bx - px;
    const double ey = by - py;
    const double fx = cx - px;
    const double fy = cy - py;

    const double ap = dx * dx + dy * dy;
    const double bp = ex * ex + ey * ey;
    const double cp = fx * fx + fy * fy;

    return (dx * (ey * cp - bp * fy) -
            dy * (ex * cp - bp * fx) +
            ap * (ex * fy - ey * fx)) < 0.0;
}

constexpr double EPSILON = std::numeric_limits<double>::epsilon();

 bool check_pts_equal(double x1, double y1, double x2, double y2) {
    return std::fabs(x1 - x2) <= EPSILON &&
           std::fabs(y1 - y2) <= EPSILON;
}

double pseudo_angle(const double dx, const double dy) {
    const double p = dx / (std::abs(dx) + std::abs(dy));
    return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0; // [0..1)
}


deltri_cpp::deltri_cpp(std::vector<double> const& in_coords)
    : coords(in_coords), m_points(in_coords){
    std::size_t n = coords.size() >> 1;

    std::vector<std::size_t> ids(n);
    std::iota(ids.begin(), ids.end(), 0);

    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double min_x = (std::numeric_limits<double>::max)();
    double min_y = (std::numeric_limits<double>::max)();
    for (const site& p : m_points)
    {
        min_x = std::min(p.x(), min_x);
        min_y = std::min(p.y(), min_y);
        max_x = std::max(p.x(), max_x);
        max_y = std::max(p.y(), max_y);
    }
    double width = max_x - min_x;
    double height = max_y - min_y;
    double span = width * width + height * height; // Everything is square dist.

    site center((min_x + max_x) / 2, (min_y + max_y) / 2);

    std::size_t i0 = INVALID_INDEX;
    std::size_t i1 = INVALID_INDEX;
    std::size_t i2 = INVALID_INDEX;

    double min_dist = (std::numeric_limits<double>::max)();
    for (size_t i = 0; i < m_points.size(); ++i)
    {
        const site& p = m_points[i];
        const double d = site::dist2(center, p);
        if (d < min_dist) {
            i0 = i;
            min_dist = d;
        }
    }

    const site& p0 = m_points[i0];

    min_dist = (std::numeric_limits<double>::max)();

    for (std::size_t i = 0; i < n; i++) {
        if (i == i0) continue;
        const double d = site::dist2(p0, m_points[i]);
        if (d < min_dist && d > 0.0) {
            i1 = i;
            min_dist = d;
        }
    }

    const site& p1 = m_points[i1];

    double min_radius = (std::numeric_limits<double>::max)();

    for (std::size_t i = 0; i < n; i++) {
        if (i == i0 || i == i1) continue;

        const double r = circumradius(p0, p1, m_points[i]);
        if (r < min_radius) {
            i2 = i;
            min_radius = r;
        }
    }

    if (!(min_radius < (std::numeric_limits<double>::max)())) {
        throw std::runtime_error("not triangulation");
    }

    const site& p2 = m_points[i2];

    if (counterclockwise(p0, p1, p2))
        std::swap(i1, i2);

    double i0x = p0.x();
    double i0y = p0.y();
    double i1x = m_points[i1].x();
    double i1y = m_points[i1].y();
    double i2x = m_points[i2].x();
    double i2y = m_points[i2].y();

    m_center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);

    std::vector<double> dists;
    dists.reserve(m_points.size());
    for (const site& p : m_points)
        dists.push_back(dist(p.x(), p.y(), m_center.x(), m_center.y()));

     std::sort(ids.begin(), ids.end(),
        [&dists](std::size_t i, std::size_t j)
            { return dists[i] < dists[j]; });

    m_hash_size = static_cast<std::size_t>(std::ceil(std::sqrt(n)));
    m_hash.resize(m_hash_size);
    std::fill(m_hash.begin(), m_hash.end(), INVALID_INDEX);

    hull_prev.resize(n);
    hull_next.resize(n);
    hull_tri.resize(n);

    hull_start = i0;

    size_t hull_size = 3;

    hull_next[i0] = hull_prev[i2] = i1;
    hull_next[i1] = hull_prev[i0] = i2;
    hull_next[i2] = hull_prev[i1] = i0;

    hull_tri[i0] = 0;
    hull_tri[i1] = 1;
    hull_tri[i2] = 2;

    m_hash[hash_key(i0x, i0y)] = i0;
    m_hash[hash_key(i1x, i1y)] = i1;
    m_hash[hash_key(i2x, i2y)] = i2;

    std::size_t max_triangles = n < 3 ? 1 : 2 * n - 5;
    triangles.reserve(max_triangles * 3);
    halfedges.reserve(max_triangles * 3);
    add_triangle(i0, i1, i2, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
    double xp = std::numeric_limits<double>::quiet_NaN();
    double yp = std::numeric_limits<double>::quiet_NaN();

     for (std::size_t k = 0; k < n; k++) {
        const std::size_t i = ids[k];
        const double x = coords[2 * i];
        const double y = coords[2 * i + 1];

        if (k > 0 && check_pts_equal(x, y, xp, yp))
            continue;
        xp = x;
        yp = y;

        if (check_pts_equal(x, y, i0x, i0y) ||
            check_pts_equal(x, y, i1x, i1y) ||
            check_pts_equal(x, y, i2x, i2y)) continue;

        std::size_t start = 0;

        size_t key = hash_key(x, y);
        for (size_t j = 0; j < m_hash_size; j++) {
            start = m_hash[fast_mod(key + j, m_hash_size)];

             if (start != INVALID_INDEX && start != hull_next[start])
                break;
        }

        assert(hull_prev[start] != start);
        assert(hull_prev[start] != INVALID_INDEX);

        start = hull_prev[start];
        size_t e = start;
        size_t q;

        while (true){
            q = hull_next[e];
            if (site::equal(m_points[i], m_points[e], span) ||
                site::equal(m_points[i], m_points[q], span))
            {
                e = INVALID_INDEX;
                break;
            }
            if (counterclockwise(x, y, coords[2 * e], coords[2 * e + 1],
                coords[2 * q], coords[2 * q + 1]))
                break;
            e = q;
            if (e == start) {
                e = INVALID_INDEX;
                break;
            }
        }

        if (e == INVALID_INDEX)
            continue;

        std::size_t t = add_triangle(
            e,
            i,
            hull_next[e],
            INVALID_INDEX,
            INVALID_INDEX,
            hull_tri[e]);

        hull_tri[i] = legalize(t + 2);
        hull_tri[e] = t;
        hull_size++;

        std::size_t next = hull_next[e];
        while (true){
            q = hull_next[next];
            if (!counterclockwise(x, y, coords[2 * next], coords[2 * next + 1],
                coords[2 * q], coords[2 * q + 1]))
                break;
            t = add_triangle(next, i, q,
                hull_tri[i], INVALID_INDEX, hull_tri[next]);
            hull_tri[i] = legalize(t + 2);
            hull_next[next] = next; // mark as removed
            hull_size--;
            next = q;
        }

        if (e == start) {
            while (true)
            {
                q = hull_prev[e];
                if (!counterclockwise(x, y, coords[2 * q], coords[2 * q + 1],
                    coords[2 * e], coords[2 * e + 1]))
                    break;
                t = add_triangle(q, i, e,
                    INVALID_INDEX, hull_tri[e], hull_tri[q]);
                legalize(t + 2);
                hull_tri[q] = t;
                hull_next[e] = e;
                hull_size--;
                e = q;
            }
        }

        hull_prev[i] = e;
        hull_start = e;
        hull_prev[next] = i;
        hull_next[e] = i;
        hull_next[i] = next;

        m_hash[hash_key(x, y)] = i;
        m_hash[hash_key(coords[2 * e], coords[2 * e + 1])] = e;
    }

    hull.resize(hull_size);
    size_t e = hull_start;
    for (int i = 0; i < hull_size; i++) {
        hull[i] = e;
        e = hull_next[e];
    }

}

std::vector<double> deltri_cpp::get_hull_sites(){
    std::vector<double> hull_sites;
    size_t e = hull_start;
    size_t cnt = 1;
    do {
        hull_sites.push_back(coords[2 * e]);
        hull_sites.push_back(coords[2 * e + 1]);
        cnt++;
        e = hull_next[e];
    } while (e != hull_start);
    return hull_sites;
}

double deltri_cpp::get_hull_area(){
    std::vector<double> hull_area;
    size_t e = hull_start;
    size_t cnt = 1;
    do {
        hull_area.push_back((coords[2 * e] - coords[2 * hull_prev[e]]) *
            (coords[2 * e + 1] + coords[2 * hull_prev[e] + 1]));
        cnt++;
        e = hull_next[e];
    } while (e != hull_start);
    return sum(hull_area);
}

double deltri_cpp::get_triangle_area(){
    std::vector<double> vals;
    for (size_t i = 0; i < triangles.size(); i += 3)
    {
        const double ax = coords[2 * triangles[i]];
        const double ay = coords[2 * triangles[i] + 1];
        const double bx = coords[2 * triangles[i + 1]];
        const double by = coords[2 * triangles[i + 1] + 1];
        const double cx = coords[2 * triangles[i + 2]];
        const double cy = coords[2 * triangles[i + 2] + 1];
        double val = std::fabs((by - ay) * (cx - bx) - (bx - ax) * (cy - by));
        vals.push_back(val);
    }
    return sum(vals);
}

std::size_t deltri_cpp::legalize(std::size_t a) {
    std::size_t i = 0;
    std::size_t ar = 0;
    m_edge_stack.clear();

    while (true) {
        const size_t b = halfedges[a];

        const size_t a0 = 3 * (a / 3);
        ar = a0 + (a + 2) % 3;

        if (b == INVALID_INDEX) {
            if (i > 0) {
                i--;
                a = m_edge_stack[i];
                continue;
            } else {
                //i = INVALID_INDEX;
                break;
            }
        }

        const size_t b0 = 3 * (b / 3);
        const size_t al = a0 + (a + 1) % 3;
        const size_t bl = b0 + (b + 2) % 3;

        const std::size_t p0 = triangles[ar];
        const std::size_t pr = triangles[a];
        const std::size_t pl = triangles[al];
        const std::size_t p1 = triangles[bl];

        const bool illegal = in_circle(
            coords[2 * p0],
            coords[2 * p0 + 1],
            coords[2 * pr],
            coords[2 * pr + 1],
            coords[2 * pl],
            coords[2 * pl + 1],
            coords[2 * p1],
            coords[2 * p1 + 1]);

        if (illegal) {
            triangles[a] = p1;
            triangles[b] = p0;

            auto hbl = halfedges[bl];

            if (hbl == INVALID_INDEX) {
                std::size_t e = hull_start;
                do {
                    if (hull_tri[e] == bl) {
                        hull_tri[e] = a;
                        break;
                    }
                    e = hull_prev[e];
                } while (e != hull_start);
            }
            link(a, hbl);
            link(b, halfedges[ar]);
            link(ar, bl);
            std::size_t br = b0 + (b + 1) % 3;

            if (i < m_edge_stack.size()) {
                m_edge_stack[i] = br;
            } else {
                m_edge_stack.push_back(br);
            }
            i++;

        } else {
            if (i > 0) {
                i--;
                a = m_edge_stack[i];
                continue;
            } else {
                break;
            }
        }
    }
    return ar;
}

std::size_t deltri_cpp::hash_key(const double x, const double y) const {
    const double dx = x - m_center.x();
    const double dy = y - m_center.y();
    return fast_mod(
        static_cast<std::size_t>(std::llround(std::floor(pseudo_angle(dx, dy) * static_cast<double>(m_hash_size)))),
        m_hash_size);
}

std::size_t deltri_cpp::add_triangle(
    std::size_t i0,
    std::size_t i1,
    std::size_t i2,
    std::size_t a,
    std::size_t b,
    std::size_t c) {
    std::size_t t = triangles.size();
    triangles.push_back(i0);
    triangles.push_back(i1);
    triangles.push_back(i2);
    link(t, a);
    link(t + 1, b);
    link(t + 2, c);
    return t;
}

void deltri_cpp::link(const std::size_t a, const std::size_t b) {
    std::size_t s = halfedges.size();
    if (a == s) {
        halfedges.push_back(b);
    } else if (a < s) {
        halfedges[a] = b;
    } else {
        throw std::runtime_error("Cannot link edge");
    }
    if (b != INVALID_INDEX) {
        std::size_t s2 = halfedges.size();
        if (b == s2) {
            halfedges.push_back(a);
        } else if (b < s2) {
            halfedges[b] = a;
        } else {
            throw std::runtime_error("Cannot link edge");
        }
    }
}

}

//[[Rcpp::export]]
Rcpp::List deltri_cpp(std::vector<double> coords) {

  deltri::deltri_cpp del(coords);

  return Rcpp::List::create(Rcpp::_["coords"] = del.coords,
                            Rcpp::_["triangles"] = del.triangles,
                            Rcpp::_["halfedges"] = del.halfedges,
                            Rcpp::_["convexhull"] = del.hull,
                            Rcpp::_["convexhull.area"] = del.get_hull_area());

}
