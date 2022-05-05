#pragma once

#ifdef DELTRI_HEADER_ONLY
#define INLINE inline
#else
#define INLINE
#endif

#include <limits>
#include <vector>
#include <ostream>

namespace deltri {

constexpr std::size_t INVALID_INDEX =
    (std::numeric_limits<std::size_t>::max)();

class site
{
public:
    site(double x, double y) : m_x(x), m_y(y)
    {}
    site() : m_x(0), m_y(0)
    {}

    double x() const
    { return m_x; }

    double y() const
    { return m_y; }

    double magnitude2() const
    { return m_x * m_x + m_y * m_y; }

    static double determinant(const site& p1, const site& p2)
    {
        return p1.m_x * p2.m_y - p1.m_y * p2.m_x;
    }

    static site vector(const site& p1, const site& p2)
    {
        return site(p2.m_x - p1.m_x, p2.m_y - p1.m_y);
    }

    static double dist2(const site& p1, const site& p2)
    {
        site vec = vector(p1, p2);
        return vec.m_x * vec.m_x + vec.m_y * vec.m_y;
    }

    static bool equal(const site& p1, const site& p2, double span)
    {
        double dist = dist2(p1, p2) / span;

        // ABELL - This number should be examined to figure how how
        // it correlates with the breakdown of calculating determinants.
        return dist < 1e-20;
    }

private:
    double m_x;
    double m_y;
};

inline std::ostream& operator<<(std::ostream& out, const site& p)
{
    out << p.x() << "/" << p.y();
    return out;
}


class sites
{
public:
    using const_iterator = site const *;

    sites(const std::vector<double>& coords) : m_coords(coords)
    {}

    const site& operator[](size_t offset)
    {
        return reinterpret_cast<const site&>(
            *(m_coords.data() + (offset * 2)));
    };

    sites::const_iterator begin() const
        { return reinterpret_cast<const site *>(m_coords.data()); }
    sites::const_iterator end() const
        { return reinterpret_cast<const site *>(
            m_coords.data() + m_coords.size()); }
    size_t size() const
        { return m_coords.size() / 2; }

private:
    const std::vector<double>& m_coords;
};

class deltri_cpp {

public:
    std::vector<double> const& coords;
    sites m_points;

    // 'triangles' stores the indices to the 'X's of the input
    // 'coords'.
    std::vector<std::size_t> triangles;

    // 'halfedges' store indices into 'triangles'.  If halfedges[X] = Y,
    // It says that there's an edge from X to Y where a) X and Y are
    // both indices into triangles and b) X and Y are indices into different
    // triangles in the array.  This allows you to get from a triangle to
    // its adjacent triangle.  If the a triangle edge has no adjacent triangle,
    // its half edge will be INVALID_INDEX.
    std::vector<std::size_t> halfedges;

    std::vector<std::size_t> hull_prev;
    std::vector<std::size_t> hull_next;

    // vector for hull sites
    std::vector<std::size_t> hull;

    // This contains indexes into the triangles array.
    std::vector<std::size_t> hull_tri;
    std::size_t hull_start;

    INLINE deltri_cpp(std::vector<double> const& in_coords);
    INLINE std::vector<double> get_hull_sites();
    INLINE double get_hull_area();
    INLINE double get_triangle_area();

private:
    std::vector<std::size_t> m_hash;
    site m_center;
    std::size_t m_hash_size;
    std::vector<std::size_t> m_edge_stack;

    INLINE std::size_t legalize(std::size_t a);
    INLINE std::size_t hash_key(double x, double y) const;
    INLINE std::size_t add_triangle(
        std::size_t i0,
        std::size_t i1,
        std::size_t i2,
        std::size_t a,
        std::size_t b,
        std::size_t c);
    INLINE void link(std::size_t a, std::size_t b);
};

} //namespace deltri

#undef INLINE