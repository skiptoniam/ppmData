#ifndef deltri_hh
#define deltri_hh

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

    std::vector<std::size_t> triangles;

    std::vector<std::size_t> halfedges;

    std::vector<std::size_t> hull_prev;
    std::vector<std::size_t> hull_next;

    std::vector<std::size_t> hull;

    std::vector<std::size_t> hull_tri;
    std::size_t hull_start;

    deltri_cpp(std::vector<double> const& in_coords);
    std::vector<double> get_hull_sites();
    double get_hull_area();
    double get_triangle_area();

private:
    std::vector<std::size_t> m_hash;
    site m_center;
    std::size_t m_hash_size;
    std::vector<std::size_t> m_edge_stack;

    std::size_t legalize(std::size_t a);
    std::size_t hash_key(double x, double y) const;
    std::size_t add_triangle(
        std::size_t i0,
        std::size_t i1,
        std::size_t i2,
        std::size_t a,
        std::size_t b,
        std::size_t c);
    void link(std::size_t a, std::size_t b);
};

}

#endif

