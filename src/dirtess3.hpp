#pragma once

#include "deltri.hpp"

namespace dirtess {

  class polygon
  {
  private:
    std::vector<deltri::site> points;

  public:
    polygon(){};
    polygon(const deltri::site & one);
    polygon(const std::vector<deltri::site> & points);

    double x(int i) const{ return points[i].x(); }
    double y(int i) const{ return points[i].y(); }

    void push(deltri::site newpoints) { // note: not a pointer
      this->points.push_back(newpoints);  // note: not address-of
    }

    size_t size() const { return points.size(); }

    void resize(int sz) {
      points.resize(sz);
    }

  };

}
