#pragma once

#include "mcd.h"

#include <vector>

class Vertex {
    const double x;
    const double y;
  public:
    Vertex(double x_, double y_)
      : x(x_)
      , y(y_)
    {}
};

using VertexList = std::vector<Vertex>;
