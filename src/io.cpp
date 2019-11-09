#include "io.h"
#include "geom.h"

#include <sstream>
#include <iostream>

/** Parse a file and return a vertex list
 */
std::shared_ptr<VertexList>
load_vertices(std::istream& in) {
  std::string line;

  auto res = std::make_shared<VertexList>();

  int idx = 0;
  while (std::getline(in, line)) {
    if (line.rfind("#", 0) == 0) continue; // comment
    std::stringstream ss(line);
    double x, y;
    ss >> x >> y;
    res->emplace_back(Vertex(x, y, idx++));
  }

  return res;
}
