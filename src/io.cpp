#include "io.h"
#include "geom.h"

#include <sstream>
#include <iostream>

/** Parse a file and return a vertex list
 */
VertexList
load_vertices(std::istream& in) {
  std::string line;

  VertexList res;

  int idx = 0;
  while (std::getline(in, line)) {
    if (line.rfind("#", 0) == 0) continue; // comment
    std::stringstream ss(line);
    double x, y;
    ss >> x >> y;
    res.emplace_back(Vertex(x, y, idx++));
  }

  return res;
}

std::pair<VertexList, std::vector<std::pair<unsigned,unsigned>>>
load_obj(std::istream& in) {
  std::string line;

  VertexList vertex_list;
  std::vector<std::pair<unsigned,unsigned>> edge_list;

  int idx = 0;
  while (std::getline(in, line)) {
    if (line.rfind("#", 0) == 0) continue; // comment
    if (line.rfind("v", 0) == 0) {
      std::stringstream ss(line);
      char type;
      double x, y;
      ss >> type >> x >> y;
      vertex_list.emplace_back(Vertex(x, y, idx++));
    } else if ((line.rfind("f", 0) == 0) ||
                line.rfind("l", 0) == 0) {
      std::stringstream ss(line);
      char type;
      unsigned first_vertex;
      ss >> type;
      ss >> first_vertex;
      unsigned prev_vertex = first_vertex;
      unsigned v;
      while (!ss.fail()) {
        ss >> v;
        edge_list.emplace_back(std::make_pair(prev_vertex, v));
      }
      if (type == 'f') {
        edge_list.emplace_back(std::make_pair(v, first_vertex));
      }
    }
  }

  return std::make_pair(std::move(vertex_list), std::move(edge_list));
}
