/* CG:SHOP 2020: Minimum Convex Decomposition -- Recursor Tool
*
*  Copyright 2019, 2020 Peter Palfraader
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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

std::pair<VertexList, InputEdgeSet>
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
        edge_list.emplace_back(sorted_pair(prev_vertex-1, v-1));
      }
      if (type == 'f') {
        edge_list.emplace_back(sorted_pair(v-1, first_vertex-1));
      }
    }
  }
  /*
  sort(edge_list.begin(), edge_list.end());
  edge_list.erase(unique(edge_list.begin(), edge_list.end()), edge_list.end());
  return std::make_pair(std::move(vertex_list), std::move(edge_list));
  */

  InputEdgeSet edge_set(edge_list.size());
  for (auto& i: edge_list)
      edge_set.insert(i);

  return std::make_pair(std::move(vertex_list), std::move(edge_set));
}
