#pragma once

#include "geom.h"

#include <istream>

VertexList load_vertices(std::istream& in);
std::pair<VertexList, InputEdgeSet> load_obj(std::istream& in);
