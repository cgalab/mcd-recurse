#pragma once

#include "geom.h"

#include <istream>

VertexList load_vertices(std::istream& in);
std::pair<VertexList, std::vector<std::pair<unsigned,unsigned>>> load_obj(std::istream& in);
