#pragma once

#include "geom.h"

#include <istream>

std::unique_ptr<VertexList>
load_vertices(std::istream& in);
