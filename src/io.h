#pragma once

#include "geom.h"

#include <istream>

std::shared_ptr<VertexList>
load_vertices(std::istream& in);
