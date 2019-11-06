#pragma once

#include "mcd.h"
#include "tools.h"

#include <vector>

class Vertex {
  public:
  const double x;
  const double y;
  const int idx;
  public:
  Vertex(double x_, double y_, int idx_)
    : x(x_)
    , y(y_)
    , idx(idx_)
  {}
};
using VertexList = std::vector<Vertex>;

class Node : public Vertex {
  public:
  Node(double x_, double y_, int idx_)
    : Vertex(x_, y_, idx_)
  {}
};

class Edge;
class DECL;

/** A multi-layer half-edge data structure.
 *
 * Edges come in two forms:  constraints and triangulation edges.
 *
 * The triangulation edges themselves make up a classic DECL,
 * with buddy-pointer, next and prev pointer.
 *
 * Some of the triangulation edges are also constraints.  These will,
 * additionally have next/prev links to the next constraint-edge.
 */
class Edge {
  //friend class DECL;
  friend std::ostream& operator<<(std::ostream&, const DECL&);

  bool is_constrained;
  bool on_ch;
  Edge *opposite, *prev, *next;
  Edge *prev_constrained, *next_constrained;
  Node *n;
public:
  Edge()
    : is_constrained(0)
    , on_ch(0)
    , opposite(0)
    , prev(0)
    , next(0)
    , prev_constrained(0)
    , next_constrained(0)
  {}
  Edge(Edge *prev_, Edge *next_, Edge *buddy)
    : is_constrained(1)
    , on_ch(buddy == NULL)
    , opposite(buddy)
    , prev(prev_)
    , next(next_)
    , prev_constrained(prev_)
    , next_constrained(next_)
  {}
};

class DECL {
  FixedVector<Edge> edges;

  static void decl_triangulate_prepare(const VertexList& vertices, struct triangulateio& tin);
  void decl_triangulate_process(const VertexList& vertices, const struct triangulateio& tout);
  void decl_triangulate(const VertexList& vertices);

  public:
  /** Initialize the DECL with the vertices and a triangulation of their CH */
  DECL(const VertexList& vertices);

  friend std::ostream& operator<<(std::ostream&, const DECL&);
};

std::ostream& operator<<(std::ostream&, const DECL&);
