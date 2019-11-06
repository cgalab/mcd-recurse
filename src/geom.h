#pragma once

#include "mcd.h"
#include "tools.h"

#include <vector>

class Vertex {
  public:
  const double x;
  const double y;
#ifndef NDEBUG
  const int idx;
  Vertex(double x_, double y_, int idx_)
    : x(x_)
    , y(y_)
    , idx(idx_)
  {}
#else
  Vertex(double x_, double y_, int)
    : x(x_)
    , y(y_)
  {}
#endif
};
using VertexList = std::vector<Vertex>;

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
 * additionally have next/prev links to the next constrained edge.
 */
class Edge {
  //friend class DECL;
  friend std::ostream& operator<<(std::ostream&, const DECL&);

  bool is_constrained;    /** Whether this is a constrained edge, i.e., one that is a boundary of a face in our convex decomposition */
  Edge *opposite;         /** Pointer to the buddy of this edge.  NULL indicates that this edge is on the CH. */
  Edge *next;             /** Pointer to the next edge of this triangle. This edge will start at Vertex v. */
  Edge *prev;             /** Pointer to the previous face of this triangle. */
  Edge *next_constrained; /** If constrained, pointer to the next constrained edge of this (decomposition) face.  This edge will start at Vertex v. */
  Edge *prev_constrained; /** If constrained, pointer to the previous constrained edge of this (decomposition) face. */
  const Vertex *v;        /** Vertex this edge points to.  Tail of prev and prev_constrained. */
public:
  Edge()
    : is_constrained(0)
    , opposite(0)
    , next(0)
    , prev(0)
    , next_constrained(0)
    , prev_constrained(0)
    , v(NULL)
  {}
  Edge(Edge *next_, Edge *prev_, Edge *opposite_, const Vertex *v_)
    : is_constrained(1)
    , opposite(opposite_)
    , next(next_)
    , prev(prev_)
    , next_constrained(next_)
    , prev_constrained(prev_)
    , v(v_)
  {}

protected:
  bool check_unconstrain_at_tip() const {
    const Vertex *tip = v;
    const Vertex *left = next_constrained->v;
    const Vertex *right = opposite->prev_constrained->prev_constrained->v;
    return false;
  }
public:
  /** Check whether this edge can be removed/unconstrained
   *
   * Edges can be removed if they are not on the CH,
   * and if removing them will not result in a reflex vertex
   * at their tip or tail.
   */
  bool can_unconstrain() const {
    if (opposite == NULL) return false;
  }

#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {};
#endif
};

class DECL {
  FixedVector<Edge> edges;

  static void decl_triangulate_prepare(const VertexList& vertices, struct triangulateio& tin);
  void decl_triangulate_process(const VertexList& vertices, const struct triangulateio& tout);
  void decl_triangulate(const VertexList& vertices);

public:
  /** Initialize the DECL with the vertices and a triangulation of their CH */
  DECL(const VertexList& vertices);

#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {};
#endif
  friend std::ostream& operator<<(std::ostream&, const DECL&);
};

std::ostream& operator<<(std::ostream&, const DECL&);
