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

  /** determine Vertices' a, b, c's relative orientation.
   *
   * returns 1 if they are counterclockwise, 0 if collinear, -1 if clockwise.
   */
  static int orientation(const Vertex& a, const Vertex &b, const Vertex &c) {
    double det1 = (a.x - c.x) * (b.y - c.y);
    double det2 = (a.y - c.y) * (b.x - c.x);
    double det = det1 - det2;
    return signum(det);
  }
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

public:
  void print_tip_and_neighbors() const {
    assert(opposite);

    const Vertex &tip = *v;
    const Vertex &left = *next_constrained->v;
    const Vertex &right = *opposite->prev_constrained->prev_constrained->v;

    DBG(DBG_GENERIC) << "right (" << right.x << ", " << right.y << ")";
    DBG(DBG_GENERIC) << "tip   (" << tip.x << ", " << tip.y << "); tail (" << get_tail()->x << ", " << get_tail()->y << ")";
    DBG(DBG_GENERIC) << "left  (" << left.x << ", " << left.y << ")";
    DBG(DBG_GENERIC) << "det (r,t,l)" << Vertex::orientation(right, tip, left);
  }

protected:
  bool check_unconstrain_at_tip() const {
    assert(opposite);

    const Vertex &tip = *v;
    const Vertex &left = *next_constrained->v;
    const Vertex &right = *opposite->prev_constrained->prev_constrained->v;

    return ((&left != &right) && Vertex::orientation(right, tip, left) >= 0);
  }
public:
  /** Check whether this edge can be removed/unconstrained
   *
   * Edges can be removed if they are not on the CH,
   * and if removing them will not result in a reflex vertex
   * at their tip or tail.
   */
  bool can_unconstrain() const {
    assert(is_constrained);
    return (opposite != NULL &&
      check_unconstrain_at_tip() &&
      opposite->check_unconstrain_at_tip());
  }

  bool get_is_constrained() const { return is_constrained; };
  Edge* get_opposite() const { return opposite; };
  const Vertex* get_head() const { return v; };
  const Vertex* get_tail() const { return prev->v; };


  /** Mark this edge as not a constrained edge.
   *
   * Update pointers to next constrained in neighbors.
   */
  void unconstrain() {
    assert(can_unconstrain());
    assert_valid();

    next_constrained->prev_constrained = opposite->prev_constrained;
    prev_constrained->next_constrained = opposite->next_constrained;

    opposite->next_constrained->prev_constrained = prev_constrained;
    opposite->prev_constrained->next_constrained = next_constrained;

    is_constrained = false;
    prev_constrained = NULL;
    next_constrained = NULL;
    opposite->is_constrained = false;
    opposite->prev_constrained = NULL;
    opposite->next_constrained = NULL;

    assert_valid();
  };

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

  const Vertex * const vertex_base_ptr;
  int num_faces = 0;
public:
  DECL(const VertexList& vertices);

  void unconstrain_all();

#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {};
#endif

  void write_obj_segments(const VertexList * vertices, std::ostream &o) const;
  int get_num_faces() const { return num_faces; };

  friend std::ostream& operator<<(std::ostream&, const DECL&);
};

std::ostream& operator<<(std::ostream&, const DECL&);
