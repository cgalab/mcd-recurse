#pragma once

#include "mcd.h"
#include "tools.h"

#include <vector>

class Edge;
class DECL;

class Vertex {

public:
  const double x;
  const double y;

private:
  /* When shooting holes in the DECL, we store information here */
  bool vertex_to_be_removed = false;
  bool vertex_on_removal_boundary = false;

public:
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

  friend class DECL;
  friend std::ostream& operator<<(std::ostream&, const Vertex&);
};
using VertexList = std::vector<Vertex>;

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
  friend class DECL;
  friend std::ostream& operator<<(std::ostream&, const DECL&);

  bool is_constrained = 1; /** Whether this is a constrained edge, i.e., one that is a boundary of a face in our convex decomposition */
  Edge *opposite;          /** Pointer to the buddy of this edge.  NULL indicates that this edge is on the CH. */
  Edge *next;              /** Pointer to the next edge of this triangle. This edge will start at Vertex v. */
  Edge *prev;              /** Pointer to the previous face of this triangle. */
  Edge *next_constrained;  /** If constrained, pointer to the next constrained edge of this (decomposition) face.  This edge will start at Vertex v. */
  Edge *prev_constrained;  /** If constrained, pointer to the previous constrained edge of this (decomposition) face. */
  Vertex *v;               /** Vertex this edge points to.  Tail of prev and prev_constrained. */
  bool triangle_to_be_removed = false; /** During hole shooting, we use this as a flag. */
#ifndef NDEBUG
  const int idx;
#endif
public:
  Edge(Edge *next_, Edge *prev_, Edge *opposite_, Vertex *v_, [[maybe_unused]] int idx_)
    : opposite(opposite_)
    , next(next_)
    , prev(prev_)
    , next_constrained(next_)
    , prev_constrained(prev_)
    , v(v_)
#ifndef NDEBUG
    , idx(idx_)
#endif
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

  Vertex* get_tail() const { return prev->v; }

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
  }

  /** This is run during DECL reset which will mark all edges as constrained.
   *
   * Mark this edge constrained, and set next_constrained/prev_constrained to
   * next/prev as they soon will be constrained even if they aren't yet.
   */
  void reset_all_constraints() {
    is_constrained = true;
    prev_constrained = prev;
    next_constrained = next;
  }

#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {}
#endif

  friend std::ostream& operator<<(std::ostream&, const Edge&);
};

class DECL {
  FixedVector<Edge> edges;

  static void decl_triangulate_prepare(const VertexList& vertices, struct triangulateio& tin);
  void decl_triangulate_process(VertexList& vertices, const struct triangulateio& tout);
  void decl_triangulate(VertexList& vertices);

  const Vertex * const vertex_base_ptr;
  int num_vertices = 0;
  int num_triangles = 0;
  int num_faces = 0;

private:
  std::vector<Edge*> vertices_to_remove;
  std::vector<Edge*> removal_boundary_vertices;
  std::vector<Edge*> halfedges_to_remove;
  int vertices_to_remove_on_ch = 0;
  int faces_removed = 0;

  void shoot_hole_select_vertices(unsigned size);
  void shoot_hole_identify_affected_elements_around_vertex(Edge* const e_vertex);
  void shoot_hole_list_triangles_in_face(Edge *e);
  void shoot_hole_identify_enclosed_boundary_vertices();
  void shoot_hole_identify_affected_elements();

public:
  DECL(VertexList& vertices);

  void find_convex_decomposition() {
    unconstrain_all();
    shoot_hole(2);
    //shoot_hole(sqrt(num_vertices));
  }
  void unconstrain_all();
  void reset_constraints();

  void shoot_hole(unsigned size);

#ifndef NDEBUG
  void assert_valid() const;
  void assert_hole_shooting_reset() const {
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.triangle_to_be_removed == false
                                                                          && e.v->vertex_to_be_removed == false
                                                                          && e.v->vertex_on_removal_boundary == false;}) );
    assert(vertices_to_remove.size() == 0);
    assert(removal_boundary_vertices.size() == 0);
    assert(halfedges_to_remove.size() == 0);
    assert(vertices_to_remove_on_ch == 0);
    assert(faces_removed == 0);
  }
#else
  void assert_valid() const {}
  void assert_hole_shooting_reset() const {}
#endif

  void write_obj_segments(const VertexList * vertices, std::ostream &o) const;
  int get_num_faces() const { return num_faces; }

  friend std::ostream& operator<<(std::ostream&, const DECL&);
};

std::ostream& operator<<(std::ostream&, const Vertex&);
std::ostream& operator<<(std::ostream&, const Edge&);
std::ostream& operator<<(std::ostream&, const DECL&);
