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
  bool vertex_on_outer_removal_boundary = false;

public:
#ifndef NDEBUG
  const int idx_;
  Vertex(double x_, double y_, int idx)
    : x(x_)
    , y(y_)
    , idx_(idx)
  {}
  int idx() const { return idx_; };
#else
  Vertex(double x_, double y_, int)
    : x(x_)
    , y(y_)
  {}
  int idx() const { return 0; };
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
  const int idx_;
#endif
public:
  Edge(Edge *next_, Edge *prev_, Edge *opposite_, Vertex *v_, [[maybe_unused]] int idx)
    : opposite(opposite_)
    , next(next_)
    , prev(prev_)
    , next_constrained(next_)
    , prev_constrained(prev_)
    , v(v_)
#ifndef NDEBUG
    , idx_(idx)
#endif
  {}
#ifndef NDEBUG
  int idx() const { return idx_; };
#else
  int idx() const { return 0; };
#endif

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
   * Upd/It
   * ate pointers to next constrained in neighbors.
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
  /** Iterate around v once, returning a handle for each face.
   *
   * In general, we'll walk around a vertex in clockwise order.
   *
   * However, if we hit the Convex Hull we restart at the initial vertex
   * and continue from there counter-clockwise.
   */
  class AroundVertexFacesIterator {
  protected:
    bool hit_ch_ = false;

    Edge * const e_start;
    Edge * e_cur;

  public:
    AroundVertexFacesIterator(Edge* e_vertex)
      : e_start(e_vertex)
      , e_cur(e_vertex)
    {
      assert(e_vertex);
      assert(e_vertex->is_constrained);
    }

    Edge* operator*() const { return e_cur; }
    Edge* operator->() const { return e_cur; }

    AroundVertexFacesIterator& operator++() {
      assert(e_cur != NULL);
      if (!hit_ch_) {
        assert(e_cur->v == e_start->v);
        e_cur = e_cur->next_constrained->opposite;
        if (e_cur == e_start) {
          e_cur = NULL;
          return *this;
        } else if (e_cur == NULL) {
          DBG(DBG_ITERATORS) << "Smacked into the CH.";
          hit_ch_ = true;
          e_cur = e_start;
        }
      }
      if (hit_ch_) {
        if (!e_cur->opposite) {
          e_cur = NULL;
          return *this;
        }
        e_cur = e_cur->opposite->prev_constrained;
        assert(e_cur->v == e_start->v);
      }
      return *this;
    }

    bool hit_ch() const { return hit_ch_; }
  };

  /** Iterate around v once, returning a handle for each face.
   *
   * This iterator returns all faces in strictly clockwise order, even
   * when running into the CH.
   */
  class AroundVertexFacesStriclyClockwiseIterator : public AroundVertexFacesIterator {
  private:
    using Base = AroundVertexFacesIterator;
    Edge *edge_before_ch_ = NULL;
    Edge *edge_after_ch_ = NULL;
  public:
    AroundVertexFacesStriclyClockwiseIterator(Edge* e_vertex)
      : Base(e_vertex)
    {}

    AroundVertexFacesStriclyClockwiseIterator& operator++() {
      Edge *edge_before = e_cur;
      Base::operator++();

      if (! hit_ch_) {
        return *this;
      } else {
        assert(edge_before_ch_ == NULL);
        assert(edge_after_ch_ == NULL);
        hit_ch_ = false;

        edge_before_ch_ = edge_before;
        if (e_cur == NULL) {
          // The face/edge we started at was already on the CH.
          edge_after_ch_ = e_start;
        } else {
          // Search the extremal edge on the counter-clockwise side of the start edge.
          while (e_cur->opposite) {
            e_cur = e_cur->opposite->prev_constrained;
            assert(e_cur->v == e_start->v);
          }
          edge_after_ch_ = e_cur;
        }
        return *this;
      }
    }
    bool hit_ch() const { return (edge_after_ch_ != NULL); }
    Edge * edge_after_ch() const { return edge_after_ch_; }
    Edge * edge_before_ch() const { return edge_before_ch_; }
  };

  /** State of a region of the DECL.
   *
   * We store things in threes, so elements for triangle 0 are at pos 0, 1, 2,
   * and for triangle n at pos 3n, 3n+1, 3n+2.
   *
   * edge 0 is opposite vertex 0.
   *
   * This is just like triangle handles things.
   */
  class SavedState {
  public:
    std::vector<Vertex*> triangle_vertices; /** For each triangle, the three vertices making up its boundary */
    std::vector<int> buddy;                 /** For each triangle edge, the index of its buddy edge. */
    std::vector<int> constraints;           /** For each triangle edge, if it is constrained, the index of the next and prev constrained edge.  -1 otherwise. */
    int num_faces;
    int num_vertices;
    int num_vertices_on_boundary;
    SavedState()
      : num_faces(0)
      , num_vertices(0)
      , num_vertices_on_boundary(0)
      { }
    SavedState& operator= (const SavedState&) = default;
    SavedState& operator= (SavedState&&) = default;
    SavedState(SavedState&& o) = default;
    SavedState(const FixedVector<Edge> &parentlist, const std::vector<Edge*>& edge_pointers, int num_faces_, int num_vertices_, int num_vertices_on_boundary_);
    SavedState(const FixedVector<Edge> &alledges, int num_faces_, int num_vertices_, int num_vertices_on_boundary_);
  };
  #if 0
  /** State of a region cut out from a DECL.
   *
   * the boundary vertices are in order, witnessed by edges on that boundary,
   * except if we skipped over removed vertices that were on the CH.
   *
   * We can tell whether a vertex (witnessed by e) is such one by checking
   * whether the previous vertex in the list is e's tail.
   */
  class SavedStateFromParent : public SavedState {
  public:
    std::vector<Edge*> boundary_vertices;   /** List of boundary vertices, in order. */
    std::vector<Edge*> edge_pointers;       /** Pointers into the parent's edge list of where things come from */
    SavedStateFromParent(const FixedVector<Edge> &parentlist, const std::vector<Edge*>& edge_pointers,
  };
  #endif


private:
  Edge *get_next_face_cw_around_vertex(Edge *e) const;
  Edge *get_next_face_ccw(Edge *e) const;

private:
  static void decl_triangulate_prepare(const VertexList& vertices, struct triangulateio& tin);
  void decl_triangulate_process(VertexList& vertices, const struct triangulateio& tout);
  void decl_triangulate(VertexList& vertices);

  std::shared_ptr<VertexList> all_vertices; /* potentially more than this decl handles */

  FixedVector<Edge> edges;
  int num_vertices;
  int num_vertices_on_boundary;
  int num_triangles;
  int num_faces;

private:
  std::vector<Edge*> vertices_to_remove;
  std::vector<Edge*> removal_boundary_vertices;
  std::vector<Edge*> halfedges_to_remove;
  int vertices_to_remove_on_ch = 0;
  int faces_removed = 0;

  void shoot_hole_select_vertices(unsigned size);
  void shoot_hole_identify_affected_elements_around_vertex(Edge* const e_vertex);
  void shoot_hole_list_triangles_in_face(Edge *e);
  Edge* shoot_hole_identify_boundary_vertices_start();
  bool shoot_hole_identify_boundary_vertices();
  bool shoot_hole_identify_affected_elements();
  void shoot_hole_inject_saved_state(const SavedState &state);
  void shoot_hole(unsigned size, int num_iterations, int max_recurse);
  void shoot_holes(int max_recurse);

#ifndef NDEBUG
  void assert_hole_shooting_vertices_clean() const {
    #if 1
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_to_be_removed == false;} ) );
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_on_removal_boundary == false;} ) );
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_on_outer_removal_boundary == false;} ) );
    #else
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_to_be_removed == false
                                                                          && e.v->vertex_on_removal_boundary == false
                                                                          && e.v->vertex_on_outer_removal_boundary == false;}) );
    #endif
  }
  void assert_hole_shooting_reset() const {
    #if 1
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.triangle_to_be_removed == false;}) );
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_to_be_removed == false;} ) );
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_on_removal_boundary == false;} ) );
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_on_outer_removal_boundary == false;} ) );
    #else
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.triangle_to_be_removed == false
                                                                          && e.v->vertex_to_be_removed == false
                                                                          && e.v->vertex_on_removal_boundary == false
                                                                          && e.v->vertex_on_outer_removal_boundary == false;}) );
    #endif
    assert(vertices_to_remove.size() == 0);
    assert(removal_boundary_vertices.size() == 0);
    assert(halfedges_to_remove.size() == 0);
    assert(vertices_to_remove_on_ch == 0);
    assert(faces_removed == 0);
  }
  bool vertex_is_on_ch(Edge* e) const { /* expensive */
    auto it = AroundVertexFacesIterator(e);
    for (; *it; ++it) { }
    return it.hit_ch();
  }
  void assert_vertex_is_on_ch(Edge* e) const {
    assert(vertex_is_on_ch(e));
  }
#else
  void assert_hole_shooting_vertices_clean() const {}
  void assert_hole_shooting_reset() const {}
  void assert_vertex_is_on_ch(Edge *) const {}
#endif


public:
  DECL(std::shared_ptr<VertexList> all_vertices);
  DECL(std::shared_ptr<VertexList> all_vertices, const SavedState &state);

  void find_convex_decomposition(int max_recurse=1);
  SavedState find_convex_decomposition_many(int num_iterations, int num_faces_to_beat, int max_recurse);
  void unconstrain_all();
  void reset_constraints();

#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {}
#endif

  void write_obj_segments(bool dump_vertices, std::ostream &o) const;
  int get_num_faces() const { return num_faces; }

  friend std::ostream& operator<<(std::ostream&, const DECL&);
};

std::ostream& operator<<(std::ostream&, const Vertex&);
std::ostream& operator<<(std::ostream&, const Edge&);
std::ostream& operator<<(std::ostream&, const DECL&);
