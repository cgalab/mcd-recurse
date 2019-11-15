#pragma once

#include "mcd.h"
#include "tools.h"

#include <vector>

class Edge;
class DECL;

// {{{ Vertex
class Vertex {
public:
  const double x;
  const double y;

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
// }}}

// {{{ Edge
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
class BaseEdge {
  protected:
  bool is_constrained = 1; /** Whether this is a constrained edge, i.e., one that is a boundary of a face in our convex decomposition */
  Edge *opposite;          /** Pointer to the buddy of this edge.  NULL indicates that this edge is on the CH. */
  Edge *next;              /** Pointer to the next edge of this triangle. This edge will start at Vertex v. */
  Edge *prev;              /** Pointer to the previous face of this triangle. */
  Edge *next_constrained;  /** If constrained, pointer to the next constrained edge of this (decomposition) face.  This edge will start at Vertex v. */
  Edge *prev_constrained;  /** If constrained, pointer to the previous constrained edge of this (decomposition) face. */
  Vertex *v;               /** Vertex this edge points to.  Tail of prev and prev_constrained. */
  unsigned working_set_depth = 0; /** The number of the working set this edge belongs to. */
  bool triangle_marked = false; /** During hole shooting, we use this as a flag of whether this triangle's face belongs to the next working set. */

  BaseEdge() {};
  BaseEdge(Edge *next_, Edge *prev_, Edge *opposite_, Vertex *v_)
    : opposite(opposite_)
    , next(next_)
    , prev(prev_)
    , next_constrained(next_)
    , prev_constrained(prev_)
    , v(v_)
  {}
};
class Edge : public BaseEdge {
  using Base = BaseEdge;

  friend class DECL;
  friend std::ostream& operator<<(std::ostream&, const DECL&);

#ifndef NDEBUG
  const int idx_;
#endif
public:
  Edge(Edge *next_, Edge *prev_, Edge *opposite_, Vertex *v_, [[maybe_unused]] int idx)
    : Base(next_, prev_, opposite_, v_)
#ifndef NDEBUG
    , idx_(idx)
#endif
  {}
  /* Only support setting state from Edges with the same idx. */
  Edge& operator=(const Edge& o) noexcept {
    assert(idx_ == o.idx_);
    static_cast<Base&>(*this) = o;
    return *this;
  }

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
    assert(opposite == NULL || working_set_depth == opposite->working_set_depth);
    return (opposite != NULL &&
      check_unconstrain_at_tip() &&
      opposite->check_unconstrain_at_tip());
  }

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

  /** Check whether this edge can be flipped.
   *
   * It assumes the graph is triangulated and everything is still constrained. */
  bool can_flip() const {
    if (opposite == NULL) return false;

    assert(working_set_depth == opposite->working_set_depth);

    const Vertex &tip = *v;
    const Vertex &tail = *opposite->v;
    const Vertex &left = *next->v;
    const Vertex &right = *opposite->next->v;

    assert(&left != &right);
    assert(Vertex::orientation(tail, tip, left) > 0);
    assert(Vertex::orientation(tip, tail, right) > 0);

    return ((Vertex::orientation(right, left, tail) > 0) &&
            (Vertex::orientation(left, right, tip) > 0));
  }

  /** Flip this triangulation edge.
   *
   * Note that this only works on triangulation edges and completely
   * ignores next/prev constrained pointers.  Those are assumed to
   * be fixed with reset_all_constraints() soon.
   *
   * As such, assert_valid() will fail until reset_all_constraints() is run.
   */
  void flip() {
    assert(can_flip());
    const auto triangle = [](Edge* a, Edge *b, Edge *c) {
      const auto link_edges = [](Edge* first, Edge *second) {
        first->next = second;
        second->prev = first;
      };
      link_edges(a, b);
      link_edges(b, c);
      link_edges(c, a);
    };

    Edge * l = this;     /* This edge points upwards */
    Edge * r = opposite; /* Our buddy is to our right, pointing downwards */
    Edge * tl = l->next; /* top left */
    Edge * bl = l->prev; /* bottom left */
    Edge * tr = r->prev; /* top right */
    Edge * br = r->next; /* bottom right */

    /* After the flip, this edge points from right to left and is "below" our buddy which points from left to right */
    l->v = tl->v;
    r->v = br->v;
    triangle(l, bl, br);
    triangle(r, tr, tl);
  }

  Vertex* get_tail() const { return prev->v; }

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

  /** Returns if this vertex is of degree 4 or higher.
   *
   * Maybe also, at some point, if it's a collinear point and of degree >=3,
   * and also if this is a vertex on the convex hull with degree > 2;
   *
   * Basically. this is to know if we should bother trying to improve things here.
   *
   * Requires edge to be a constrained edge.
   */
  bool vertex_is_of_higher_degree() const {
    assert(is_constrained);
    if (opposite &&
        next_constrained->opposite &&
        next_constrained->opposite->next_constrained == opposite) {
      /* Degree 2 vertex */
      return false;
    }

    if (opposite &&
        next_constrained->opposite &&
        next_constrained->opposite->next_constrained->opposite &&
        next_constrained->opposite->next_constrained->opposite->next_constrained == opposite) {
      /* Degree 3 vertex */
      /* Could still have an incident angle of exactly pi; check for that */
      const Vertex &v0 = *v;
      const Vertex &v1 = *opposite->v;
      const Vertex &v2 = *next_constrained->opposite->next_constrained->v;
      const Vertex &v3 = *next_constrained->v;
      if (Vertex::orientation(v0, v1, v3) == 0) {
        return true;
      };
      if (Vertex::orientation(v0, v2, v3) == 0) {
        return true;
      };
      if (Vertex::orientation(v0, v3, v1) == 0) {
        return true;
      };

      return false;
    }

    if (!opposite &&
        !next_constrained->opposite) {
      /* Degree 2 vertex on the convex hull */
      return false;
    }

    return true;
  }

#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {}
#endif

  friend std::ostream& operator<<(std::ostream&, const Edge&);
};
// }}}

class DECL {
  /** Iterate around v once, returning a handle for each face.
   *
   * In general, we'll walk around a vertex in clockwise order.
   *
   * However, if we hit the Convex Hull we restart at the initial vertex
   * and continue from there counter-clockwise.
   */
  // {{{ Iterators
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
  // }}}}}}

  class WorkingSet {
  public:
    unsigned depth = 0; /** A recursion depth marker to match with each edge's working_set_depth */
    std::vector<Edge*> my_edges; /** The list of edges we are working on right now, including the boundary edges */
    std::vector<Edge*> shuffled_edges; /** The list of edges we are working on right now, excluding the boundary edges, in random order */
    unsigned num_my_triangles; /** The number of triangles incident to my_edges.
                                   This is identical to the number of faces when all
                                   of my_edges are constrained.  Note that not all
                                   these triangles need to be surrounded by my_edges;
                                   other (constrained) edges can border them too. */
    unsigned num_faces_mine_constrained; /** The number of faces, in the entire DECL, when all of my triangles are their own face.  */
    WorkingSet() = default;

    WorkingSet(const WorkingSet&) = delete;
    WorkingSet(WorkingSet&&) noexcept = default;
    WorkingSet& operator=(WorkingSet&&) noexcept = default;

    WorkingSet(unsigned depth_, std::vector<Edge*>&& my_edges_, std::vector<Edge*>&& shuffled_edges_, int num_my_triangles_, int num_faces_mine_constrained_)
      : depth(depth_)
      , my_edges(std::move(my_edges_))
      , shuffled_edges(std::move(shuffled_edges_))
      , num_my_triangles(num_my_triangles_)
      , num_faces_mine_constrained(num_faces_mine_constrained_)
    {}
  };

  class SavedDecomposition {
    public:
    std::vector<Edge> edge_content; /** What was in the edges */
    unsigned saved_num_faces; /** The number of faces of the entire graph with this edge_conent. */

    SavedDecomposition(const SavedDecomposition*) = delete;
    SavedDecomposition(const WorkingSet& ws, unsigned num_faces);
  };


  /* Helper functions */
  private:
    Edge *get_next_face_cw_around_vertex(Edge *e) const;
    Edge *get_next_face_ccw(Edge *e) const;

  /* Helper functions */
  private:
    /* setup */
    static void decl_triangulate_prepare(const VertexList& vertices, struct triangulateio& tin);
    static std::pair<FixedVector<Edge>, std::vector<Edge*>> decl_triangulate_process(VertexList& vertices, const struct triangulateio& tout);
    static std::pair<FixedVector<Edge>, std::vector<Edge*>> decl_triangulate(VertexList& vertices);

    /* Decomposition */
    void flip_random_edges_and_reset_constraints();
    void unconstrain_random_edges();
    void unconstrain_random_edges_initial_improvement();
    void reinject_saved_decomposition(SavedDecomposition&& saved_state);

  /* The state of the DECL */
  private:
    /* These stay fixed over all iterations.
     *
     * (not necessarily their content, but * at least the set and order) */
    VertexList all_vertices; /** The list of all vertices. */

    FixedVector<Edge> all_edges; /** The list of all edges */

    WorkingSet working_set;
    unsigned num_faces; /** The number of faces of the entire graph right now. */

  /* The state of DECL hole finding */
  private:
    std::vector<Edge*> marked_halfedges;
    std::vector<Edge*> marking_candidates;
    unsigned num_marked_faces = 0;
    unsigned num_marked_triangles = 0;

  /* DECL hole finding functions and state mgmt */
  private:
    unsigned shoot_hole_mark_triangles_in_face(Edge * const e);
    bool shoot_hole_select_triangles(unsigned num_triangles);
    std::vector<Edge*> shoot_hole_interior_edges() const;
    bool shoot_hole(unsigned size);
    void shoot_holes();

  private:
    /* private constructor to make the public one use decl_triangulate's result. */
    DECL(VertexList&& vertices, std::pair<FixedVector<Edge>, std::vector<Edge*>>&& triangulation_result);

    /* optimization in punched hole */
    void find_convex_decomposition_many(unsigned num_iterations);

    std::geometric_distribution<unsigned> geometric_distribution;
  /* public interface */
  public:
    const unsigned hole_size_base = 7;
    const double hole_size_geometric_param = 0.4;
    const double flip_nums_exponent = 1./5;

    /** Initialize the DECL with the vertices and a triangulation of their CH */
    DECL(VertexList&& vertices)
    : DECL(std::forward<VertexList>(vertices), decl_triangulate(vertices)) {}

    void reset_constraints();
    void find_convex_decomposition();

    void write_obj_segments(bool dump_vertices, std::ostream &o) const;
    unsigned get_num_faces() const { return num_faces; }

    friend std::ostream& operator<<(std::ostream&, const DECL&);

  // Debugging things:
public:
#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {}
#endif

private:
#ifndef NDEBUG
  void assert_hole_shooting_reset() const {
    assert( std::all_of(all_edges.begin(), all_edges.end(), [](const Edge& e){return e.triangle_marked == false;}) );
    // assert( std::all_of(working_set.my_edges.begin(), working_set.my_edges.end(), [](const Edge* e){return e->triangle_marked == false;}) );

    assert(marked_halfedges.size() == 0);
    assert(marking_candidates.size() == 0);
    assert(num_marked_faces == 0);
    assert(num_marked_triangles == 0);
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
  void assert_hole_shooting_reset() const {}
  void assert_vertex_is_on_ch(Edge *) const {}
#endif


};

std::ostream& operator<<(std::ostream&, const Vertex&);
std::ostream& operator<<(std::ostream&, const Edge&);
std::ostream& operator<<(std::ostream&, const DECL&);

/* vim: set fdm=marker: */
