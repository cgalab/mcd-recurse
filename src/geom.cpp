#include "geom.h"
#include "triangle.h"

#include <numeric>
#include <cmath> 

#define NUMBER_OF_HOLE_PUNCHES_SCALE 10

static const int next_edge_offset[] = {1,2,0};
static const int prev_edge_offset[] = {2,0,1}; /* also same as vertex index that an edge points towards */

/** Find our index (0,1,2) in a neighboring triangle as returned by triangle.
 */
static inline int
tnidx_by_nidx(const int *t, int needle) {
  int r = t[0] == needle ? 0 :
      t[1] == needle ? 1 :
      t[2] == needle ? 2 : -1;
  assert(r >= 0);
  return r;
}


#ifndef NDEBUG
/** Local sanity check. */
void
Edge::
assert_valid() const {
  assert(!opposite || this == opposite->opposite);
  assert(this == prev->next);
  assert(this == next->prev);
  assert(!opposite || v == opposite->prev->v);
  if (is_constrained) {
    assert(this == prev_constrained->next_constrained);
    assert(this == next_constrained->prev_constrained);
    assert(!opposite || v == opposite->prev_constrained->v);
  } else {
    assert(prev_constrained == NULL);
    assert(next_constrained == NULL);
  }
}
#endif

DECL::
DECL(VertexList&& vertices, std::pair<FixedVector<Edge>, unsigned>&& triangulation_result)
  : all_vertices(std::forward<VertexList>(vertices))
  , number_of_ch_vertices(triangulation_result.second)
  , all_edges(std::forward<FixedVector<Edge>>(triangulation_result.first))
{
  num_faces = 1 + (all_edges.size() + number_of_ch_vertices) / 2 - all_vertices.size();

  working_set.my_edges.resize(all_edges.size());
  std::iota(std::begin(working_set.my_edges), std::end(working_set.my_edges), all_edges.data());
  working_set.num_my_triangles = num_faces;
  working_set.num_faces_mine_constrained = num_faces;
}

/** Prepare triangle's in/out data structure with the vertex list
 */
void
DECL::
decl_triangulate_prepare(const VertexList& vertices, struct triangulateio& tin) {
  int num_v = vertices.size();

  /* Load vertex coordinates */
  tin.numberofpoints = num_v;
  tin.pointlist = (double *) my_malloc_c(num_v*2 * sizeof(double));
  double *dp = tin.pointlist;
  for (int i=0; i<num_v; ++i) {
    *(dp++) = vertices[i].x;
    *(dp++) = vertices[i].y;
  }
  assert(dp == tin.pointlist + tin.numberofpoints*2);
}

/** Process triangle's in/out data structure and create the DECL
 */
std::pair<FixedVector<Edge>, unsigned>
DECL::
decl_triangulate_process(VertexList& vertices, const struct triangulateio& tout) {
  FixedVector<Edge> edges;

  unsigned num_t = tout.numberoftriangles;
  unsigned num_e = tout.numberofedges;
  unsigned num_ch_v = tout.numberofsegments;

  int num_halfedges = num_e*2 - num_ch_v;
  assert(num_halfedges > 0);
  edges.reserve(num_halfedges);
  Edge * edge_end = edges.data();
  int *tptr = tout.trianglelist;
  int *nptr = tout.neighborlist;
  unsigned edges_on_ch = 0;
  for (unsigned i=0; i<num_t; ++i) {
    for (unsigned j=0; j<3; ++j) {
      Edge *buddy = NULL;
      if (*nptr >= 0) {
        int our_index_in_neighbor = tnidx_by_nidx(tout.neighborlist + 3*(*nptr), i);
        buddy = &edges[3*(*nptr) + our_index_in_neighbor];
      } else {
        ++edges_on_ch;
      }
      int edge_points_to_vertex_idx = tptr[ prev_edge_offset[j] ];
      Vertex *edge_points_to_vertex = &vertices[edge_points_to_vertex_idx];
      edges.emplace_back(Edge(edge_end+next_edge_offset[j], edge_end+prev_edge_offset[j], buddy, edge_points_to_vertex, edges.size()));
      ++nptr;
    }
    edge_end += 3;
    tptr += 3;
  }
  assert(num_ch_v == edges_on_ch);
  assert(tptr == tout.trianglelist + 3*num_t);
  assert(nptr == tout.neighborlist + 3*num_t);
  assert(edge_end == &edges[num_halfedges]);

  assert((edges.size() + num_ch_v) % 2 == 0);
  assert(num_t == 1 + (edges.size() + num_ch_v)/2 - vertices.size());

  // num_faces = num_t;
  unsigned number_of_ch_vertices = num_ch_v;
  return std::make_pair(std::move(edges), number_of_ch_vertices);
}

/** Triangulare the pointset and create the DECL
 */
std::pair<FixedVector<Edge>, unsigned>
DECL::
decl_triangulate(VertexList& vertices) {
  struct triangulateio tin, tout;
  memset(&tin, 0, sizeof(tin));
  memset(&tout, 0, sizeof(tout));

  decl_triangulate_prepare(vertices, tin);

  /* triangle options:
   * //N: no nodes output
   * //P: no poly output
   * Q: quiet
   * c: triangulate the convex hull
   * n: output neighbors list
   * // p: operate on a PSLG, with vertices, segments, and more.  We want
   *  segments part of that.
   * z: start indexes at 0
   */
  char trioptions[] = "Qcnz";
  triangulate(trioptions, &tin, &tout, NULL);
  auto res = decl_triangulate_process(vertices, tout);

  my_free_c(tin.pointlist);
  my_free_c(tin.segmentlist);
  my_free_c(tout.trianglelist);
  my_free_c(tout.neighborlist);
  my_free_c(tout.pointlist);
  my_free_c(tout.pointmarkerlist);
  my_free_c(tout.segmentlist);
  my_free_c(tout.segmentmarkerlist);
  return res;
}


/** unconstrain all edges for which this is possible. */
void
DECL::
unconstrain_all() {
  //DBG_FUNC_BEGIN(DBG_GENERIC);
  DBG_INDENT_INC();
  std::vector<Edge*> shuffled_edges = working_set.my_edges;
  std::shuffle(std::begin(shuffled_edges), std::end(shuffled_edges), random_engine);

  int old_num_faces = num_faces;
  for (Edge * const e : shuffled_edges) {
    /* Only check one edge of each half-edge-pair */
    if (e->opposite == NULL) continue;
    if (e->opposite < e) continue;

    if (! e->can_unconstrain()) continue;
    e->unconstrain();
    --num_faces;
  }
  //DBG(DBG_GENERIC) << "Now have " << num_faces << "/" << old_num_faces << " faces";
  DBG_INDENT_DEC();
  //DBG_FUNC_END(DBG_GENERIC);
}

/** Mark all edges as constrained again, resetting everything. */
void
DECL::
reset_constraints() {
  //DBG_FUNC_BEGIN(DBG_GENERIC);
  DBG_INDENT_INC();
  assert_hole_shooting_reset();
  for (Edge * const e : working_set.my_edges) {
    e->reset_all_constraints();
  }
  num_faces = working_set.num_faces_mine_constrained;
  DBG_INDENT_DEC();
  //DBG_FUNC_END(DBG_GENERIC);
}


/* Returns the next face cw of e, even if v is on the CH.
 *
 * e needs to be a constrained edge.
 */
Edge *
DECL::
get_next_face_cw_around_vertex(Edge *e) const {
  assert(e->is_constrained);

  //DBG(DBG_GENERIC) << " Looping around " << *e;
  Edge* r = e->next_constrained->opposite;
  if (!r) {
    r = e;
    //DBG(DBG_GENERIC) << "   r is NULL, restarting with " << *r;
    while (r->opposite) {
      r = r->opposite->prev;
      //DBG(DBG_GENERIC) << "   r now is " << *r;
      assert(r != e);
    }
  }
  assert(r->v == e->v);
  return r;
}

/* Returns the next face ccw of e.  e must not be on the boundary.
 *
 * e needs to be a constrained edge.
 */
Edge *
DECL::
get_next_face_ccw(Edge *e) const {
  assert(e->is_constrained);
  assert(e->opposite);

  Edge* r = e->opposite->prev_constrained;
  assert(r->v == e->v);
  return r;
}

/** Save the decomposition of the edges in a working set.
 */
DECL::SavedDecomposition::
SavedDecomposition(const WorkingSet& ws, int num_faces) :
  saved_num_faces(num_faces)
{
  edge_content.reserve(ws.my_edges.size());
  for (Edge * const e : ws.my_edges) {
    edge_content.emplace_back(*e);
  }
}

/** Re-inject a saved decomposition.
 */
void
DECL::
reinject_saved_decomposition(SavedDecomposition&& saved_decomposition) {
  assert(saved_decomposition.edge_content.size() == working_set.my_edges.size());

  unsigned size = working_set.my_edges.size();
  for (unsigned i = 0; i<size; ++i) {
    std::swap(*working_set.my_edges[i], saved_decomposition.edge_content[i]);
  }
  num_faces = saved_decomposition.saved_num_faces;
}


#if 0
/** Save the state of a region
 */
DECL::SavedState::
SavedState(std::vector<Edge*>&& my_edges_, int num_my_triangles_, int num_faces_mine_constrained_, int num_faces_)
  : Base(std::forward<std::vector<Edge*>>(my_edges_), num_my_triangles_, num_faces_mine_constrained_, num_faces_)
{
  edge_content.reserve(my_edges.size());
  for (Edge * const e : my_edges) {
    edge_content.emplace_back(*e);
  }
}

#endif



/** Find a simply-connected set of vertices, and make a hole
 */
void
DECL::
shoot_hole(unsigned size, int num_iterations, int max_recurse) {
  DBG_INDENT_INC();

  assert_hole_shooting_reset();

  #if 0
  shoot_hole_select_vertices(size);
  SavedState state(edges, halfedges_to_remove, faces_removed, vertices_to_remove.size() + removal_boundary_vertices.size(), removal_boundary_vertices.size() + vertices_to_remove_on_ch );

  /* in-place recursion here */
  find_convex_decomposition(num_iterations, state.num_faces, max_recurse);
  xx
  #endif

  DBG_INDENT_DEC();
}

void
DECL::
shoot_holes(int max_recurse) {
  DBG_INDENT_INC();

  int hole_size = working_set.num_my_triangles;
  while (hole_size >= 10) {
    hole_size = int(std::pow(double(hole_size), 2./3));
    int number_of_decompositions_per_hole = hole_size;
    int number_of_hole_punches = working_set.num_my_triangles/hole_size * NUMBER_OF_HOLE_PUNCHES_SCALE;

    DBG(DBG_GENERIC)
      << "Calling shoot_hole " << number_of_hole_punches << " times"
      << "; hole_size: " << hole_size
      << "; number_of_decompositions_per_hole: " << number_of_decompositions_per_hole
      << "; max_recurse: " << max_recurse;
    for (int i=0; i<number_of_hole_punches; ++i) {
      DEBUG_STMT({
        if ( (number_of_decompositions_per_hole >  100 && i % 100 == 0)
          || (number_of_decompositions_per_hole <= 100 && i % 500 == 0)
              ) {
          DBG(DBG_GENERIC) << "i: " << i << "; current num faces: " << num_faces;
        }
      });
      shoot_hole(hole_size, number_of_decompositions_per_hole, max_recurse);
      assert_hole_shooting_reset();
    }
  };

  DBG_INDENT_DEC();
}
void
DECL::
find_convex_decomposition(unsigned num_iterations, unsigned initial_num_faces_to_beat, unsigned max_recurse) {
  DBG_FUNC_BEGIN(DBG_GENERIC);

  unsigned num_faces_to_beat = initial_num_faces_to_beat ? initial_num_faces_to_beat : num_faces;
  bool have_solution = false;
  bool current_is_best = false;
  SavedDecomposition best = SavedDecomposition(working_set, num_faces);
  unsigned iter = 0;
  while (1) {
    unconstrain_all();
    if (max_recurse > 0) {
      shoot_holes(max_recurse-1);
    }

    current_is_best = (num_faces < num_faces_to_beat);
    if (current_is_best) {
      DBG(DBG_GENERIC) << "Iteration " << iter << "/" << num_iterations << ": This solution: " << num_faces << "; NEW BEST; previous best: " << num_faces_to_beat;
      have_solution = true;
      num_faces_to_beat = num_faces;
      best = SavedDecomposition(working_set, num_faces);
    } else {
      DBG(DBG_GENERIC) << "Iteration " << iter << "/" << num_iterations << ": This solution: " << num_faces << "; current best: " << num_faces_to_beat;
    }
    ++iter;

    if (iter >= num_iterations) break;
    DBG(DBG_GENERIC) << "Resetting constraints";
    reset_constraints();
  }

  if (have_solution) {
    DBG(DBG_GENERIC) << "Done " << num_iterations << " iterations.  Best now is " << num_faces_to_beat << " from " << initial_num_faces_to_beat;
  } else {
    DBG(DBG_GENERIC) << "Done " << num_iterations << " iterations.  We failed to improve on the bound of " << num_faces_to_beat << " faces";
  };
  if (!current_is_best) {
    DBG(DBG_GENERIC) << "Re-injecting saved decomposition with " << best.saved_num_faces << " faces";
    reinject_saved_decomposition(std::move(best));
  };
  assert_hole_shooting_reset();
  assert(num_faces_to_beat == num_faces);

  DBG_FUNC_END(DBG_GENERIC);
}


/** Runs validity checks for all the edges */
#ifndef NDEBUG
void
DECL::
assert_valid() const {
  for (const auto &e : all_edges) {
    e.assert_valid();
  }
}
#endif

/** Write the (constraint) segments to the output stream in obj format
 *
 * If we have a vertex list, also print those.
 */
void
DECL::
write_obj_segments(bool dump_vertices, std::ostream &o) const {
  if (dump_vertices) {
    for (const auto &v : all_vertices) {
      o << "v " << v.x << " " << v.y << " 0" << std::endl;
    }
  }
  for (const auto &e : all_edges) {
    if (!e.is_constrained) continue;
    if (e.opposite && e.opposite < &e) continue;

    int tail_idx = (e.get_tail() - all_vertices.data())+1;
    int head_idx = (e.v - all_vertices.data())+1;

    o << "l "
      << tail_idx
      << " "
      << head_idx
      << std::endl;
  }
}

std::ostream& operator<<(std::ostream& os, const Vertex& v) {
  os << "v"
     << v.idx()
     << "(" << v.x << "; " << v.y << ")";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Edge& e) {
  os << "e"
     << e.idx()
     << "(" << *e.get_tail()
     << "->" << *e.v;
  os << "; n/p/o";
  if (e.is_constrained) {
    os << "/N/P";
  }
  os << ": ";

  os << (e.next->idx())
     << "/" << (e.prev->idx())
     << "/" << (e.opposite ? e.opposite->idx() : -1);
  if (e.is_constrained) {
    os << "/" << (e.next_constrained->idx())
       << "/" << (e.prev_constrained->idx());
  }
  os << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const DECL& d) {
  os << "DECL" << std::endl;
  os << " #edges: " << d.all_edges.size() << std::endl;
  for (unsigned i=0; i<d.all_edges.size(); ++i) {
    const Edge &e = d.all_edges[i];

    os << "   edge #" << e
       << std::endl;
  }
  return os;
}

/* vim: set fdm=marker: */
