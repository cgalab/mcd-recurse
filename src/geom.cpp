#include "geom.h"
#include "triangle.h"

#include <numeric>

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

/** Initialize the DECL with the vertices and a triangulation of their CH */
DECL::
DECL(const VertexList& vertices)
: vertex_base_ptr(&vertices[0]) {
  decl_triangulate(vertices);
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
void
DECL::
decl_triangulate_process(const VertexList& vertices, const struct triangulateio& tout) {
  const int next_edge_offset[] = {1,2,0};
  const int prev_edge_offset[] = {2,0,1}; /* also same as vertex index that an edge points towards */

  int num_t = tout.numberoftriangles;
  int num_e = tout.numberofedges;
  int num_ch_v = tout.numberofsegments;

  int num_halfedges = num_e*2 - num_ch_v;
  edges.reserve(num_halfedges);
  Edge * edge_end = &edges[0];
  int *tptr = tout.trianglelist;
  int *nptr = tout.neighborlist;
  int edges_on_ch = 0;
  for (int i=0; i<num_t; ++i) {
    for (int j=0; j<3; ++j) {
      Edge *buddy = NULL;
      if (*nptr >= 0) {
        int our_index_in_neighbor = tnidx_by_nidx(tout.neighborlist + 3*(*nptr), i);
        buddy = &edges[3*(*nptr) + our_index_in_neighbor];
      } else {
        ++edges_on_ch;
      }
      int edge_points_to_vertex_idx = tptr[ prev_edge_offset[j] ];
      const Vertex *edge_points_to_vertex = &vertices[edge_points_to_vertex_idx];
      edges.emplace_back(Edge(edge_end+next_edge_offset[j], edge_end+prev_edge_offset[j], buddy, edge_points_to_vertex));
      ++nptr;
    }
    edge_end += 3;
    tptr += 3;
  }
  assert(num_ch_v == edges_on_ch);
  assert(tptr == tout.trianglelist + 3*num_t);
  assert(nptr == tout.neighborlist + 3*num_t);
  assert(edge_end == &edges[num_halfedges]);

  num_faces = num_t;
}

/** Triangulare the pointset and create the DECL
 */
void
DECL::
decl_triangulate(const VertexList& vertices) {
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
  decl_triangulate_process(vertices, tout);

  my_free_c(tin.pointlist);
  my_free_c(tin.segmentlist);
  my_free_c(tout.trianglelist);
  my_free_c(tout.neighborlist);
  my_free_c(tout.pointlist);
  my_free_c(tout.pointmarkerlist);
  my_free_c(tout.segmentlist);
  my_free_c(tout.segmentmarkerlist);
}


/** unconstrain all edges for which this is possible. */
void
DECL::
unconstrain_all() {
  std::vector<int> edge_index(edges.size());
  std::iota (std::begin(edge_index), std::end(edge_index), 0);

  std::shuffle(std::begin(edge_index), std::end(edge_index), random_engine);

  for (const int idx : edge_index) {
    Edge& e = edges[idx];

    /* Only check one edge of each half-edge-pair */
    if (e.get_opposite() == NULL) continue;
    if (e.get_opposite() < &e) continue;

    if (! e.can_unconstrain()) continue;
    e.unconstrain();
    --num_faces;
  }
}

/** Runs validity checks for all the edges */
#ifndef NDEBUG
void
DECL::
assert_valid() const {
  for (const auto &e : edges) {
    e.assert_valid();
  }
}
#endif

/** Write the (constraint) segments to the output stream in obj format
 */
void
DECL::
write_obj_segments(std::ostream &o) const {
  for (const auto &e : edges) {
    if (!e.get_is_constrained()) continue;
    if (e.get_opposite() && e.get_opposite() < &e) continue;

    int tail_idx = (e.get_tail() - vertex_base_ptr)+1;
    int head_idx = (e.get_head() - vertex_base_ptr)+1;

    o << "l "
      << tail_idx
      << " "
      << head_idx
      << std::endl;
  }
}

std::ostream& operator<<(std::ostream& os, const DECL& d) {
  os << "DECL" << std::endl;
  os << " #edges: " << d.edges.size() << std::endl;
  for (unsigned i=0; i<d.edges.size(); ++i) {
    const Edge &e = d.edges[i];

    os << "   edge #" << i
       << ": points-to: ("
#ifndef NDEBUG
       << e.v->idx << ": "
#endif
       << e.v->x << ", " << e.v->y << ")"
       << "; n/p/o: " << (e.next - &d.edges[0])
       << "; " << (e.prev - &d.edges[0])
       << "; " << (e.opposite ? (e.opposite - &d.edges[0]) : -1 )
       << std::endl;
  }
  return os;
}
