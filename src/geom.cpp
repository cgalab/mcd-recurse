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
DECL(VertexList& vertices)
: vertex_base_ptr(&vertices[0])
, num_vertices(vertices.size())
{
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
decl_triangulate_process(VertexList& vertices, const struct triangulateio& tout) {
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

  num_triangles = num_t;
  num_faces = num_t;
}

/** Triangulare the pointset and create the DECL
 */
void
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
  DBG_FUNC_BEGIN(DBG_GENERIC);
  std::vector<int> edge_index(edges.size());
  std::iota (std::begin(edge_index), std::end(edge_index), 0);

  std::shuffle(std::begin(edge_index), std::end(edge_index), random_engine);

  for (const int idx : edge_index) {
    Edge& e = edges[idx];

    /* Only check one edge of each half-edge-pair */
    if (e.opposite == NULL) continue;
    if (e.opposite < &e) continue;

    if (! e.can_unconstrain()) continue;
    e.unconstrain();
    --num_faces;
  }
  DBG_FUNC_END(DBG_GENERIC);
}

/** Mark all edges as constrained again, resetting everything. */
void
DECL::
reset_constraints() {
  DBG_FUNC_BEGIN(DBG_GENERIC);
  assert_hole_shooting_reset();
  for (auto &e : edges) {
    e.reset_all_constraints();
  }
  num_faces = num_triangles;
  DBG_FUNC_END(DBG_GENERIC);
}


/** Pick a random vertex and select size many vertices around it to form a hole
 */
void
DECL::
shoot_hole_select_vertices(unsigned size) {
  DBG_FUNC_BEGIN(DBG_GENERIC);
  /* We pick a random vertex by picking a random edge that points to it.
   * This induces a bias towards higher-degree vertices, but that's
   * probably actually good.
   */
  vertices_to_remove.reserve(2*size); /* Slightly larger since we might need
                                         the extra bit during the BFS, and to
                                         take care of the case of holes which
                                         we then have to fill up */

  Edge *random_edge = &*random_element(std::begin(edges), std::end(edges), random_engine);
  /* Find a local constraint edge */
  while (! random_edge->is_constrained) {
    random_edge = random_edge->next->opposite;
  }
  vertices_to_remove.emplace_back(random_edge);
  random_edge->v->vertex_to_be_removed = true;

  DBG(DBG_GENERIC) << " At start of find-loop.";
  unsigned find_neighbors_of_vertex_idx = 0;
  while (vertices_to_remove.size() < size) {
    assert(find_neighbors_of_vertex_idx < vertices_to_remove.size());

    Edge *vertex_to_remove = vertices_to_remove[find_neighbors_of_vertex_idx];
    assert(vertex_to_remove);
    assert(vertex_to_remove->is_constrained);

    DBG(DBG_GENERIC) << " At head of find-loop body.  vertex_to_remove is: " << *vertex_to_remove;

    for (auto e_it = AroundVertexFacesIterator(vertex_to_remove); *e_it && vertices_to_remove.size() < size; ++e_it) {
      Edge * vertex_candidate = e_it->next_constrained;
      DBG(DBG_GENERIC) << "  Iterating around vertex.  Current vertex_candidate " << *vertex_candidate;
      if (! vertex_candidate->v->vertex_to_be_removed) {
        /* Which is not yet marked for removal.  Do that now. */
        vertices_to_remove.emplace_back(vertex_candidate);
        vertex_candidate->v->vertex_to_be_removed = true;
        DBG(DBG_GENERIC) << "   Added vertex " << *vertex_candidate;
      }
    }
    ++find_neighbors_of_vertex_idx;
  }
  DBG_FUNC_END(DBG_GENERIC);
}

/** Identify all triangles in the face left of e_start.
 */
void
DECL::
shoot_hole_list_triangles_in_face(Edge * const e_start) {
  DBG_FUNC_BEGIN(DBG_GENERIC);
  if (e_start->triangle_to_be_removed) {
    DBG(DBG_GENERIC) << "Face already marked for removal.";
    DBG(DBG_GENERIC) << "  t: " << e_start->v->idx << ", " << e_start->next->v->idx << ", " << e_start->prev->v->idx;
  } else {
    DBG(DBG_GENERIC) << "Marking Face for removal.";
    Edge *e = e_start;
    int cnt_triangles = 0;
    int cnt_constrained_halfedges = 0;
    do {
      assert(e->next->next->next == e);

      DBG(DBG_GENERIC) << "  e: " << *e;
      DBG(DBG_GENERIC) << "  t: " << e->v->idx << ", " << e->next->v->idx << ", " << e->prev->v->idx;
      if (!e->triangle_to_be_removed) {
        halfedges_to_remove.emplace_back(e);
        halfedges_to_remove.emplace_back(e->next);
        halfedges_to_remove.emplace_back(e->prev);
        e->triangle_to_be_removed = true;
        e->next->triangle_to_be_removed = true;
        e->prev->triangle_to_be_removed = true;
        DEBUG_STMT(++cnt_triangles);
      }
      if (!e->v->vertex_to_be_removed) {
        if (!e->v->vertex_on_removal_boundary) {
          DBG(DBG_GENERIC) << "Touched a boundary vertex " << *e->v;
          e->v->vertex_on_removal_boundary = true;
          removal_boundary_vertices.emplace_back(e);
        } else {
          DBG(DBG_GENERIC) << "Touched an already marked boundary vertex " << *e->v;
        }
      } else {
        DBG(DBG_GENERIC) << "Touched a vertex to be removed " << *e->v;
      }

      e = e->next;
      if (!e->is_constrained) {
        e = e->opposite;
        assert(e);
      } else {
        DEBUG_STMT(++cnt_constrained_halfedges);
      };
    } while (e != e_start);
    DBG(DBG_GENERIC) << "Hit " << cnt_triangles << " triangles.";
    DBG(DBG_GENERIC) << "Hit " << cnt_constrained_halfedges << " constrained halfedges.";
    assert(cnt_triangles == cnt_constrained_halfedges-2);
    ++faces_removed;
  }
  DBG_FUNC_END(DBG_GENERIC);
}


/** Identify all faces that are incident to a vertex
 *
 * and all all their triangles' halfedges to the halfedges_to_remove vector.
 */
void
DECL::
shoot_hole_identify_affected_elements_around_vertex(Edge* const e_vertex) {
  DBG_FUNC_BEGIN(DBG_GENERIC);
  DBG(DBG_GENERIC) << "Working on vertex pointed to by " << *e_vertex;
  assert(e_vertex->v->vertex_to_be_removed);
  assert(!e_vertex->v->vertex_on_removal_boundary);
  auto e_it = AroundVertexFacesIterator(e_vertex);
  for (; *e_it; ++e_it) {
    shoot_hole_list_triangles_in_face(*e_it);
  }

  if (e_it.hit_ch()) {
    ++vertices_to_remove_on_ch;
  }
  DBG_FUNC_END(DBG_GENERIC);
}

/** Identify vertices which we think are on the removal boundary but which are actually fully enclosed by removed faces.
 */
void
DECL::
shoot_hole_identify_enclosed_boundary_vertices() {
  DBG_FUNC_BEGIN(DBG_GENERIC);
  /* Check if any of the "boundary" vertices have been entirely enclosed by faces we want to remove */
  auto e_bd_v_it = removal_boundary_vertices.begin();
  while (e_bd_v_it != removal_boundary_vertices.end()) {
    Edge * const e_bd_v = *e_bd_v_it;

    DBG(DBG_GENERIC) << "v: " << e_bd_v->v;
    assert(!e_bd_v->v->vertex_to_be_removed);
    assert(e_bd_v->v->vertex_on_removal_boundary);

    assert(e_bd_v->is_constrained);
    bool enclosed = true;
    auto e_it = AroundVertexFacesIterator(e_bd_v);
    for (; *e_it; ++e_it) {
      if (!e_it->triangle_to_be_removed) {
        enclosed = false;
        break;
      }
    }
    if (enclosed) {
      DBG(DBG_GENERIC) << "Found enclosed vertex " << *e_bd_v->v;
      vertices_to_remove_on_ch += e_it.hit_ch();
      vertices_to_remove.emplace_back(e_bd_v);
      e_bd_v->v->vertex_to_be_removed = true;
      e_bd_v->v->vertex_on_removal_boundary = false;

      shoot_hole_identify_affected_elements_around_vertex(e_bd_v);

      std::swap(*e_bd_v_it, removal_boundary_vertices.back());
      removal_boundary_vertices.pop_back();
      continue;
    }
    ++e_bd_v_it;
  }

  DBG_FUNC_END(DBG_GENERIC);
}

/** Identify all elements that are in any face incident to a vertex of vertices_to_remove
 */
void
DECL::
shoot_hole_identify_affected_elements() {
 /* We iterate over all the vertices.
  *
  * Around each vertex, we iterate over all incident faces that have not yet
  * been handled.
  *
  * In each face we identify all relevant triangles.
  *
  * While we do this, we also build up a list of vertices that are not to be
  * removed but which are incident to any of the faces we are removing.
  *
  * If there is more than one boundary, this means we created a non-simply
  * connected set of faces to remove.  In this case, we also want to remove
  * any faces (and vertices) that are enclosed by the removed set.
  *
  * We want to produce several things here:
  *  - The list of vertices on the boundary.
  *  - Add any extra vertices to remove to the vertices_to_remove vector.
  *  - The number of faces removed.
  *  - A "serialized" copy of state within the hole so that it can
  *    be re-injected at a later time if we decide we want to keep it
  *    after all.
  *  - A list of edge locations that held that state so we can play-back
  *    said stored state, or a new decomposition should we prefer it.
  */
  DBG_FUNC_BEGIN(DBG_GENERIC);
  for (auto e_vertex_it : vertices_to_remove) {
    shoot_hole_identify_affected_elements_around_vertex(e_vertex_it);
  }
  shoot_hole_identify_enclosed_boundary_vertices();

  DBG(DBG_GENERIC) << "Removing " << faces_removed << " faces.";
  DBG(DBG_GENERIC) << "Removing " << vertices_to_remove.size() << " vertices.";
  DBG(DBG_GENERIC) << "Of those, " << vertices_to_remove_on_ch << " are on the CH.";
  DBG(DBG_GENERIC) << "Removing " << halfedges_to_remove.size()/3 << " triangles.";
  DBG(DBG_GENERIC) << "Removing " << halfedges_to_remove.size() << " halfedges.";
  DBG(DBG_GENERIC) << "Found " << removal_boundary_vertices.size() << " boundary vertices:";


  /* We deal with a planar graph, so v-e+f == 2, including the outer face.
   *
   * This also must hold for the area we want to remove.  So we can figure out
   * whether it has enclosed a region without marking it for removal by a
   * counting argument.
   *
   * v is simply the number of all vertices, i.e. the ones on the removal boundary + the ones to remove.
   * f is the number of triangles, plus the outer face, plus any enclosed faces.
   * e is tricky, since we don't know it directly.
   *   Each internal triangulation half-edge *pairs* counts one towards e.
   *   And each triangulation half-edge *singleton* also counts one.  So
   *   how many singletons are there?  As many as are edges on the boundary of the
   *   area we want to remove.  Which is the same as the number of vertices
   *   on that boundary.  Note that this is the boundary of the area we want to
   *   remove, not just the boundary with the remaining DECL.  So we have to
   *   count both removal_boundary_vertices *and* vertices_to_remove_on_ch.
   */
  {
    assert(halfedges_to_remove.size() % 3 ==0);

    int total_num_v  = vertices_to_remove.size() + removal_boundary_vertices.size();
    int total_num_2e = halfedges_to_remove.size() + removal_boundary_vertices.size() + vertices_to_remove_on_ch;
    int total_num_f_accounted_for = halfedges_to_remove.size()/3;
    assert(total_num_2e % 2 == 0);
    int inner_faces = 1 + total_num_2e/2 - total_num_v - total_num_f_accounted_for;
    DBG(DBG_GENERIC) << "Number of enclosed faces: " << inner_faces;
  }

  DBG_FUNC_END(DBG_GENERIC);
}


/** Find a simply-connected set of vertices, and make a hole
 */
void
DECL::
shoot_hole(unsigned size) {
  DBG_FUNC_BEGIN(DBG_GENERIC);

  assert_hole_shooting_reset();

  shoot_hole_select_vertices(size);
  shoot_hole_identify_affected_elements();

  DBG(DBG_GENERIC) << " Iterating over vertices marked for removal:";
  for (auto e : vertices_to_remove) {
    DBG(DBG_GENERIC) << *e;
    e->v->vertex_to_be_removed = false;
  }
  DBG_FUNC_END(DBG_GENERIC);
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
 *
 * If we have a vertex list, also print those.
 */
void
DECL::
write_obj_segments(const VertexList * vertices, std::ostream &o) const {
  if (vertices) {
    for (const auto &v : *vertices) {
      o << "v " << v.x << " " << v.y << " 0" << std::endl;
    }
  }
  for (const auto &e : edges) {
    if (!e.is_constrained) continue;
    if (e.opposite && e.opposite < &e) continue;

    int tail_idx = (e.get_tail() - vertex_base_ptr)+1;
    int head_idx = (e.v - vertex_base_ptr)+1;

    o << "l "
      << tail_idx
      << " "
      << head_idx
      << std::endl;
  }
}

std::ostream& operator<<(std::ostream& os, const Vertex& v) {
  os << "v"
#ifndef NDEBUG
     << v.idx
#endif
     << "(" << v.x << "; " << v.y << ")";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Edge& e) {
  os << "e"
#ifndef NDEBUG
     << e.idx
#endif
     << "(" << *e.get_tail()
     << "->" << *e.v;
#ifndef NDEBUG
  os << "; n/p/o";
  if (e.is_constrained) {
    os << "/N/P";
  }
  os << ": ";

  os << (e.next->idx)
     << "/" << (e.prev->idx)
     << "/" << (e.opposite ? e.opposite->idx : -1);
  if (e.is_constrained) {
    os << "/" << (e.next_constrained->idx)
       << "/" << (e.prev_constrained->idx);
  }
#endif
  os << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const DECL& d) {
  os << "DECL" << std::endl;
  os << " #edges: " << d.edges.size() << std::endl;
  for (unsigned i=0; i<d.edges.size(); ++i) {
    const Edge &e = d.edges[i];

    os << "   edge #" << e
       << std::endl;
  }
  return os;
}
