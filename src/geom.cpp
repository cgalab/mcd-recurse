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

/** Initialize the DECL with the vertices and a triangulation of their CH */
DECL::
DECL(std::shared_ptr<VertexList> vertices)
: all_vertices(vertices)
, num_vertices(vertices->size())
{
  decl_triangulate(*vertices);
}

/** Initialize the DECL with the saved structure, ignoring constraints */
DECL::
DECL(std::shared_ptr<VertexList> vertices, const SavedState &state)
: all_vertices(vertices)
, num_vertices(state.num_vertices)
, num_vertices_on_boundary(state.num_vertices_on_boundary)
, num_triangles(state.triangle_vertices.size() / 3)
, num_faces(state.triangle_vertices.size() / 3)
{
  DBG_INDENT_INC();
  //DBG(DBG_GENERIC) << "Setting up a DECL with " << num_vertices << " vertices";
  int num_halfedges = state.triangle_vertices.size();
  assert(num_halfedges % 3 == 0);
  int num_t = num_halfedges/3;
  edges.reserve(num_halfedges);

  Edge * edge_end = &edges[0];
  Vertex * const * vptr = state.triangle_vertices.data();
  int const * bptr = state.buddy.data();
  int edges_on_boundary = 0;
  DBG(DBG_STATERESTORE) << "num_t: " << num_t;
  for (int i=0; i<num_t; ++i ) {
    DBG(DBG_STATERESTORE) << " t: " << i;
    DBG(DBG_STATERESTORE) << "  v   : " << *vptr[0] << "; " << *vptr[1] << "; " << *vptr[2];
    DBG(DBG_STATERESTORE) << "  bidx: " <<  bptr[0] << "; " <<  bptr[1] << "; " <<  bptr[2];

    for (int j=0; j<3; ++j, ++bptr) {
      Edge *buddy = NULL;
      if (*bptr >= 0) {
        buddy = &edges[*bptr];
      } else {
        DBG(DBG_STATERESTORE) << "  Edge " << j << " is on boundary";
        ++edges_on_boundary;
      }

      Vertex* edge_points_to_vertex = vptr[ prev_edge_offset[j] ];
      edges.emplace_back(Edge(edge_end+next_edge_offset[j], edge_end+prev_edge_offset[j], buddy, edge_points_to_vertex, edges.size()));
    }
    edge_end += 3;
    vptr += 3;
  }
  DBG(DBG_STATERESTORE) << "edges_on_boundary             : " << edges_on_boundary;
  DBG(DBG_STATERESTORE) << "state.num_vertices_on_boundary: " << state.num_vertices_on_boundary;
  assert(state.num_vertices_on_boundary == edges_on_boundary);
  assert(vptr == state.triangle_vertices.data() + num_halfedges);
  assert(bptr == state.buddy.data() + num_halfedges);
  assert(edge_end == &edges[num_halfedges]);

  assert_hole_shooting_reset();
  DBG_INDENT_DEC();
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
  num_vertices_on_boundary = num_ch_v;
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
  //DBG_FUNC_BEGIN(DBG_GENERIC);
  DBG_INDENT_INC();
  std::vector<int> edge_index(edges.size());
  std::iota (std::begin(edge_index), std::end(edge_index), 0);

  std::shuffle(std::begin(edge_index), std::end(edge_index), random_engine);
  int old_num_faces = num_faces;

  for (const int idx : edge_index) {
    Edge& e = edges[idx];

    /* Only check one edge of each half-edge-pair */
    if (e.opposite == NULL) continue;
    if (e.opposite < &e) continue;

    if (! e.can_unconstrain()) continue;
    e.unconstrain();
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
  for (auto &e : edges) {
    e.reset_all_constraints();
  }
  num_faces = num_triangles;
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


/** Pick a random vertex and select size many vertices around it to form a hole
 */
void
DECL::
shoot_hole_select_vertices(unsigned size) {
  DBG_FUNC_BEGIN(DBG_SHOOTHOLE);
  vertices_to_remove.reserve(2*size); /* Slightly larger since we might need
                                         the extra bit during the BFS, and to
                                         take care of the case of holes which
                                         we then have to fill up */

    /* We pick a random vertex by picking a random edge that points to it.
     * Then we look for the next constraint edge.  We don't do this by
     * walking cyclicly, but by walking the edge vector.  We hope
     * to avoid some kind of bias here.
     */
  auto random_edge_it = random_element(std::begin(edges), std::end(edges), random_engine);
  DBG(DBG_SHOOTHOLE) << &*random_edge_it;
  DBG(DBG_SHOOTHOLE) << &*std::end(edges);
  /* Find a local constraint edge */
  while (! random_edge_it->is_constrained) {
    ++random_edge_it;
    if (random_edge_it == std::end(edges)) {
      random_edge_it = std::begin(edges);
    };
  }
  Edge *random_edge = &*random_edge_it;
  vertices_to_remove.emplace_back(random_edge);
  random_edge->v->vertex_to_be_removed = true;

  DBG(DBG_SHOOTHOLE) << " At start of find-loop.";
  unsigned find_neighbors_of_vertex_idx = 0;
  while (vertices_to_remove.size() < size) {
    assert(find_neighbors_of_vertex_idx < vertices_to_remove.size());

    Edge *vertex_to_remove = vertices_to_remove[find_neighbors_of_vertex_idx];
    assert(vertex_to_remove);
    assert(vertex_to_remove->is_constrained);

    DBG(DBG_SHOOTHOLE) << " At head of find-loop body.  vertex_to_remove is: " << *vertex_to_remove;

    for (auto e_it = AroundVertexFacesIterator(vertex_to_remove); *e_it && vertices_to_remove.size() < size; ++e_it) {
      Edge * vertex_candidate = e_it->next_constrained;
      DBG(DBG_SHOOTHOLE) << "  Iterating around vertex.  Current vertex_candidate " << *vertex_candidate;
      if (! vertex_candidate->v->vertex_to_be_removed) {
        /* Which is not yet marked for removal.  Do that now. */
        vertices_to_remove.emplace_back(vertex_candidate);
        vertex_candidate->v->vertex_to_be_removed = true;
        DBG(DBG_SHOOTHOLE) << "   Added vertex " << *vertex_candidate;
      }
    }
    ++find_neighbors_of_vertex_idx;
  }
  DBG_FUNC_END(DBG_SHOOTHOLE);
}

/** Identify all triangles in the face left of e_start.
 */
void
DECL::
shoot_hole_list_triangles_in_face(Edge * const e_start) {
  DBG_FUNC_BEGIN(DBG_SHOOTHOLE);
  if (e_start->triangle_to_be_removed) {
    DBG(DBG_SHOOTHOLE) << "Face already marked for removal.";
    DBG(DBG_SHOOTHOLE) << "  t: " << e_start->v->idx() << ", " << e_start->next->v->idx() << ", " << e_start->prev->v->idx();
  } else {
    DBG(DBG_SHOOTHOLE) << "Marking Face for removal.";
    Edge *e = e_start;
    int cnt_triangles = 0;
    int cnt_constrained_halfedges = 0;
    do {
      assert(e->next->next->next == e);

      DBG(DBG_SHOOTHOLE) << (e->is_constrained ? "E: (" : "  e(")
                         << "v" << e->get_tail()->idx() << "->" << "v" << e->v->idx() << ");"
                         << "  t: (" << e->v->idx() << ", " << e->next->v->idx() << ", " << e->prev->v->idx() << ")"
                         << (e->is_constrained ? "; constrained" : "");
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
          DBG(DBG_SHOOTHOLE) << "  " << *e->v << ": Touched a boundary vertex";
          e->v->vertex_on_removal_boundary = true;
          removal_boundary_vertices.emplace_back(e);
        } else {
          DBG(DBG_SHOOTHOLE) << "  Touched an already marked boundary vertex " << *e->v;
        }
      } else {
        DBG(DBG_SHOOTHOLE) << "  Touched a vertex to be removed " << *e->v;
      }

      e = e->next;
      if (!e->is_constrained) {
        e = e->opposite;
        assert(e);
      } else {
        DEBUG_STMT(++cnt_constrained_halfedges);
      };
    } while (e != e_start);
    DBG(DBG_SHOOTHOLE) << "Hit " << cnt_triangles << " triangles.";
    DBG(DBG_SHOOTHOLE) << "Hit " << cnt_constrained_halfedges << " constrained halfedges.";
    assert(cnt_triangles == cnt_constrained_halfedges-2);
    faces_removed += 1;
  }
  DBG_FUNC_END(DBG_SHOOTHOLE);
}


/** Identify all faces that are incident to a vertex
 *
 * and all all their triangles' halfedges to the halfedges_to_remove vector.
 */
void
DECL::
shoot_hole_identify_affected_elements_around_vertex(Edge* const e_vertex) {
  DBG_FUNC_BEGIN(DBG_SHOOTHOLE);
  DBG(DBG_SHOOTHOLE) << "Working on vertex " << *e_vertex->v << " pointed to by " << *e_vertex;
  assert(e_vertex->v->vertex_to_be_removed);
  assert(!e_vertex->v->vertex_on_removal_boundary);
  auto e_it = AroundVertexFacesIterator(e_vertex);
  for (; *e_it; ++e_it) {
    shoot_hole_list_triangles_in_face(*e_it);
  }

  if (e_it.hit_ch()) {
    DBG(DBG_SHOOTHOLE) << "We hit the CH during working on vertex " << *e_vertex->v;
    ++vertices_to_remove_on_ch;
  }
  DBG_FUNC_END(DBG_SHOOTHOLE);
}

/** Identify a start vertex to walk around the boundary
 *
 * We start with a vertex which is likely to be on the eventual boundary, such
 * as the one with minimal x-coordinate.  It might still not end up on the
 * eventual boundary if it is a CH vertex that is incident, on its CH edges, by
 * marked faces.  (With marked faces we mean faces marked to be removed.)
 * Note that there could be non-marked faces incident still, in which those would
 * be enclosed by marked faces and should also be removed.)
 */
Edge*
DECL::
shoot_hole_identify_boundary_vertices_start() {
  DBG_FUNC_BEGIN(DBG_SHOOTHOLE);
  Edge* res;

  /* Find a reasonable vertex to start with */
  bool tried_max = false;
  auto extreme_x_vertex_edge_it = std::min_element(removal_boundary_vertices.begin(), removal_boundary_vertices.end(),
    [](const Edge* a, const Edge* b) {return a->v->x < b->v->x;});
    /* Note on the use of min_element: We could do std::minmax_element to get
     * both the minimum and the maximum at the same time for still linear cost.
     * However, constants matter and it is somewhat unlikely we'll ever need
     * more than one extremal point.  So let's not.
     */

  /* Now verify we like that vertex. */
  while(1) {
    Edge* extreme_x_vertex_edge = *extreme_x_vertex_edge_it;
    DBG(DBG_SHOOTHOLE) << "extreme_x_vertex: " << *extreme_x_vertex_edge->v << " via " << **extreme_x_vertex_edge_it;

    /* Find out if the marked area around this vertex is a single continuous set.
     */

    int changes_seen = 0;
    bool in_marked_face = true;
    bool hit_ch = false;

    auto e_it = AroundVertexFacesStriclyClockwiseIterator(extreme_x_vertex_edge);
    assert(e_it->triangle_to_be_removed);
    for (; *e_it; ++e_it) {
      DBG(DBG_SHOOTHOLE) << **e_it;
      if (in_marked_face != e_it->triangle_to_be_removed) {
        ++changes_seen;
        in_marked_face = e_it->triangle_to_be_removed;
      }
    }
    changes_seen += (!in_marked_face);
    /* We started in a marked face, and the last face in our walk around
     * the vertex was not marked.  So we need to count one more. */
    DBG(DBG_SHOOTHOLE) << "Changes seen: " << changes_seen;
    assert(changes_seen % 2 == 0);

    /* If we saw exactly two changes, then this is a good candidate.  */
    if (changes_seen == 2) {
      /*
       * The only reason to disqualify it would be if it was on the CH, and the
       * two faces incident to the vertex and the CH were both marked for
       * removal.
       */
      if (!e_it.hit_ch()) break;
      DBG(DBG_SHOOTHOLE) << "We did hit the CH however.";
      DBG(DBG_SHOOTHOLE) << "  Are the edges before/after the CH marked for removal? "
        << (e_it.edge_before_ch()->triangle_to_be_removed ? "yes/" : "no/")
        << (e_it.edge_after_ch()->triangle_to_be_removed ? "yes." : "no.");
      if (!e_it.edge_before_ch()->triangle_to_be_removed) break;
      if (!e_it.edge_after_ch()->triangle_to_be_removed) break;
      DBG(DBG_SHOOTHOLE) << "So we don't like this one after all.";
    }

    /* Else we really don't want to start with this vertex.  Let's try the one with maximum x-coordinate instead. */
    if (tried_max) {
      /* We already tried that too.  We could handle this case and try to still
       * get it working, but at some point it's just not worth it.  So give up,
       * throw an exception and retry with different randomness. */
      res = NULL;
      goto done;
    }
    /* Try the max one */
    tried_max = true;
    extreme_x_vertex_edge_it = std::max_element(removal_boundary_vertices.begin(), removal_boundary_vertices.end(),
      [](const Edge* a, const Edge* b) {return a->v->x < b->v->x;});
  };
  DBG(DBG_SHOOTHOLE) << "Returning vertex " << *(*extreme_x_vertex_edge_it)->v << " via " << **extreme_x_vertex_edge_it;

  res = *extreme_x_vertex_edge_it;
 done:
  DBG_FUNC_END(DBG_SHOOTHOLE);
  return res;
}

/** Identify boundary vertices and enclosed vertices.
 *
 * We put the boundary vertices on removal_boundary_vertices,
 * and identify and mark as removed any enclosed regious.
 *
 * returns true if ok and false on failure.
 */
bool
DECL::
shoot_hole_identify_boundary_vertices() {
  DBG_FUNC_BEGIN(DBG_SHOOTHOLE);

  assert( std::all_of(removal_boundary_vertices.begin(), removal_boundary_vertices.end(),
    [](const Edge* e){return e->v->vertex_on_removal_boundary;}) );

  /* Find a reasonable vertex to start with */
  Edge* extreme_x_vertex_edge = shoot_hole_identify_boundary_vertices_start();
  if (!extreme_x_vertex_edge) {
    DBG_FUNC_END(DBG_SHOOTHOLE);
    return false;
  }

  /* Go to the removal boundary */
  std::vector<Edge *> boundary;

  Edge *e_start = extreme_x_vertex_edge;
  assert(e_start->triangle_to_be_removed);
  while (e_start->triangle_to_be_removed) e_start = get_next_face_cw_around_vertex(e_start);

  int num_vertices_on_boundary_list_handled = 0;
  Edge *e = e_start;
  do {
    assert(!e->v->vertex_on_outer_removal_boundary);
    assert(e->v->vertex_on_removal_boundary || e->v->vertex_to_be_removed);
    if (e->triangle_to_be_removed) {
      DBG(DBG_SHOOTHOLE) << "vertex " << *e->v << " is on the CH and the faces incident to vertex and the CH are all getting removed.  So should this vertex.";
      assert_vertex_is_on_ch(e);
      if (e->v->vertex_to_be_removed) {
        assert(!e->v->vertex_on_removal_boundary);
        assert(e->v->vertex_to_be_removed);

        DBG(DBG_SHOOTHOLE) << "We already had marked vertex to be removed.";
      } else {
        assert(e->v->vertex_on_removal_boundary);
        assert(!e->v->vertex_to_be_removed);

        e->v->vertex_on_removal_boundary = false;
        e->v->vertex_to_be_removed = true;
        vertices_to_remove.emplace_back(e);
        vertices_to_remove_on_ch += 1;

        num_vertices_on_boundary_list_handled += 1;
      }
    } else {
      assert(e->v->vertex_on_removal_boundary);
      assert(!e->v->vertex_to_be_removed);

      DBG(DBG_SHOOTHOLE) << "Adding vertex " << *e->v << " to boundary; via " << *e;
      DEBUG_STMT({
        bool on_ch = vertex_is_on_ch(e);
        DBG(DBG_SHOOTHOLE) << "vertex is on CH: " << (on_ch ? "yes" : "no");
      });
      e->v->vertex_on_outer_removal_boundary = true;
      boundary.push_back(e);

      num_vertices_on_boundary_list_handled += 1;

      while (!e->triangle_to_be_removed) e = get_next_face_cw_around_vertex(e);
    }

    e = e->prev_constrained;
    assert(e->triangle_to_be_removed);
    e = get_next_face_cw_around_vertex(e);
  } while (e != e_start);


  /* Now deal with vertices not on the outer boundary */
  int num_vertices_on_boundary_list_unaccouned = removal_boundary_vertices.size() - num_vertices_on_boundary_list_handled;
  assert(num_vertices_on_boundary_list_unaccouned >= 0);
  if (num_vertices_on_boundary_list_unaccouned == 0) {
    DBG(DBG_SHOOTHOLE) << "All boundary vertices accounted for";
  } else {
    DBG(DBG_SHOOTHOLE) << "Some boundary vertices unaccounted for.  Count: " << num_vertices_on_boundary_list_unaccouned;
    for (auto e_it : removal_boundary_vertices) {
      if (e_it->v->vertex_on_outer_removal_boundary) continue;
      if (!e_it->v->vertex_on_removal_boundary) continue;
      Edge *e_to_remove = &*e_it;
      DBG(DBG_SHOOTHOLE) << "Unaccounted for vertex " << *e_to_remove->v << "; via " << *e_to_remove;

      e_to_remove->v->vertex_on_removal_boundary = false;
      e_to_remove->v->vertex_to_be_removed = true;
      vertices_to_remove.emplace_back(e_to_remove);

      int old_vertices_to_remove_on_ch = vertices_to_remove_on_ch;
      shoot_hole_identify_affected_elements_around_vertex(e_to_remove);
      assert(vertices_to_remove_on_ch >= old_vertices_to_remove_on_ch);

      if (vertices_to_remove_on_ch != old_vertices_to_remove_on_ch) {
        assert(vertex_is_on_ch(e_to_remove));
        DBG(DBG_GENERIC) << "We have an accounted-for vertex on the CH.  Calling this one failed.";
        DBG_FUNC_END(DBG_SHOOTHOLE);
        return false;
      }
    }
  }

  for (auto i : removal_boundary_vertices) {
    assert(i->v->vertex_on_removal_boundary != i->v->vertex_to_be_removed);
    i->v->vertex_on_removal_boundary = false;
  }
  assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_on_removal_boundary == false;}) );

  DEBUG_STMT( {
    /* Ensure that an edge has vertex_on_outer_removal_boundary set if and only iff it is in the boundary vector. */
    assert( std::all_of(boundary.begin(), boundary.end(),
      [](const Edge* e){return e->v->vertex_on_outer_removal_boundary;}) );
    for (auto i : boundary) i->v->vertex_on_outer_removal_boundary = false;
    assert( std::all_of(edges.begin(), edges.end(), [](const Edge& e){return e.v->vertex_on_outer_removal_boundary == false;}) );
    for (auto i : boundary) i->v->vertex_on_outer_removal_boundary = true;
  });

  removal_boundary_vertices.swap(boundary);
  DBG_FUNC_END(DBG_SHOOTHOLE);
  return true;
}

/** Identify all elements that are in any face incident to a vertex of vertices_to_remove
 *
 * returns true if ok and false on failure.
 */
bool
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
  DBG_FUNC_BEGIN(DBG_SHOOTHOLE);
  bool res = false;

  for (auto e_vertex_it : vertices_to_remove) {
    shoot_hole_identify_affected_elements_around_vertex(e_vertex_it);
  }
  bool success = shoot_hole_identify_boundary_vertices();
  if (!success) {
    DBG(DBG_GENERIC | DBG_SHOOTHOLE) << "Hole shooting failed.";
    goto done;
  }

  DBG(DBG_SHOOTHOLE) << "Removing " << faces_removed << " faces.";
  DBG(DBG_SHOOTHOLE) << "Removing " << vertices_to_remove.size() << " vertices out of " << num_vertices;
  DBG(DBG_SHOOTHOLE) << "Of those, " << vertices_to_remove_on_ch << " are on the CH.";
  DBG(DBG_SHOOTHOLE) << "Removing " << halfedges_to_remove.size()/3 << " triangles.";
  DBG(DBG_SHOOTHOLE) << "Removing " << halfedges_to_remove.size() << " halfedges.";
  DBG(DBG_SHOOTHOLE) << "Found " << removal_boundary_vertices.size() << " boundary vertices:";


  /* We deal with a planar graph, so v-e+f == 2, including the outer face.
   *
   * This also must hold for the area we want to remove.  So we can check
   * if it has enclosed a region without marking it for removal by a counting
   * argument.  It should not have, as shoot_hole_identify_boundary_vertices should
   * have eaten away any enclosed areas.
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
  DEBUG_STMT({
    assert(halfedges_to_remove.size() % 3 == 0);
    DBG(DBG_SHOOTHOLE) << "vertices to remove:";
    for (auto v : vertices_to_remove) {
      DBG(DBG_SHOOTHOLE) << " v" << v->v->idx();
    }
    DBG(DBG_SHOOTHOLE) << "vertices on removal boundary:";
    for (auto v : removal_boundary_vertices) {
      DBG(DBG_SHOOTHOLE) << " v" << v->v->idx() << "; via " << *v;
    }
  });

  {
    int total_num_v  = vertices_to_remove.size() + removal_boundary_vertices.size();
    int total_num_2e = halfedges_to_remove.size() + removal_boundary_vertices.size() + vertices_to_remove_on_ch;
    int total_num_f_accounted_for = halfedges_to_remove.size()/3;
    DBG(DBG_SHOOTHOLE) << " total_num_v: " << total_num_v;
    DBG(DBG_SHOOTHOLE) << " total_num_2e: " << total_num_2e;
    DBG(DBG_SHOOTHOLE) << " total_num_f_accounted_for: " << total_num_f_accounted_for;
    int inner_faces = 1 + total_num_2e/2 - total_num_v - total_num_f_accounted_for;
    DBG(DBG_SHOOTHOLE) << " number of remaining enclosed faces: " << inner_faces;

    if (inner_faces != 0) {
      DBG(DBG_SHOOTHOLE) << "We ate up everything except a small set in the middle.";
    } else {
      assert(total_num_2e % 2 == 0);
      res = true;
    }
  }
 done:
  DBG_FUNC_END(DBG_SHOOTHOLE);
  return res;
}

/** Makes a trianglelist for all the edges in child_edge_pointers.
 *
 * The parentlist we need to figure out neighborhood relations.  This requires
 * additional memory and runtime linear in the size of parentlist.
 */
DECL::SavedState::
SavedState(const FixedVector<Edge> &parentlist, const std::vector<Edge*>& child_edge_pointers, int num_faces_, int num_vertices_, int num_vertices_on_boundary_)
  : num_faces(num_faces_)
  , num_vertices(num_vertices_)
  , num_vertices_on_boundary(num_vertices_on_boundary_)
{
  DBG_FUNC_BEGIN(DBG_STATESAVE);
  std::vector<int> edge_idx_in_child_edge_pointers(parentlist.size(), -1);
  int num_halfedges = child_edge_pointers.size();
  assert(num_halfedges % 3 == 0);
  int num_2edges = num_halfedges + num_vertices_on_boundary;
  DBG(DBG_STATESAVE) << "num_vertices " << num_vertices;
  DBG(DBG_STATESAVE) << "num_vertices_on_boundary " << num_vertices_on_boundary;
  DBG(DBG_STATESAVE) << "num_2edges " << num_2edges;
  DBG(DBG_STATESAVE) << "num_halfedges " << num_halfedges;
  assert(num_2edges % 2 == 0);
  assert(num_vertices - num_2edges/2 + num_halfedges/3 == 1); /* v - e + f == <number of components> */

  DBG(DBG_STATESAVE) << "Saving State of " << num_faces << " faces, backed by " << num_halfedges /3 << " triangles.";
  DBG(DBG_STATESAVE) << "  Half edges:";
  {
    int i=0;
    for (auto e_it : child_edge_pointers) {
      DBG(DBG_STATESAVE) << "   " << *e_it;
      assert(e_it->triangle_to_be_removed);
      assert(e_it->v->vertex_to_be_removed != e_it->v->vertex_on_outer_removal_boundary);

      int idx = (&*e_it - parentlist.data());
      assert(idx == e_it->idx());
      assert(edge_idx_in_child_edge_pointers[idx] == -1);
      edge_idx_in_child_edge_pointers[idx] = i;
      ++i;
    }
  }

  triangle_vertices.resize(num_halfedges);
  constraints.resize(num_halfedges*2);
  buddy.resize(num_halfedges);

  {
    Vertex** tv = triangle_vertices.data();
    int* c = constraints.data();
    int* b = buddy.data();

    DBG(DBG_STATESAVE) << "copying data";
    int i=0;
    for (auto e_it : child_edge_pointers) {
      Edge *e = &*e_it;
      DBG(DBG_STATESAVE) << "  " << *e;
      assert(e == child_edge_pointers[i]);
      assert(child_edge_pointers[3*(i/3) + next_edge_offset[i%3]] == e->next);
      *(tv++) = e->next->v;
      if (e->is_constrained) {
        int idx;
        idx = e->next_constrained - &parentlist[0];
        assert(edge_idx_in_child_edge_pointers.at(idx) != -1);
        *(c++) = edge_idx_in_child_edge_pointers[idx];

        idx = e->prev_constrained - &parentlist[0];
        assert(edge_idx_in_child_edge_pointers.at(idx) != -1);
        *(c++) = edge_idx_in_child_edge_pointers[idx];
      } else {
        *(c++) = -1;
        *(c++) = -1;
      };
      Edge *o = e->opposite;
      int idx = o - &parentlist[0];
      assert(o == NULL || (idx == o->idx()));
      *(b++) = o ? edge_idx_in_child_edge_pointers[ idx ] : -1;
      DBG(DBG_STATESAVE) << "   o: " << o << "; idx: " << (o ? idx : -1) << "; setting b to " << *(b-1);
      ++i;
    }
  }
  DBG_FUNC_END(DBG_STATESAVE);
}

/** Makes a trianglelist for all the edges.
 *
 * This basically dumps the state of an entire edgelist.
 *
 * It is very similar to the other SavedState constructor, but doesn't require additional memory
 * to obtain edge indices for the buddy list.
 */
DECL::SavedState::
SavedState(const FixedVector<Edge> &alledges, int num_faces_, int num_vertices_, int num_vertices_on_boundary_)
  : num_faces(num_faces_)
  , num_vertices(num_vertices_)
  , num_vertices_on_boundary(num_vertices_on_boundary_)
{
  DBG_FUNC_BEGIN(DBG_STATESAVE);
  int num_halfedges = alledges.size();
  assert(num_halfedges % 3 == 0);
  {
    int num_2edges = num_halfedges + num_vertices_on_boundary;
    DBG(DBG_STATESAVE) << "num_vertices " << num_vertices;
    DBG(DBG_STATESAVE) << "num_vertices_on_boundary " << num_vertices_on_boundary;
    DBG(DBG_STATESAVE) << "num_2edges " << num_2edges;
    DBG(DBG_STATESAVE) << "num_halfedges " << num_halfedges;
    assert(num_2edges % 2 == 0);
    assert(num_vertices - num_2edges/2 + num_halfedges/3 == 1); /* v - e + f == <number of components> */
  }

  DBG(DBG_STATESAVE) << "Saving State of " << num_faces << " faces, backed by " << num_halfedges /3 << " triangles.";

  triangle_vertices.resize(num_halfedges);
  constraints.resize(num_halfedges*2);
  buddy.resize(num_halfedges);

  {
    Vertex** tv = triangle_vertices.data();
    int* c = constraints.data();
    int* b = buddy.data();

    DBG(DBG_STATESAVE) << "copying data";
    int i=0;
    for (const auto& e : alledges) {
      DBG(DBG_STATESAVE) << "  " << e;
      assert(&e == &alledges[i]);
      assert(&alledges[3*(i/3) + next_edge_offset[i%3]] == e.next);
      assert(&alledges[3*(i/3) + prev_edge_offset[i%3]] == e.prev);
      *(tv++) = e.next->v;
      if (e.is_constrained) {
        int idx;
        idx = e.next_constrained - &alledges[0];
        assert(idx == e.next_constrained->idx());
        *(c++) = idx;

        idx = e.prev_constrained - &alledges[0];
        assert(idx == e.prev_constrained->idx());
        *(c++) = idx;
      } else {
        *(c++) = -1;
        *(c++) = -1;
      };
      const Edge *o = e.opposite;
      int idx = o ? (o - &alledges[0]) : -1;
      assert(o == NULL || (idx == o->idx()));
      *(b++) = idx;
      DBG(DBG_STATESAVE) << "   o: " << o << "; idx: " << idx;
      if (idx >= 0) {
        DBG(DBG_STATESAVE) << "    that idx is " << alledges[ idx ];
      }
      ++i;
    }
  }
  DBG_FUNC_END(DBG_STATESAVE);
}

/** Inject saved new state into this DECL at the locations of halfedges_to_remove
 */
void
DECL::
shoot_hole_inject_saved_state(const SavedState &state) {
  DBG_FUNC_BEGIN(DBG_STATERESTORE);

  assert_valid();
  assert(halfedges_to_remove.size() == state.triangle_vertices.size());

  Vertex* const * tv = state.triangle_vertices.data();
  const int* c = state.constraints.data();
  const int* b = state.buddy.data();
  std::vector<Edge*> boundary_from_within;

  /* Do the interior */
  int i=0;
  for (auto e_it : halfedges_to_remove) {
    Edge *e = &*e_it;
    #if 0
    if (i % 3 == 0) {
      DBG(DBG_STATERESTORE) << " b " << *(halfedges_to_remove[i  ]);
      DBG(DBG_STATERESTORE) << " b " << *(halfedges_to_remove[i+1]);
      DBG(DBG_STATERESTORE) << " b " << *(halfedges_to_remove[i+2]);

      //DBG(DBG_STATERESTORE) << "   v:" << (*(tv  ))->idx();
      //DBG(DBG_STATERESTORE) << "   v:" << (*(tv+1))->idx();
      //DBG(DBG_STATERESTORE) << "   v:" << (*(tv+2))->idx();

      DBG(DBG_STATERESTORE) << "   c:" << *(c  ) << "; " << *(c+1);
      DBG(DBG_STATERESTORE) << "   c:" << *(c+2) << "; " << *(c+3);
      DBG(DBG_STATERESTORE) << "   c:" << *(c+4) << "; " << *(c+5);
    }
    #endif

    assert(e == halfedges_to_remove[i]);

    assert(halfedges_to_remove[3*(i/3) + next_edge_offset[i%3]] == e->next);
    assert(halfedges_to_remove[3*(i/3) + prev_edge_offset[i%3]] == e->prev);
    e->next->v = *(tv++);
    e->is_constrained = (*c != -1);
    if (e->is_constrained) {
      assert(*(c+1) != -1);
      e->next_constrained = halfedges_to_remove[ *(c++) ];
      e->prev_constrained = halfedges_to_remove[ *(c++) ];
    } else {
      assert(*(c+1) == -1);
      c += 2;
      e->next_constrained = NULL;
      e->prev_constrained = NULL;
    };

    if (*b != -1) {
      e->opposite = halfedges_to_remove[ *(b++) ];
    } else {
      e->opposite = NULL;
      boundary_from_within.emplace_back(e);
      ++b;
    }
    e->triangle_to_be_removed = false;
    ++i;
    #if 0
    if (i % 3 == 0) {
      DBG(DBG_STATERESTORE) << " n " << *(halfedges_to_remove[i-3]);
      DBG(DBG_STATERESTORE) << " n " << *(halfedges_to_remove[i-2]);
      DBG(DBG_STATERESTORE) << " n " << *(halfedges_to_remove[i-1]);
      DBG(DBG_STATERESTORE);
    }
    #endif
  }

  vertices_to_remove.clear();
  halfedges_to_remove.clear();

  /* And now stitch them together along the boundary */
  {
    assert(boundary_from_within.size() >= removal_boundary_vertices.size());

    Edge *prev_vertex = removal_boundary_vertices.back();
    Edge *edge_in_region = NULL;
    DEBUG_STMT({
      DBG(DBG_STATERESTORE) << "Boundary from within:";
      for (auto e : boundary_from_within) {
        DBG(DBG_STATERESTORE) << " " << *e;
      };
      DBG(DBG_STATERESTORE) << "Boundary outside:";
      for (auto e : removal_boundary_vertices) {
        DBG(DBG_STATERESTORE) << " " << *e;
      };
    });

    DBG(DBG_STATERESTORE) << "Stitching:";
    auto e_it = removal_boundary_vertices.begin();
    while ((*e_it)->get_tail() != prev_vertex->v) { /* see comment in for loop below on when an entry is an actual edge */
      DBG(DBG_STATERESTORE) << " Skipping over initial CH vertex " << *(*e_it)->v << "  via " << **e_it;
      prev_vertex = *e_it;
      ++e_it;
      assert(e_it != removal_boundary_vertices.end());
    }

    DBG(DBG_STATERESTORE) << " First boundary edge from outside is " << *(*e_it) << "; Looking for that in the boundary from within";
    Edge *e_inner = NULL;
    {
      Edge *eo = *e_it;
      for (auto e_within_it = boundary_from_within.rbegin(); e_within_it != boundary_from_within.rend(); ++e_within_it) {
        Edge *ei = *e_within_it;
        if (ei->v == eo->get_tail() &&
            eo->v == ei->get_tail()) {
          e_inner = ei;
          DBG(DBG_STATERESTORE) << "   found it: " << *e_inner;
          break;
        }
      }
    }
    assert(e_inner != NULL);

    [[maybe_unused]] bool hit_a_ch_vertex = false;
    for (; e_it != removal_boundary_vertices.end(); ++e_it) {
      Edge *e = *e_it;
      /* This is an edge with the removed and re-injected area if this 
       * edge's previous vertex is the vertex of the previous entry
       * in the list.
       *
       * Otherwise this entry in the list is just a witness for a vertex that
       * is on the CH, and does not represent an edge shared with the region.
       */
      if (e->get_tail() == prev_vertex->v) {
        DBG(DBG_STATERESTORE) << " have edge (" << e->get_tail()->idx() << ", " << e->v->idx() << ") via " << *e;

        if (e_inner->v != e->get_tail()) {
          assert(hit_a_ch_vertex);
          while (e_inner->v != e->get_tail()) {
            DBG(DBG_STATERESTORE) << "   Looking for match in the interior.  This is nowhere close: " << *e_inner;
            while (e_inner->opposite) e_inner = get_next_face_ccw(e_inner);
            e_inner = e_inner->prev_constrained;
            hit_a_ch_vertex = false;
          }
        }
        assert(!hit_a_ch_vertex);

        while (e->v != e_inner->get_tail()) {
          DBG(DBG_STATERESTORE) << "   Looking for match in the interior.  This isn't it: " << *e_inner;
          e_inner = get_next_face_ccw(e_inner);
        }
        DBG(DBG_STATERESTORE) << "   joining up " << *e << " and " << *e_inner;
        e->opposite = e_inner;
        e_inner->opposite = e;
        DBG(DBG_STATERESTORE) << "      are now " << *e << " and " << *e_inner;

        e_inner = e_inner->prev_constrained;
      } else {
        DBG(DBG_STATERESTORE) << " just a CH vertex " << *e->v;
        hit_a_ch_vertex = true;
      }
      prev_vertex = e;
    }
  }

  num_faces += -faces_removed + state.num_faces;

  removal_boundary_vertices.clear();
  vertices_to_remove_on_ch = 0;
  faces_removed = 0;

  assert_valid();
  assert_hole_shooting_reset();

  DBG_FUNC_END(DBG_STATERESTORE);
}

/** Find a simply-connected set of vertices, and make a hole
 */
void
DECL::
shoot_hole(unsigned size, int num_iterations, int max_recurse) {
  //DBG_FUNC_BEGIN(DBG_GENERIC);
  DBG_INDENT_INC();
  DBG(DBG_GENERIC) << "Shoot hole of size " << size << " and try to improve " << num_iterations << " times";

  assert_hole_shooting_reset();

  shoot_hole_select_vertices(size);
  bool res = shoot_hole_identify_affected_elements();
  int affected_vertices = vertices_to_remove.size() + removal_boundary_vertices.size();
  DBG(DBG_GENERIC) << "Affected vertices: " << affected_vertices << " out of " << num_vertices;
  if (!res || affected_vertices >= num_vertices) {
    DBG(DBG_GENERIC) << "Hole shooting failed";
    for (auto e : vertices_to_remove) e->v->vertex_to_be_removed = false;
    for (auto e : halfedges_to_remove) e->triangle_to_be_removed = false;
    for (auto e : removal_boundary_vertices) {
      e->v->vertex_on_outer_removal_boundary = false;
      e->v->vertex_on_removal_boundary = false;
    }

    vertices_to_remove.clear();
    removal_boundary_vertices.clear();
    halfedges_to_remove.clear();
    vertices_to_remove_on_ch = 0;
    faces_removed = 0;

    assert_hole_shooting_reset();
  } else {
    SavedState state(edges, halfedges_to_remove, faces_removed, vertices_to_remove.size() + removal_boundary_vertices.size(), removal_boundary_vertices.size() + vertices_to_remove_on_ch );
    for (auto e : vertices_to_remove) e->v->vertex_to_be_removed = false;
    for (auto e : removal_boundary_vertices) e->v->vertex_on_outer_removal_boundary = false;

    assert_hole_shooting_vertices_clean();

    DECL child(all_vertices, state);

    /* Recurse here */
    SavedState result = child.find_convex_decomposition_many(num_iterations, faces_removed, max_recurse-1);

    if (result.num_faces > 0) {
      /* We have improved */
      int old_num_faces = num_faces;
      shoot_hole_inject_saved_state(result);
      DBG(DBG_GENERIC) << "Moving from " << old_num_faces << " to " << num_faces << " faces";
    } else {
      for (auto e : halfedges_to_remove) e->triangle_to_be_removed = false;
      vertices_to_remove.clear();
      removal_boundary_vertices.clear();
      halfedges_to_remove.clear();
      vertices_to_remove_on_ch = 0;
      faces_removed = 0;
      assert_hole_shooting_reset();
    }
  }
  DBG_INDENT_DEC();
  //DBG_FUNC_END(DBG_GENERIC);
}

void
DECL::
shoot_holes(int max_recurse) {
  DBG_FUNC_BEGIN(DBG_GENERIC);
  assert_hole_shooting_reset();
  int hole_size = num_vertices;
  while (hole_size >= 10) {
    hole_size = int(std::pow(double(hole_size), 2./3));
    int number_of_decompositions_per_hole = hole_size;
    int number_of_hole_punches = num_vertices/hole_size * NUMBER_OF_HOLE_PUNCHES_SCALE;

    DBG(DBG_GENERIC) << "Calling shoot_hole " << number_of_hole_punches << " times with hole_size: " << hole_size << "; number_of_decompositions_per_hole: " << number_of_decompositions_per_hole << "; max_recurse: " << max_recurse;
    for (int i=0; i<number_of_hole_punches; ++i) {
      if ( (number_of_decompositions_per_hole >  100 && i % 100 == 0)
        || (number_of_decompositions_per_hole <= 100 && i % 500 == 0)
            ) {
        DBG(DBG_GENERIC) << "i: " << i << "; current num faces: " << num_faces;
      }
      shoot_hole(hole_size, number_of_decompositions_per_hole, max_recurse);
      assert_hole_shooting_reset();
    }
  };
  DBG_FUNC_END(DBG_GENERIC);
}

void
DECL::
find_convex_decomposition(int max_recurse) {
  DBG_FUNC_BEGIN(DBG_GENERIC);
  unconstrain_all();
  shoot_holes(max_recurse);
  DBG_FUNC_END(DBG_GENERIC);
}

DECL::SavedState
DECL::
find_convex_decomposition_many(int num_iterations, int initial_num_faces_to_beat, int max_recurse) {
  //DBG_FUNC_BEGIN(DBG_GENERIC);
  DBG_INDENT_INC();
  assert_hole_shooting_vertices_clean();
  bool have_solution = false;
  SavedState result;

  int iter = 0;
  int num_faces_to_beat = initial_num_faces_to_beat;
  while (1) {
    unconstrain_all();
    if (max_recurse > 0) {
      shoot_holes(max_recurse);
    }

    int this_num_faces = get_num_faces();

    if (this_num_faces < num_faces_to_beat || num_faces_to_beat == 0) {
      DBG(DBG_GENERIC) << "Iteration " << iter << "/" << num_iterations << ": Improved solution: " << num_faces_to_beat << "->" << this_num_faces;
      have_solution = true;
      num_faces_to_beat = this_num_faces;
      result = SavedState(edges, this_num_faces, num_vertices, num_vertices_on_boundary);
    }
    ++iter;

    if (iter >= num_iterations) break;
    // DBG(DBG_GENERIC) << "Resetting constraints";
    reset_constraints();
  }
  if (have_solution) {
    DBG(DBG_GENERIC) << "Done " << num_iterations << " iterations.  Best now is " << num_faces_to_beat << " from " << initial_num_faces_to_beat;
    assert(num_faces_to_beat == result.num_faces);
  } else {
    DBG(DBG_GENERIC) << "Done " << num_iterations << " iterations.  We failed to improve on the bound of " << num_faces_to_beat << " faces";
  };
  assert_hole_shooting_reset();
  assert(have_solution == (result.num_faces != 0));
  DBG_INDENT_DEC();
  //DBG_FUNC_END(DBG_GENERIC);
  return result;
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
write_obj_segments(bool dump_vertices, std::ostream &o) const {
  if (dump_vertices) {
    for (const auto &v : *all_vertices) {
      o << "v " << v.x << " " << v.y << " 0" << std::endl;
    }
  }
  for (const auto &e : edges) {
    if (!e.is_constrained) continue;
    if (e.opposite && e.opposite < &e) continue;

    int tail_idx = (e.get_tail() - all_vertices->data())+1;
    int head_idx = (e.v - all_vertices->data())+1;

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
  os << " #edges: " << d.edges.size() << std::endl;
  for (unsigned i=0; i<d.edges.size(); ++i) {
    const Edge &e = d.edges[i];

    os << "   edge #" << e
       << std::endl;
  }
  return os;
}
