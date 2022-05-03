#pragma once
#include <base/cgal_typedefs.h>

using namespace std;

namespace learning{

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////

int nextGtCells(Delaunay& Dt,
                  Cell_handle& current_cell, Ray ray, Vertex_handle vit, int oppositeVertex,
                  double sigma, int image, bool inside);

void firstGtCell(Delaunay& Dt,
               Delaunay::Finite_vertices_iterator& vit, Ray ray, int image,
               bool inside);

void countGtRays(dataHolder& data);

// end of namespace rayTracing
}
