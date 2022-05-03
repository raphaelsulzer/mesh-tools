#pragma once

#include <base/cgal_typedefs.h>

using namespace std;

namespace learning{

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////

int traverseCells(Delaunay& Dt,
                  Cell_handle& current_cell, Ray ray, double firstDist, double maxTraversalDist, int oppositeVertex,
                  int image, bool inside);

void firstCell(Delaunay& Dt,
               Delaunay::Finite_vertices_iterator& vit, Ray ray, int image,
               bool inside);

void rayTracing(Delaunay& Dt, runningOptions ro);

void aggregateScoreAndLabel(Delaunay& Dt);


// end of namespace rayTracing
}
