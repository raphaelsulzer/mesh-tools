#pragma once

#include <base/cgal_typedefs.h>

using namespace std;

namespace processing{

bool rayTriangleIntersection(const Point& rayOrigin,
                           const Vector& rayVector,
                           const Triangle& inTriangle,
                           Point& outIntersectionPoint);

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
pair<double, double> wasureScore(double dist2, Vertex_handle vh, bool inside);
pair<double, double> labatutScore(double dist2, Vertex_handle vh, bool inside);


int traverseCells(Delaunay& Dt,
                  Cell_handle& current_cell, Ray ray, Vertex_handle vit, int oppositeVertex,
                  double sigma, int image, bool inside,
                  string scoreType);

void firstCell(Delaunay& Dt,
               Delaunay::Finite_vertices_iterator& vit, Ray ray, int image,
               bool inside,
               string scoreType);

void rayTracing(Delaunay& Dt, runningOptions ro);


// end of namespace rayTracing
}



