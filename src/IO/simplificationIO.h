#pragma once

#include <base/cgal_typedefs.h>
#include <CGAL/structure_point_set.h>

#include <exe/simplification.h>

void exportStructured(dirHolder& dir, CGAL::Point_set_with_structure<EPICK>& pss);


void exportConvexHull2d(dirHolder& dir,DT2& as,CGAL::Color& col);
void exportPolygon(dirHolder& dir, Polygon_2& poly, EPICK::Plane_3& plane);
void exportEdges(dirHolder& dir, vector<Segment>& segs);

void exportAlphaShape(dirHolder& dir, vector<Point>& points);

