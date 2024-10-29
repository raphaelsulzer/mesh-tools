#pragma once
#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <util/geometricOperations.h>

#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <boost/foreach.hpp>
#include <random>

typedef CGAL::Side_of_triangle_mesh<SurfaceMesh, EPICK> Point_inside;
typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh> Primitive;
typedef CGAL::AABB_traits_3<EPICK, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

int makeFieldOnDelaunayVertices(dataHolder& data, string mode);
