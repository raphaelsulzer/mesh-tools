#pragma once



#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <exe/labatut.h>
#include <processing/rayTracingTet.h>
#include <learning/learningRayTracing.h>
#include <util/geometricOperations.h>

#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <boost/foreach.hpp>
#include <random>


//typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
//typedef CGAL::AABB_traits<K, Primitive> Traits;
//typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::Side_of_triangle_mesh<SurfaceMesh, EPICK> Point_inside;
typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

using namespace std;


struct tetFeatures{

    double vol = 0;
    double longest_edge = 0;
    double shortest_edge = 0;
    double radius = 0;
};



int labelObjectWithClosedGroundTruth(dirHolder dir,dataHolder& data, int sampling_points);
int labelObjectWithOpenGroundTruth(dirHolder dir,dataHolder& data, int sampling_points);
int labelObjectWithImplicit(dataHolder& data, int sampling_points);

//int checkLabelWithLidarScan(dataHolder& data);


