#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H



#include <CGAL/property_map.h>
#include <CGAL/intersections.h>
#include <CGAL/intersections_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Line_3.h>




using namespace std;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

/////////////////////////////////////////////////
//////////////////// Kernel /////////////////////
/////////////////////////////////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel         EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel           EPECK;



/////////////////////////////////////////////////
/////////////////// Point (set) /////////////////
/////////////////////////////////////////////////
#include <CGAL/Point_set_3.h>
typedef EPICK::Vector_3                                            Vector;
typedef EPICK::Point_3                                             Point;

// for estimate normals of a point set
typedef pair<Point, Vector> PointVectorPair; // the only thing where this is still used is in the estimateNormals() function
// not used, but could be nice
typedef CGAL::Point_set_3<Point> Point_set;



/////////////////////////////////////////////////
/////////////// Other geometry //////////////////
/////////////////////////////////////////////////

typedef EPICK::Iso_cuboid_3                                        Iso_cuboid;
typedef CGAL::Bbox_3                                               Bbox;
typedef EPICK::Ray_3                                               Ray;
typedef EPICK::Triangle_3                                          Triangle;
typedef EPICK::Intersect_3                                         Intersect;
typedef EPICK::Segment_3                                           Segment;
typedef EPICK::Line_3                                              Line;
typedef EPICK::Sphere_3                                            Sphere;
typedef EPICK::Tetrahedron_3                                       Tetrahedron;
typedef CGAL::Cartesian_converter<EPICK,EPECK>                     IK_to_EK;
typedef CGAL::Cartesian_converter<EPECK,EPICK>                     EK_to_IK;

#include <CGAL/Polygon_2.h>
typedef CGAL::Polygon_2<EPICK> Polygon_2;


/////////////////////////////////////////////////
/////////// 3D Delaunay triangulation ///////////
/////////////////////////////////////////////////
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
struct ray_intersection_info{

    int cell_idx;
    double first_intersection;
    double second_intersection;
    int outside;
};


struct vertex_info{
    // point/vertex information
    int global_idx = -1;
    int finite_idx = -1;
    double sigma = 0.01; // noise of the point;
    int alpha = 1; // how many points this vertex represents, will be increased by adaptive Delaunay triangulation
    CGAL::Color color;
    Vector normal;
    Vector gt_normal;
    // sensor information
    Vector sensor_vec;
    vector<Point> sensor_positions;
    int sensor_idx; // the index of the sensor which recorded this point.
    // TODO: for LiDAR this should be running same as the point idx; and for MVS this is simply the camera idx.
    // thus at some point maybe this index should be incooperated in the sensor_positions vector.
//    vector<int> sensor_tet; // see fileOP.cpp for explanation
    // learning
    vector<ray_intersection_info> ray; // the int here is the index of one of the cells that this ray has passed
    int outlier;

    // iso surface extraction
    double field_value;

};
typedef CGAL::Triangulation_vertex_base_with_info_3<vertex_info, EPICK>    VB;

struct cell_info{
    // indices
    int global_idx = -1; // this is an ordered, contiguous zero-based index for the All_cells_iterator;
    int finite_idx = -1; // this is an ordered, contiguous zero-based index for the Finite_cells_iterator;
    // graph cut, label
    int gc_label = 2;   // 0 meaning inside; 1 meaning outside !!, see exportColoredFacetsPLY function in fileIO
//    double gt_inside = 0.0;
//    double gt_outside = 0.0;


    // graph_cut, cell based
    double inside_score = 0.0; // used for occupancy and ray scoring, depending on mode
    double outside_score = 0.0;
//    int inside_count = 0;
//    int outside_count = 0;
//    double outside_dist = 0.0;
//    double inside_dist = 0.0;

    // graph_cut, facet based
    vector<double> facet_weights = vector<double>(4);


    // learning, cell based features
    vector<double> cb_vertex_inside;
    vector<double> cb_vertex_outside;
    vector<double> cb_vertex_last;
    vector<double> cb_facet_inside_first;
    vector<double> cb_facet_inside_second;
    vector<double> cb_facet_outside_first;
    vector<double> cb_facet_outside_second;
    vector<double> cb_facet_last_first;
    vector<double> cb_facet_last_second;


    // learning, facet based features
    vector<vector<double>> fb_vertex_inside = vector<vector<double>>(4);
    vector<vector<double>> fb_vertex_outside = vector<vector<double>>(4);
    vector<vector<double>> fb_vertex_last = vector<vector<double>>(4);
    vector<vector<double>> fb_facet_inside = vector<vector<double>>(4);
    vector<vector<double>> fb_facet_outside = vector<vector<double>>(4);
    vector<vector<double>> fb_facet_last = vector<vector<double>>(4);



//    // for traversing photogrammetry Dt with Lidar rays, i.e. count how many ground truth rays traverse the tets build with MVS data
//    int gt_inside_first = 0;
//    int gt_outside_first = 0;
//    int use_for_learning = 1;

    // manifoldness
    int manifold_label = 2; // 2 meaning no fixed label exists
    vector<int> nm_edge_ids;


};
typedef CGAL::Triangulation_cell_base_with_info_3<cell_info, EPICK>        CB;         // cell base

// Delaunay triangulation data structure
typedef CGAL::Triangulation_data_structure_3<VB, CB>                Tds;        // triangulation data structure
typedef CGAL::Delaunay_triangulation_3<EPICK, Tds>                  Delaunay;   // delaunay triangulation based on triangulation data structure
typedef Delaunay::Edge                                              Edge;
typedef Delaunay::Facet                                             Facet;
typedef Delaunay::Cell_handle                                       Cell_handle;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
typedef Delaunay::Simplex                                           Simplex;




////////////////////////////////////////////////////////////////
//////////////////////////// Polyhedron ////////////////////////
////////////////////////////////////////////////////////////////
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
namespace PMP = CGAL::Polygon_mesh_processing;


typedef CGAL::Polyhedron_3<EPICK, CGAL::Polyhedron_items_with_id_3>     Polyhedron;
typedef CGAL::Polyhedron_3<EPECK, CGAL::Polyhedron_items_with_id_3>     Polyhedron_Exact;


////////////////////////////////////////////////////////////////
////////////////////////// Surface Mesh ////////////////////////
////////////////////////////////////////////////////////////////
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>
#include <CGAL/Unique_hash_map.h>
typedef CGAL::Surface_mesh<EPECK::Point_3> SurfaceMeshExact;
typedef CGAL::Surface_mesh<Point> SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor             vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor               edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_iterator                 edge_iterator;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor           halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor               face_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;


////////////////////////////////////////////////////////////////
///////////////////////////// Ransac ///////////////////////////
////////////////////////////////////////////////////////////////
#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Regularization/regularize_planes.h>

// Efficient RANSAC types
typedef CGAL::Shape_detection::Efficient_RANSAC_traits
<EPICK, Point_set, Point_set::Point_map, Point_set::Vector_map>     RansacTraits;

typedef  CGAL::Shape_detection::Efficient_RANSAC<RansacTraits>      Efficient_ransac;
typedef  CGAL::Shape_detection::Plane<RansacTraits>                 RansacPlane;
typedef Efficient_ransac::Plane_range                               RansacPlaneRange;

#include <CGAL/structure_point_set.h>
typedef CGAL::Point_set_with_structure<EPICK> PSS;

#endif


