#pragma once

#include <exe/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

#include <util/helper.h>
#include <util/geometricOperations.h>
#include <processing/meshProcessing.h>
#include <processing/evaluation.h>

#include "boost/tuple/tuple.hpp"
#include <CGAL/jet_estimate_normals.h>


#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Regularization/regularize_planes.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>


// from the MeshPolygonization code
//// GRAPH //
//struct GraphVertex {
//    unsigned int segment;
//};

//typedef boost::adjacency_list<boost::setS,
//                              boost::vecS,
//                              boost::undirectedS,
//                              GraphVertex> Graph;

//typedef Graph::vertex_descriptor Graph_vertex;
//typedef Graph::vertex_iterator   Graph_vertex_iterator;
//typedef Graph::edge_descriptor   Graph_edge;
//typedef Graph::edge_iterator     Graph_edge_iterator;
// GRAPH //


#include "boost/tuple/tuple.hpp"
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/structure_point_set.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Optimal_bounding_box/oriented_bounding_box.h>

enum Vertex_status { PLANE, EDGE, CORNER, FREE };

struct index_dist_map{
    pair<int, double> corner;
    pair<int, double> edge;
    pair<int, double> plane;
    double free;
};


//void ransac(dirHolder& dir, dataHolder& data,  meshProcessingOptions options);

//typedef CGAL::Point_set_with_structure<EPICK> PSS;
//PSS structure(dirHolder& dir, dataHolder& data,
//               meshProcessingOptions options);
//void snapPoint_set(dirHolder& dir, dataHolder& data,
//          PSS& pss,
//          meshProcessingOptions options);
//void snapMesh(dirHolder& dir, dataHolder& data,
//          PSS& pss,
//          meshProcessingOptions options);
//void straighten(dirHolder& dir, dataHolder& data, meshProcessingOptions options);


#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
typedef CGAL::Triangulation_vertex_base_with_info_2<Point,EPICK> Vb2;
typedef CGAL::Triangulation_face_base_2<EPICK> Fb2;
typedef CGAL::Triangulation_data_structure_2<Vb2,Fb2> Tds2;
typedef CGAL::Triangulation_2<EPICK,Tds2> DT2;


bool sort_by_area(pair<Polygon_2, double>& a,
                    pair<Polygon_2, double>& b);

bool sort_by_dist(pair<int, double>& a,
                    pair<int, double>& b);

