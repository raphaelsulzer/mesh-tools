#ifndef IO_H
#define IO_H

#include <base/cgal_typedefs.h>
#include <exe/exeOptions.h>
#ifdef RECONBENCH
#include <modeling/ImplicitFunction.h>
#endif

//#include "open3d/Open3D.h"
//#include "open3d/geometry/TetraMesh.h"

struct dirHolder{

    // TODO: change all of these strings to boost::filesystem::path
    // should not be thaat much work, simply need to replace all the path concatenation from
    // string addition to boost path concatenation

    string path;

    string read_file;
    string read_file_type;
    string write_file;

    string gt_poly_file;
    string gt_scan_file;
    string prediction_file;
    string occ_file;

    string transformation_file;
    string crop_file;

    string suffix = "";

};
#include <xtensor/xbuilder.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor-io/xnpz.hpp>
struct dataHolder{

//    MVS::Scene omvs_scene;

    vector<Point> points;
    vector<vertex_info> infos;
    bool has_color = false, has_normal = false, has_gt_normal = false, has_sensor = false;
    vector<array<size_t,3>> facets;
    map<int, Point> sensor_map;

    // xtensor
    vector<double> xpoints;
    vector<double> xnormals;
    vector<double> xgtnormals;
    vector<double> xsensor_positions;
    vector<double> xverts; // this is for saving the points in the order they appear in the 3DT
    vector<int> xvertsi; // vert is infinite?
    vector<int> xfacets;
    vector<int> xfacetsi; // facet is infinite?
    vector<int> xnfacets; // cell neighbors of facets
    vector<int> xtets;
    vector<int> xtetsi; // tet is infinite?
    // eval
//    vector<int> xpoint2tet;
//    vector<bool> xocc; // xpoints.size() != xocc.size(); because first is the scan, second the eval points


    // learning
//    vector<Point> sampled_sensors;
//    vector<Point> sampled_points;

    vector<Point> gt_points;
    vector<vertex_info> gt_infos;

//    vector<double> xgt_points;
//    vector<bool> gt_occupancies;
    xt::xarray<double> xgt_points;
    xt::xarray<uint8_t> xgt_occupancies;
    vector<int> xgt_point2tet;

    Point translation_vector;
    double scale_factor = 0.0;

    vector<PointVectorPair> pVP;

    vector<Point> remaining_points;
    vector<array<size_t,3>> remaining_facets;


    Delaunay Dt;

    double mean_triangle_area = 0.0;
    double mean_cell_volume = 0.0;

    Polyhedron poly;
    Polyhedron gt_poly;
    Point gt_centroid;
    Polyhedron reconstructed_poly;

    SurfaceMesh smesh;
    SurfaceMesh vsa_mesh;

    Eigen::Matrix4d transformation_matrix;

    // could be nice but not supported by every algorithm
    Point_set point_set;

    // ransac
    Efficient_ransac efficient_ransac;
    vector<PointVectorPair> ransac_points;
    vector<EPICK::Plane_3> ransac_planes;

    // structuting
    vector<Point> structured_points;
    vector<Segment> pss_segments;
//    PSS pss;

    // bounding box
    array<Point, 8> bb_array;
    Polyhedron bb_poly;
    SurfaceMesh bb_surface_mesh;

    // alpha shape
    vector<vector<Point>> ass_points;

    // polygons
    vector<Polygon_2> alpha_polygons;


    #ifdef RECONBENCH
    ImplicitFunction* implicit;
    #endif

};
struct exportOptions{

    // points
    bool points = false;
    bool normals = false;
    bool sensor_position = false;
    bool sensor_vec = false;
    bool color = false;
    bool cam_index = false;
    bool score = false;

    // facets
    bool facetColor = false;

    // scan
    bool scan = false;
    bool cameras = false;

    // reconstruction
    bool convexHull = false;
    bool sampling = false;
    bool cellScore = false;
    bool coloredFacets = false;
    bool mesh = false;
    bool interface = true;
    bool isosurface = false;

    bool toply = false;
    bool tonpz = false;

    string sampling_method;
    double sampling_method_option; // specify grid spacing when mesh is sampled for evaluation
};


//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////
void concatenateData(dataHolder& data1, dataHolder& data2,int copyInfo = 1);

int toXTensor(dataHolder& data);

int importImplicit(const dirHolder dir, dataHolder& data);

int importTransformationMatrix(const dirHolder dir, dataHolder& data);

int importPLYPoints(const dirHolder dir, dataHolder& data);
int importPLYMesh(const dirHolder& dir, SurfaceMesh& import_mesh);
int importPLYMeshWithSensorStructure(const dirHolder dir, dataHolder& data);

int importMesh(const dirHolder& dir, dataHolder& data);

int importOFFMesh(const dirHolder& dir, Polyhedron& import_poly);
int importOFFMesh(const string path, Polyhedron& import_poly);
int importOFFMesh(const dirHolder& dir, SurfaceMesh& import_mesh);
int importNPZPoints(const dirHolder& dir, dataHolder& data);

int importOMVSScene(const dirHolder dir, dataHolder& data);

/////////////////////////////////////////////////////////////////////
/////////////////////////////// OUTPUT //////////////////////////////
/////////////////////////////////////////////////////////////////////
void printPLYHeader(std::fstream& fo, exportOptions& eo,
                    int nv, int nf,
                    int precision);


void exportCameraCenter(const dirHolder& dir, dataHolder& data);

void exportOFF(const dirHolder& dir, Polyhedron& poly);
void exportOFF(const dirHolder& dir, SurfaceMesh& out_mesh);

int exportNPZ(const dirHolder& dir, dataHolder& data);
int export3DT(const dirHolder dir, dataHolder& data);
void exportPLY_manual(const dirHolder& dir, SurfaceMesh& out_mesh);
void exportPLY_bin(const dirHolder& dir, SurfaceMesh& out_mesh);
int exportPLY(const dirHolder& dir, Point_set& points);
int exportPLY(const dirHolder& dir, vector<Point>& points);
int exportPLY(const dirHolder& dir, SurfaceMesh& out_mesh);
void exportPLY(const dirHolder&, vector<Point>&, vector<vertex_info>&, exportOptions&);

void fixSensorCenter();

///// custom /////

void exportMITSurface(const dirHolder& dir, Delaunay& Dt, bool mean, double iso_value);
void exportIsoValues(const dirHolder& dir, Delaunay& Dt, vector<double>& values);
void exportCellCenter(const dirHolder& dir, const Delaunay& Dt);
void exportCellScore(const dirHolder& dir, const Delaunay& Dt);


void exportColoredFacets(const dirHolder& dir, const Delaunay& Dt,
                            bool optimized);

void exportConvexHull(const dirHolder& dir, const Delaunay& Dt);

void exportInterface(const dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& eo);


#include <CGAL/structure_point_set.h>
void exportStructured(const dirHolder& dir, PSS& pss);


#endif

