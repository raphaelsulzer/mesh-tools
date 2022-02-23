#ifndef FILEIO_H
#define FILEIO_H

#include <base/cgal_typedefs.h>
#include <exe/exeOptions.h>
#ifdef RECONBENCH
#include <modeling/ImplicitFunction.h>
#endif

//#include "open3d/Open3D.h"
//#include "open3d/geometry/TetraMesh.h"

struct dirHolder{

    string path;
    string mvs_file;
    string scene;

    string gt_poly_file;
    string gt_scan_file;
    string prediction_file;

    string read_file;
    string write_file;
    string read_file_type;

    string transformation_file;
    string crop_file;

    string suffix = "";
    string rw_string;
    string ol_string = "";

};
struct dataHolder{

//    MVS::Scene omvs_scene;

    vector<Point> points;
    vector<vertex_info> infos;
    bool has_color, has_normal, has_gt_normal, has_sensor;
    vector<array<size_t,3>> facets;
    map<int, Point> sensor_map;

    // xtensor
    vector<double> xpoints;
    vector<double> xnormals;
    vector<double> xgtnormals;
    vector<double> xsensor_positions;
    vector<double> xverts; // this is for saving the points in the order they appear in the 3DT
    vector<int> xfacets;
    vector<int> xnfacets; // cell neighbors of facets
    vector<int> xtets;


    // learning
//    vector<Point> sampled_sensors;
//    vector<Point> sampled_points;

    vector<Point> gt_points;
    vector<vertex_info> gt_infos;

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

int importImplicit(dirHolder dir, dataHolder& data);

int importTransformationMatrix(dirHolder dir, dataHolder& data);

int importPLYPoints(dirHolder dir, dataHolder& data);
int importPLYMesh(dirHolder& dir, SurfaceMesh& import_mesh);
int importPLYMeshWithSensorStructure(dirHolder dir, dataHolder& data);

int importMesh(dirHolder& dir, dataHolder& data);

int importOFFMesh(dirHolder& dir, Polyhedron& import_poly);
int importOFFMesh(const string path, Polyhedron& import_poly);
int importOFFMesh(dirHolder& dir, SurfaceMesh& import_mesh);
int importNPZ(dirHolder& dir, dataHolder& data);

int importOMVSScene(dirHolder dir, dataHolder& data);

/////////////////////////////////////////////////////////////////////
/////////////////////////////// OUTPUT //////////////////////////////
/////////////////////////////////////////////////////////////////////
void printPLYHeader(std::fstream& fo, exportOptions& eo,
                    int nv, int nf,
                    int precision);


void exportCameraCenter(dirHolder& dir, dataHolder& data);

void exportOFF(dirHolder& dir, Polyhedron& poly);
void exportOFF(dirHolder& dir, SurfaceMesh& out_mesh);

int exportNPZ(dirHolder& dir, dataHolder& data);
int export3DT(const dirHolder dir, dataHolder data);
void exportPLY_manual(dirHolder& dir, SurfaceMesh& out_mesh);
void exportPLY_bin(dirHolder& dir, SurfaceMesh& out_mesh);
int exportPLY(dirHolder& dir, Point_set& points);
int exportPLY(dirHolder& dir, vector<Point>& points);
int exportPLY(dirHolder& dir, SurfaceMesh& out_mesh);
void exportPLY(dirHolder&, vector<Point>&, vector<vertex_info>&, exportOptions&);

void fixSensorCenter();

///// custom /////

void exportMITSurface(dirHolder& dir, Delaunay& Dt, bool mean, double iso_value);
void exportIsoValues(dirHolder& dir, Delaunay& Dt, vector<double>& values);
void exportCellCenter(dirHolder& dir, const Delaunay& Dt);
void exportCellScore(dirHolder& dir, const Delaunay& Dt);


void exportColoredFacets(dirHolder& dir, const Delaunay& Dt,
                            bool optimized);

void exportConvexHull(dirHolder& dir, const Delaunay& Dt);

void exportInterface(dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& eo);


#include <CGAL/structure_point_set.h>
void exportStructured(dirHolder& dir, PSS& pss);


#endif

