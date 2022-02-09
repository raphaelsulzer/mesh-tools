#ifndef MESHPROCESSING_H
#define MESHPROCESSING_H


#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>


struct meshProcessingOptions{

    int try_to_close = 0;
    int try_to_make_manifold = 0;
    int number_of_components_to_keep = 0; // 0 meaning all
    int factor_for_removing_large_faces = 0; // 0 meaning remove no faces

    // edge collapse
    double keep_edges;

    // variational shape approximation
    int max_number_of_proxies;

    // ransac
    double ransac_epsilon;
    double ransac_normal;
    int ransac_regularization;

    //structure
    double structure_epsilon;

    // snapping
    double snap_epsilon;

    // simplify
    int simplify;



};

////////////////////////////////////////////////////////////
////////////////////// Triangulation ///////////////////////
////////////////////////////////////////////////////////////
void makeDelaunayWithInfo(dataHolder& data);
void makeAdaptiveDelaunayWithInfo(dataHolder& data, double epsilon=0.1);

void indexDelaunay(dataHolder& data, runningOptions options);


bool scanObjectClosed(dataHolder& data, runningOptions& options);
bool scanObjectOpen(dataHolder& data, runningOptions& options);


////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////

void preprocessSensorMesh(std::vector<Point>& points,
                          std::vector<vertex_info>& infos,
                          std::vector<std::vector<int>>& polys);




void edgeCollapse(dataHolder& data, meshProcessingOptions mo);

void holeFilling(SurfaceMesh& poly);

void variationalShapeApproximation(dataHolder& data, meshProcessingOptions mo);



////////////////////////////////////////////////////////////
/////////////////// postprocessing functions ////////////////
////////////////////////////////////////////////////////////

void removeLargeFaces(dataHolder& data, int factor);

void createSurfaceMesh(dataHolder& data,
                       meshProcessingOptions& options);


////////////////////////////////////////////////////////////
/////////////////// info functions ////////////////
////////////////////////////////////////////////////////////
void calculateMeanTriangleArea(dataHolder& data);
void calculateMeanCellVolume(dataHolder& data);


void createVertexSet(dataHolder& data);

Point getBoundingBoxCentroid(std::array<Point, 8>& obb_points, double z_offset);
void getOrientedBoundingBox(Polyhedron& poly, std::array<Point, 8>& obb_points, SurfaceMesh& obb);


#endif // MESHPROCESSING_H
