#pragma once

#include <base/cgal_typedefs.h>
//#ifdef Open3D
//#include "open3d/Open3D.h"
//#endif
#include <IO/fileIO.h>

////////////////////////////////////////////////////////////
///////////////////////////// PCA //////////////////////////
////////////////////////////////////////////////////////////
//#ifdef Open3D
//bool pointInsideBbox(Point& p, const open3d::geometry::OrientedBoundingBox& bb);
//bool pointInsideBbox(Point& p, const open3d::geometry::AxisAlignedBoundingBox& bb);
//#endif

bool pointInsideBbox(Point& p, Bbox& bb);
bool pointInsideBbox(Point& p, SurfaceMesh& sm_obb);
int applyCropping(dirHolder dir, dataHolder& data);
int applyCropping(dataHolder& data);

int clipWithBoundingBox(vector<Point>& input_point, vector<vertex_info>& input_infos, SurfaceMesh& sm_obb);
void getOrientedBoundingBox(vector<Point>& points, SurfaceMesh& obb, runningOptions ro);

double pcaKNN(Point_set& point_set);
// PCA with kNN neighborhood
double pcaKNN(Delaunay& Dt);
// PCA with Delaunay neighborhood
void pcaDt(Delaunay& Dt);



////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////
int applyTransformationToPoint(Point& p, Eigen::Matrix4d& rt);
int applyTransformationMatrix(dataHolder& data);

void scalePointSet(dirHolder dir, dataHolder& data, double scale);
int standardizePointSet(dirHolder dir, dataHolder& data, double scale_factor);
int reStandardizePointSet(dirHolder dir, dataHolder& data);


void gridSimplify(vector<Point>& points, double cell_size);
void gridSimplify(dataHolder& data, double cell_size);

void sampleFun(vector<Point>& points_in, vector<vertex_info>& infos_in,const int percentage);
void sampleFun(vector<Point>& points_in,const int percentage);

void createPointSet(dataHolder& data);
