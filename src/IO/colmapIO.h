#pragma once

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>


////////////////////////////////////////////////
//////////////// WU TENGS CODE /////////////////
////////////////////////////////////////////////
// for colmap fused.ply.vis(binary format)
bool loadVisibilityFile(const string pszVisFile, vector<vector<int>>& vis_info);
bool loadImageFile(const string pszImageFile, map<int, Point>& image_proj_center);
bool getSensorMap(const string pszImageFile, map<int, Point>& image_proj_center);
Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d& qvec);
Eigen::Vector3d ProjectionCenterFromPose(const Eigen::Vector4d& qvec, const Eigen::Vector3d& tvec);
void readColmapFiles(dirHolder dir, dataHolder& data);

void applyTransformation(Point& p, Eigen::Matrix4d& rt);

// old: using colmap classes
void readColmapRecon(dirHolder& dir, dataHolder& data);

