#pragma once

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

// for colmap fused.ply.vis(binary format)
bool loadVisibilityFile(const char* pszVisFile, std::vector<std::vector<int>>& vis_info);

// because colmap save qvec and tvec
// we need cvec
Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d& qvec);

Eigen::Vector3d ProjectionCenterFromPose(const Eigen::Vector4d& qvec,
                                         const Eigen::Vector3d& tvec);

bool loadImageFile(const char* pszImageFile, std::map<int, Point>& image_proj_center);

void getSensorMap(const std::string path, std::map<int, Point>& sensorMap);

bool sortSensors(const std::pair<double, Point> &a,
              const std::pair<double, Point> &b);

void applyTransformation(Point& p, Eigen::Matrix4d& rt);

void readColmapRecon(dirHolder& dir, dataHolder& data);

void readColmapPoints(const dirHolder& dir, dataHolder& data);
