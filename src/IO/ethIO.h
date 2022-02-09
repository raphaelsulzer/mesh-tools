#ifndef ETHIO_H
#define ETHIO_H

#include <base/cgal_typedefs.h>
#include <exe/exeOptions.h>

#ifdef PCL
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloudPtr;
#endif

////////////////////////////////////////////////////////////
/////////////////////// FILE I/O ///////////////////////////
////////////////////////////////////////////////////////////

int importETH3D(dirHolder dir, dataHolder& data);


#endif
