#pragma once

#include <base/cgal_typedefs.h>
//#ifdef Open3D
//#include <open3d/Open3D.h>
//#endif
#include <IO/fileIO.h>

////////////////////////////////////////////////////////////
///////////////////////////// PCA //////////////////////////
////////////////////////////////////////////////////////////


// Concurrency
int estimateNormals(dataHolder& data, runningOptions ro);
int orientNormalsWithSecondCloud(dataHolder& data1, dataHolder& data2);

bool sortSensors(const std::pair<double, Point> &a, const std::pair<double, Point> &b);
void orderSensors(dataHolder& data);
