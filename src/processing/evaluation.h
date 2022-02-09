#ifndef EVALUATION_H
#define EVALUATION_H


#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

int sampleMesh(dataHolder& data, exportOptions exportO);
double calcIOU(Polyhedron& gt, SurfaceMesh& reconm, int test_points, double& iou);
double calcIOU(dataHolder& data, int test_points, double& iou);
double calcMeanEdgeLength(dataHolder& data);

void checkMeshQuality(std::vector<Point>& ori_points, Polyhedron& surface_mesh);
int printMeshEvaluation(dirHolder dir, dataHolder& data, double iou = -1);

#endif // EVALUATION_H
