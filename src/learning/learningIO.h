#ifndef LEARNINGIO_H
#define LEARNINGIO_H

#include <base/cgal_typedefs.h>
#include <learning/learning.h>

using namespace std;

// export functinos
void raySegmentTriangleTest(string out_path, int n_samples);

void segmentTriangleTest(string out_path, int n_samples);
void segmentPlaneIntersection(string out_path, int n_samples);

void segmentTriangleIntersection(string out_path, int n_samples);

void triangleTraining(string out_path, Delaunay& Dt);
void cellTraining(string out_path, Delaunay& Dt, bool cellBased=true);

void tetTrainingHandcrafted(string out_path, Delaunay& Dt);

// export functions
void exportGraph(dirHolder& dir, runningOptions& options, Delaunay& Dt);
void exportFacetRayGraph(dirHolder& dir, runningOptions& options, Delaunay& Dt);


// import functions
int importPrediction(dirHolder dir, dataHolder& data, runningOptions options);


#endif // LEARNINGIO_H
