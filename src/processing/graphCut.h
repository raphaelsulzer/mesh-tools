#pragma once



// for GCoptimization
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "GCoptimization.h"

#include <base/cgal_typedefs.h>
#include "exe/exeOptions.h"
#include "IO/fileIO.h"

#include <util/geometricOperations.h>

typedef double gtype;


////////////////////////////////////////////////////////////
////////////////////// Feature scaling /////////////////////
////////////////////////////////////////////////////////////
void standardizeScores(Delaunay& Dt);

void normalizeScores(Delaunay& Dt);

void logScore(Delaunay& Dt);

void softmax(Delaunay& Dt);

////////////////////////////////////////////////////////////
//////////////////////// Optimization //////////////////////
////////////////////////////////////////////////////////////

double smoothCost(int s1, int s2, int l1, int l2, void *data);

//// in this version, set data and smoothness terms using arrays
//// grid neighborhood is set up "manually". Uses spatially varying terms. Namely
//// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with
//// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1
//std::pair<std::map<Cell_handle, int>, std::vector<int>> GeneralGraph_DArraySArraySpatVarying(std::pair<Delaunay&, Cell_map&> dt_cells, std::map<Cell_handle, int>& cell_indexMap, std::vector<int> result, int num_iterations)
void graphCutTet(dataHolder& data, runningOptions options);
void graphCutFacet(dataHolder& data, runningOptions options);






