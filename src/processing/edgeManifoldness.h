//#ifndef MANIFOLDNESS_H
//#define MANIFOLDNESS_H

#pragma once

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <processing/rayTracingTet.h>
#include <processing/tetIntersection.h>

typedef std::pair<std::pair<Cell_handle, int>, double> Facet_score;

// get non manifold edges
int isManifoldEdge(Delaunay& Dt, Delaunay::Finite_edges_iterator& e);

// get non manifold edges
int isManifoldEdge(Delaunay& Dt, Edge e);

int isManifoldClique(Delaunay& Dt, Delaunay::Finite_edges_iterator& e);


int getCellLabel(Cell_handle& c, bool manifold);


double nonManifoldCliqueEnergy(const Delaunay& Dt, Delaunay::Finite_edges_iterator& e, double reg_weight);

#include <boost/dynamic_bitset.hpp>
typedef std::pair<double, boost::dynamic_bitset<>> Combination_score;
int makeCombinations(int n, std::vector<Combination_score>& c);

bool sortCombinations(const Combination_score &a,
              const Combination_score &b);

void fixNonManifoldEdges(Delaunay& Dt, double regularization_weight);




////////////////////////////////////
//////// Point manifoldness ////////
////////////////////////////////////


// TODO: could also implement this with kind of a recursive function that always looks in neighbourhood
// so start with a seed cell, go to next cell and see if it has a different label and if it is connected to the previous cell
// ... if I have more than two regions, point is not manifold

// get non manifold points
int isManifoldPoint(Delaunay& Dt, Vertex_handle p);

// get non manifold points
int isManifoldPoint(Delaunay& Dt, Delaunay::Finite_vertices_iterator& p, std::vector<Cell_handle>& incident_cells);

// get non manifold edges
int isPointManifoldClique(Delaunay& Dt, Delaunay::Finite_vertices_iterator& p);

int isPEManifoldClique(Delaunay& Dt, Delaunay::Finite_edges_iterator& e);

double nonManifoldCliqueEnergy(const Delaunay& Dt, std::vector<Cell_handle> all_cells, double reg_weight);

void getNonManifoldPoints(Delaunay& Dt, double regularization_weight, std::string path);


void fixNonManifoldPoints(Delaunay& Dt, double regularization_weight);



//#endif // MANIFOLDNESS_H
