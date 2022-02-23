#ifndef HELPER_H
#define HELPER_H


#include <base/cgal_typedefs.h>

using namespace std;



#ifdef RECONBENCH
const Vector3RB point2Vector3RB(const Point p);
#endif

const Eigen::Vector3d cgal2Eigen(const Point p);
const Point eigen2Cgal(const Eigen::Vector3f p);

vector<Eigen::Vector3d> cgal2Eigen(vector<Point> in);
vector<Point> eigen2Cgal(const vector<Eigen::Vector3d> in);

Polyhedron_Exact inexact2exactPolyhedron(Polyhedron in);

Polyhedron exact2inexactPolyhedron(Polyhedron_Exact ex);

size_t splitString(const string& txt, vector<string>& strs, char ch);

string double2string(double regularization_weight);

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g);

template<typename Iter>
Iter select_randomly(Iter start, Iter end);

string getCurrentPath();

string ExePath();

string whichOS();

#endif
