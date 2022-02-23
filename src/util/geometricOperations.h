#pragma once

#include <base/cgal_typedefs.h>
#include <Eigen/Core>

using namespace std;

double dot(const Vector& v1, const Vector& v2);

double dot(const Point& v1, const Point& v2);

double dot(const Vector& v1, const Point& v2);

double dot(const Point& v1, const Vector& v2);

Vector crossV(const Vector& v1, const Vector& v2);

Vector crossV(const Point& v1, const Point& v2);

Point crossP(const Vector& v1, const Vector& v2);

Point crossP(const Point& v1, const Point& v2);

Point neg(const Point& p);

Polyhedron_Exact inexact2exactPolyhedron(Polyhedron in);
Polyhedron exact2inexactPolyhedron(Polyhedron_Exact ex);

const double polyhedronVolume(vector<pair<Triangle, Vector>> polyhedron);
double pointPlaneDistance(EPICK::Plane_3 plane, Point point);

double computeCosFacetCellAngle(const Delaunay& Dt,
                                const Delaunay::Facet& facet);
double computeFacetArea(const Delaunay& Dt,
                                const Delaunay::Facet& facet);

bool rayTriangleIntersection(const Point& rayOrigin,
                           const Vector& rayVector,
                           const Triangle& inTriangle,
                           Point& outIntersectionPoint);


template <typename T>
T Percentile(const std::vector<T>& elems, const double p) {
//  CHECK(!elems.empty());
//  CHECK_GE(p, 0);
//  CHECK_LE(p, 100);

  const int idx = static_cast<int>(std::round(p / 100 * (elems.size() - 1)));
  const size_t percentile_idx =
      std::max(0, std::min(static_cast<int>(elems.size() - 1), idx));

  std::vector<T> ordered_elems = elems;
  std::nth_element(ordered_elems.begin(),
                   ordered_elems.begin() + percentile_idx, ordered_elems.end());

  return ordered_elems.at(percentile_idx);
}


const Point barycenter(const Cell_handle& ch);



