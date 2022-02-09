#ifndef LEARNINGMATH_H
#define LEARNINGMATH_H

#include <base/cgal_typedefs.h>
#include <learning/learning.h>

using namespace std;



void calculateBarycentricCoordinates(const Point& rayOrigin,
                           const Vector& rayVector,
                           const Triangle& inTriangle,
                           Point& outIntersectionPoint,
                                     int& intersection);

double pointDistance(const Tree& tree, const Point &query);

tetFeatures calcTetFeatures(const Tetrahedron& tet);

#endif // LEARNINGMATH_H
