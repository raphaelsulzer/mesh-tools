#include <base/cgal_typedefs.h>
#include <util/geometricOperations.h>
#include <Eigen/Core>
#include <CGAL/Polyhedron_copy_3.h>

using namespace std;

double dot(const Vector& v1, const Vector& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

double dot(const Point& v1, const Point& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

double dot(const Vector& v1, const Point& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

double dot(const Point& v1, const Vector& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

Vector crossV(const Vector& v1, const Vector& v2){

    return Vector(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Vector crossV(const Point& v1, const Point& v2){

    return Vector(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Point crossP(const Vector& v1, const Vector& v2){

    return Point(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Point crossP(const Point& v1, const Point& v2){

    return Point(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Point neg(const Point& p){

    return Point(-p.x(), -p.y(), -p.z());
}



const double polyhedronVolume(vector<pair<Triangle, Vector>> polyhedron){

    double volume;
    for(auto const& facet : polyhedron){
        volume+= dot(facet.first[0], facet.second);
    }
    return volume/6;

}


double pointPlaneDistance(EPICK::Plane_3 plane, Point point){

    // implemented from: http://mathworld.wolfram.com/Point-PlaneDistance.html

    double dist = (plane.a()*point.x() + plane.b()*point.y() + plane.c()*point.z() + plane.d()) /
            std::sqrt(plane.a()*plane.a() + plane.b()*plane.b() + plane.c()*plane.c());

    return dist;
}



// Implementation of geometry visualized in Figure 9 in P. Labatut, J‐P. Pons,
// and R. Keriven. "Robust and efficient surface reconstruction from range
// data." Computer graphics forum, 2009.
double computeCosFacetCellAngle(const Delaunay& Dt,
                                const Delaunay::Facet& facet) {
  if (Dt.is_infinite(facet.first)) {
    return 1.0;
  }

  const Triangle triangle = Dt.triangle(facet);

  const Vector facet_normal =
      CGAL::cross_product(triangle[1] - triangle[0], triangle[2] - triangle[0]);
  const double facet_normal_length_squared = facet_normal.squared_length();
  if (facet_normal_length_squared == 0.0) {
    return 0.5;
  }

  const Vector co_tangent = CGAL::circumcenter(Dt.tetrahedron(facet.first)) - triangle[0];
  const float co_tangent_length_squared = co_tangent.squared_length();
  if (co_tangent_length_squared == 0.0) {
    return 0.5;
  }

  return (facet_normal * co_tangent) /
         std::sqrt(facet_normal_length_squared * co_tangent_length_squared);
}


double computeFacetArea(const Delaunay& Dt,
                        const Delaunay::Facet& facet){
    if(Dt.is_infinite(facet))
        return 0.0;
    else
        return sqrt(Dt.triangle(facet).squared_area());
};


const Point barycenter(const Cell_handle& ch){
    double x = 0, y = 0, z = 0;
    for(int i = 0; i < 4; i++){
        x+=ch->vertex(i)->point().x();
        y+=ch->vertex(i)->point().y();
        z+=ch->vertex(i)->point().z();
    }
    return Point(x/4.0,y/4.0,z/4.0);
}





bool rayTriangleIntersection(const Point& rayOrigin,
                           const Vector& rayVector,
                           const Triangle& inTriangle,
                           Point& outIntersectionPoint){

    // implemented after: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    using namespace Eigen;
    const double EPSILON = 0.0000001;

    // init "Eigen" vectors
    Vector3d rayO(rayOrigin.x(), rayOrigin.y(), rayOrigin.z());
    Vector3d rayV(rayVector.x(), rayVector.y(), rayVector.z());

    Vector3d vertex0(inTriangle.vertex(0).x(), inTriangle.vertex(0).y(), inTriangle.vertex(0).z());
    Vector3d vertex1(inTriangle.vertex(1).x(), inTriangle.vertex(1).y(), inTriangle.vertex(1).z());
    Vector3d vertex2(inTriangle.vertex(2).x(), inTriangle.vertex(2).y(), inTriangle.vertex(2).z());

    Vector3d edge1, edge2, h, s, q;
    double a;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = rayV.cross(edge2);
    a = edge1.dot(h);
    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.
    double f = 1.0/a;
    s = rayO - vertex0;
    double u = f * s.dot(h);    // barycentric coordinate u
    if (u < 0.0 || u > 1.0)     // if u is not between 0 and 1, intersection point does not lie in triangle
        return false;
    q = s.cross(edge1);
    double v = f * rayV.dot(q); // barycentric coordinate v
    if (v < 0.0 || u + v > 1.0) // if v is not between 0 and 1, intersection point does not lie in triangle
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * edge2.dot(q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayOrigin + rayVector * t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}


