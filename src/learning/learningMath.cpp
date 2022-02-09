#include <base/cgal_typedefs.h>
#include <learning/learningMath.h>
#include <learning/learning.h>
#include <util/geometricOperations.h>

using namespace std;

void calculateBarycentricCoordinates(const Point& rayOrigin,
                           const Vector& rayVector,
                           const Triangle& inTriangle,
                           Point& outIntersectionPoint,
                                     int& intersection){

    // implemented after: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    // very detailed explaination of the algorithm here:
    // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection

    using namespace Eigen;
    const double EPSILON = 0.0000001;

    // init "Eigen" vectors
    Vector3d rayO(rayOrigin.x(), rayOrigin.y(), rayOrigin.z());
    Vector3d rayV(rayVector.x(), rayVector.y(), rayVector.z());

    Vector3d vertex0(inTriangle.vertex(0).x(), inTriangle.vertex(0).y(), inTriangle.vertex(0).z());
    Vector3d vertex1(inTriangle.vertex(1).x(), inTriangle.vertex(1).y(), inTriangle.vertex(1).z());
    Vector3d vertex2(inTriangle.vertex(2).x(), inTriangle.vertex(2).y(), inTriangle.vertex(2).z());

    Vector3d edge1, edge2, pvec, s, q;
    double det;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    pvec = rayV.cross(edge2);
    det = edge1.dot(pvec);

    // equation could also be expressed as
    // (it is simply written in the way above to reuse pvec and det results later on)
    // det = edge1.cross(edge2) .dot(rayV)
    // which simply means, if triangle normal and ray vector have 0 dot product
    // (see line below), they are obvioulsy parallel

    double f = 1.0/det;
    s = rayO - vertex0;
    double u = f * s.dot(pvec);    // barycentric coordinate u

    q = s.cross(edge1);
    double v = f * rayV.dot(q); // barycentric coordinate v

    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * edge2.dot(q);
    outIntersectionPoint = Point(u,v,t);


    if (det > -EPSILON && det < EPSILON) // here I could also decide what to do if ray intersect triangle from "behind"
        return;    // This ray is parallel to this triangle.
    if (u < 0.0 || u > 1.0)     // if u is not between 0 and 1, intersection point does not lie in triangle
        return;
    if (v < 0.0 || u + v > 1.0) // if v is not between 0 and 1, intersection point does not lie in triangle
        return;
    intersection = 1;

}



//double pointDistance(const Tree& tree, const Point &query){
//    // Construct AABB tree with a KdTree

//    Point_and_primitive_id primitive_id = tree.closest_point_and_primitive(query);
//    Polyhedron::Halfedge_around_facet_circulator facerunner = primitive_id.second->facet_begin();

//    Point p1;
//    Point p2;
//    Point p3;

//    p1 = facerunner->vertex()->point();
//    facerunner++;
//    p2 = facerunner->vertex()->point();
//    facerunner++;
//    p3 = facerunner->vertex()->point();

//    EPICK::Plane_3 plane(p1,p2,p3);

//    // Determine the side and return true if inside!
//    return pointPlaneDistance(plane,query);
//}



tetFeatures calcTetFeatures(const Tetrahedron& tet){

    tetFeatures tF;
    tF.radius = sqrt(CGAL::squared_distance(tet.vertex(0), CGAL::circumcenter(tet)));
    tF.vol = tet.volume();
    vector<double> edge_lengths;
    edge_lengths.push_back(sqrt(CGAL::squared_distance(tet.vertex(0),tet.vertex(1))));
    edge_lengths.push_back(sqrt(CGAL::squared_distance(tet.vertex(0),tet.vertex(2))));
    edge_lengths.push_back(sqrt(CGAL::squared_distance(tet.vertex(0),tet.vertex(3))));
    edge_lengths.push_back(sqrt(CGAL::squared_distance(tet.vertex(1),tet.vertex(2))));
    edge_lengths.push_back(sqrt(CGAL::squared_distance(tet.vertex(1),tet.vertex(3))));
    edge_lengths.push_back(sqrt(CGAL::squared_distance(tet.vertex(2),tet.vertex(3))));
    sort(edge_lengths.begin(), edge_lengths.end());
    tF.longest_edge = edge_lengths[5];
    tF.shortest_edge = edge_lengths[0];

    return tF;

}
//tetFeatures calcTetFeatures(const Tetrahedron& tet){

//    tetFeatures tF;
//    tF.radius = (CGAL::squared_distance(tet.vertex(0), CGAL::circumcenter(tet)));
//    tF.vol = tet.volume();
//    vector<double> edge_lengths;
//    edge_lengths.push_back((CGAL::squared_distance(tet.vertex(0),tet.vertex(1))));
//    edge_lengths.push_back((CGAL::squared_distance(tet.vertex(0),tet.vertex(2))));
//    edge_lengths.push_back((CGAL::squared_distance(tet.vertex(0),tet.vertex(3))));
//    edge_lengths.push_back((CGAL::squared_distance(tet.vertex(1),tet.vertex(2))));
//    edge_lengths.push_back((CGAL::squared_distance(tet.vertex(1),tet.vertex(3))));
//    edge_lengths.push_back((CGAL::squared_distance(tet.vertex(2),tet.vertex(3))));
//    sort(edge_lengths.begin(), edge_lengths.end());
//    tF.longest_edge = edge_lengths[5];
//    tF.shortest_edge = edge_lengths[0];

//    return tF;

//}






