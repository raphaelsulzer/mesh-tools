//#pragma once


//#include <base/cgal_typedefs.h>
//#include <IO/fileIO.h>
//#include <util/vectorArithmetic.h>



//// vol fun
//int tetIntersectionVol(Tetrahedron& tet,
//                            std::vector<EPICK::Plane_3>& all_planes,
//                            double& vol,
//                            int plane_count=0,
//                            std::string tet_name="");

//// vol helpers
//void tet2Planes(const Tetrahedron& tet, std::vector<EPICK::Plane_3>& planes);
////void testTetIntersectionVol();



//// test fun
//bool tetIntersectionTest(const Tetrahedron& tetrahedron_0,
//                                             const Tetrahedron& tetrahedron_1);

//// Helper function for collides_tetrahedron_tetrahedron: checks if
//// plane pv1 is a separating plane. Stores local coordinates bc
//// and the mask bit mask_edges.
//bool separating_plane_face_A_1(const std::vector<Point>& pv1,
//                  const Point& n,
//                  std::vector<double>& bc,
//                  int& mask_edges);

//// Helper function for collides_tetrahedron_tetrahedron: checks if
//// plane v1, v2 is a separating plane. Stores local coordinates bc
//// and the mask bit mask_edges.
//bool separating_plane_face_A_2(const std::vector<Point>& v1,
//                  const std::vector<Point>& v2,
//                  const Point& n,
//                  std::vector<double>& bc,
//                  int& mask_edges);

//// Helper function for collides_tetrahedron_tetrahedron: checks if
//// plane pv2 is a separating plane.
//bool separating_plane_face_B_1(const std::vector<Point>& P_V2,
//                  const Point& n);

//// Helper function for collides_tetrahedron_tetrahedron: checks if
//// plane v1, v2 is a separating plane.
//bool separating_plane_face_B_2(const std::vector<Point>& V1,
//                  const std::vector<Point>& V2,
//                  const Point& n);

//// Helper function for collides_tetrahedron_tetrahedron: checks if
//// edge is in the plane separating faces f0 and f1.
//static bool separating_plane_edge_A(const std::vector<std::vector<double> >& coord_1,
//                const std::vector<int>& masks,
//                int f0,
//                int f1);


