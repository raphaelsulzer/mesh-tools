#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <processing/tetIntersection.h>
#include <util/vectorArithmetic.h>


int tetIntersectionVol(Tetrahedron& tet,
                            std::vector<Plane_cgal>& all_planes,
                            double& vol,
                            int plane_count,
                            std::string tet_name){


//    exportOFF(tet, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dtet");
    bool exp = false;
//    bool exp = true;

    // path for intermediate tet export
    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/tet_";
//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/tet";

//    std::cout << plane_count << std::endl;

    // when tet was recursivly called it goes to plane 4
    // which doesn't exist; hence algorithm is over and
    // the volume of the current tet should be added
    if(plane_count == 4){
        double newVol = abs(tet.volume());
//        std::cout << "vol counted tet: " << tet_name << std::endl;
        vol+=newVol;
        return 1;
    }

    // get the current plane from the vector of planes
    Plane_cgal plane = all_planes[plane_count];
    // check if there are points of the Dt tet on the positive AND negative side of the
    // current plane
    std::vector<int> neg;
    std::vector<int> pos;
    std::vector<int> close;
    double epsilon = 0.000000001;
    for(int v = 0; v < 4; v++){
        // test if point is on positive or negative side
        // only if there are points on both sides of this plane I need to do sth with this plane
        Point point = tet.vertex(v);
        double dist = pointPlaneDistance(plane, point);
        if(dist < -epsilon){
            neg.push_back(v);
        }
        else if(dist > epsilon){
            pos.push_back(v);
        }
        else{
            close.push_back(v);
        }
    }

    // continue according to intersection with the current plane
    int ns = neg.size();
    int ps = pos.size();
    int cs = close.size();
    // if there was an intersection left and right of the plane
    if(ns != 0 && ps != 0){
        // initialize the new tet
        Tetrahedron newTet;
        if(ps > ns)
        {
            std::vector<Tetrahedron> splitTets;
            // first make the 3 outside split tets, which are not really needed,
            // but its intersection points are needed, for afterwards making the single inside tet
            for(int n = 0; n < ns; n++){
                for(int p = 0; p < ps; p++){
                    Point p1 = tet.vertex(neg[n]);
                    Point p2 = tet.vertex(pos[p]);
                    Line l(p1,p2);
                    CGAL::cpp11::result_of<Intersect(Line, Plane_cgal)>::type
                                  result = intersection(l, plane);
                    if(result){
                        if(const Point* intersection_point = boost::get<Point>(&*result)){
                            // make a new tet with all the negative points
                            // and replace the positive one with the new intersection point
                            std::vector<Point> newTetPoints;
                            // first add the intersection points of the previous tetrahedron
                            if(splitTets.size() > 0){
                                for(int st = 0; st < splitTets.size(); st++){
                                    newTetPoints.push_back(splitTets[st].vertex(3));
                                }
                            }
                            // add the positive point of the line that was split
                            newTetPoints.push_back(p2);
                            // add other positive points (will get less as the iteration moves on)
                            for(int pp = p+1; pp < ps; pp++)
                                newTetPoints.push_back(tet.vertex(pos[pp]));
                            // add close points
                            for(int c = 0; c < close.size(); c++)
                                newTetPoints.push_back(tet.vertex(close[c]));
                            newTetPoints.push_back(*intersection_point);
//                            Tetrahedron newTet(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            // add tet to split tets
                            // TODO; split tets do not have to be actual tets, it will be enough to save the intersection point
                            splitTets.push_back(newTet);
                        }
                    }
                }
            }
            // now form the remaining (inside) polyhedron and only call the tetIntersectionFun with that one
            std::vector<Point> newTetPoints;
            // add all the intersection points
            for(int st = 0; st < splitTets.size(); st++){
                newTetPoints.push_back(splitTets[st].vertex(3));
            }
            // need to add close points here as well,
            // because a shared edge of a split tet and Dt tet
            // leeds to only having 2 split tets, which in this positive sided
            // tetrahedron case here would not leed to a valid tetrahedron
            // because there are only two points coming from the loop above,
            // and not 3
            for(int c = 0; c < close.size(); c++)
                newTetPoints.push_back(tet.vertex(close[c]));
            // add the only negative point of this plane-cut
            newTetPoints.push_back(tet.vertex(neg[0]));
            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
            std::string tet_name = std::to_string(plane_count)+"_pos";
            if(exp)
                exportOFF(newTet, path+tet_name);
            // this is now a SINGLE tet that is completely inside the negative side of one plane
            // note: the difference here to the next else (ns > ps) is that there it is not possible
            // to form a single tet from the remaining cut lying to the inside of the one plane
            tetIntersectionVol(newTet, all_planes, vol, plane_count+1, tet_name);
        }
        else if(ns > ps)
        {
            std::vector<Tetrahedron> splitTets;
            for(int n = 0; n < ns; n++){
                for(int p = 0; p < ps; p++){
                    Point p1 = tet.vertex(neg[n]);
                    Point p2 = tet.vertex(pos[p]);
                    Line l(p1,p2);
                    CGAL::cpp11::result_of<Intersect(Line, Plane_cgal)>::type
                                  result = intersection(l, plane);
                    if(result){
                        if(const Point* intersection_point = boost::get<Point>(&*result)){
                            // make a new tet with all the negative points
                            // and always replace the positive one with the new intersection point
                            std::vector<Point> newTetPoints;
                            // first add the intersection points of the previous tetrahedron
                            // in the first iteration there will be no points yet
                            if(splitTets.size() > 0){
                                for(int st = 0; st < splitTets.size(); st++){
                                    newTetPoints.push_back(splitTets[st].vertex(3));
                                }
                            }
                            // add the negative point of the line that was split as the first point of the split tet
                            newTetPoints.push_back(p1);
                            // add other negative points (will get less as the iteration moves on)
                            for(int nn = n+1; nn < ns; nn++)
                                newTetPoints.push_back(tet.vertex(neg[nn]));
                            // add close points
                            for(int c = 0; c < close.size(); c++)
                                newTetPoints.push_back(tet.vertex(close[c]));
                            newTetPoints.push_back(*intersection_point);
//                            Tetrahedron newTet(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            // ////export
                            std::string tet_name = std::to_string(plane_count)+"_"+std::to_string(n)+std::to_string(p);
                            if(exp)
                                exportOFF(newTet, path+tet_name);
                            splitTets.push_back(newTet);
                            tetIntersectionVol(newTet, all_planes, vol, plane_count+1, tet_name);
                        }
                    }
                }
            }
        }
        else{ // meaning ns == ps
            // make a vector where all 4 intersection points are saved
            std::vector<Point> intersectionPoints;
            for(int n = 0; n < ns; n++){
                // push back the two negative points
                for(int p = 0; p < ps; p++){
                    Point p1 = tet.vertex(neg[n]);
                    Point p2 = tet.vertex(pos[p]);
                    Line l(p1,p2);
                    CGAL::cpp11::result_of<Intersect(Line, Plane_cgal)>::type
                                  result = intersection(l, plane);
                    if(result){
                        if(const Point* intersection_point = boost::get<Point>(&*result)){
                            // close points
//                            for(int c = 0; c < close.size(); c++)
//                                newTetPoints.push_back(tet.vertex(close[c]));
                            intersectionPoints.push_back(*intersection_point);
                        }
                    }
                }
            }
            for(int c = 0; c < close.size(); c++)
                intersectionPoints.push_back(tet.vertex(close[c]));
            // now that all intersection points are calculated, for the 3 inside tets
            std::vector<Tetrahedron> splitTets;
            std::vector<std::string> tetNames;
            // first make the big one with two intersection points and two negative points
            Tetrahedron newTet(intersectionPoints[0], intersectionPoints[1], tet.vertex(neg[0]), tet.vertex(neg[1]));
            //export
            tetNames.push_back(std::to_string(plane_count)+"_big");
            if(exp)
                exportOFF(newTet, path+tetNames[0]);
            splitTets.push_back(newTet);
            // now the first small one
            Tetrahedron newTet1(intersectionPoints[0], intersectionPoints[1], intersectionPoints[2], tet.vertex(neg[1]));
            tetNames.push_back(std::to_string(plane_count)+"_small1");
            if(exp)
                exportOFF(newTet1, path+tetNames[1]);
            splitTets.push_back(newTet1);
            // now the second small one
            Tetrahedron newTet2(intersectionPoints[1], intersectionPoints[2], intersectionPoints[3], tet.vertex(neg[1]));
            tetNames.push_back(std::to_string(plane_count)+"_small2");
            if(exp)
                exportOFF(newTet2, path+tetNames[2]);
            splitTets.push_back(newTet);
            // now continue with looping over planes for all three tets
            for(int t = 0; t < splitTets.size(); t++){
                tetIntersectionVol(splitTets[t], all_planes, vol, plane_count+1, tetNames[t]);
            }
        } // end of 2-2 case
    } // end of n != 0 && p != 0
    else if(ps == 0 && ns != 0){           // if positive side is empty, means no intersection with this plane
        tetIntersectionVol(tet, all_planes, vol, plane_count+1, tet_name);
    }
    // if ANY plane has negative side is empty (with the original tet!), means no intersection at all
    // of course it can be with cut tets; but then it will continue in another recursion loop
    // also works if one (ore more) points are close, and the rest empty
    else{
        return 0;
    }
    return 1;
}




void tet2Planes(const Tetrahedron& tet, std::vector<Plane_cgal>& planes){

    Point centroid = CGAL::centroid(tet[0],tet[1],tet[2],tet[3]);

//    std::vector<Plane_cgal> planes(4);
    planes[0] = Plane_cgal(tet[0],tet[2],tet[1]);
    planes[1] = Plane_cgal(tet[0],tet[1],tet[3]);
    planes[2] = Plane_cgal(tet[1],tet[2],tet[3]);
    planes[3] = Plane_cgal(tet[0],tet[3],tet[2]);
    for(int i = 0; i<4; i++){
        if(!planes[i].has_on_negative_side(centroid)){
            planes[i]=planes[i].opposite();
        }
    }
}


//void testTetIntersectionVol(){

////    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/";
//    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/tet_";
//    system("exec rm -r /home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/*");

////// Toy examples
//////    Point sp01(1.3,1.3,1.3);
//////    Point sp02(0.3,3.3,1.3);
//////    Point sp03(2.4,2,1.3);

////    Point sp0(1,1,1);
////    Point sp1(0,4,1); // toy example 1
//////    Point sp1(1,3,1);   // toy example 2
////    Point sp2(3,2,1);
////    Point sp3(2,2,4);

//////    Polyhedron six;
//////    six.make_tetrahedron(sp0,sp1,sp2,sp01);
//////    six.make_tetrahedron(sp01,sp02,sp03,sp1);
//////    six.make_tetrahedron(sp02,sp03,sp0,sp2);
//////    exportOFF(six, path+"six");

////    Point spc = CGAL::centroid(sp0,sp1,sp2,sp3);
////    Tetrahedron stet(sp0,sp1,sp2,sp3);
////    Polyhedron sp;
////    sp.make_tetrahedron(sp0,sp1,sp2,sp3);
////    exportOFF(sp, path+"sp");
////    std::vector<Plane_cgal> planes(4);
////    planes[0] = Plane_cgal(sp0,sp2,sp1);
////    planes[1] = Plane_cgal(sp0,sp1,sp3);
////    planes[2] = Plane_cgal(sp1,sp2,sp3);
////    planes[3] = Plane_cgal(sp0,sp3,sp2);
////    for(int i = 0; i<4; i++){
////        if(!planes[i].has_on_negative_side(spc)){
////            planes[i]=planes[i].opposite();
////        }
////    }

//////    Point dp0(3,1,5);
//////    Point dp1(3,3,5);
//////    Point dp2(1,2,5);
//////    Point dp3(2,2,2);
////    // example 3
////    Point dp0(0,2,5);
////    Point dp1(0,4,5);
////    Point dp2(2,3,5);
////    Point dp3(2,2,0);
////    Tetrahedron dtet(dp0,dp1,dp2,dp3);
////    Polyhedron dp;
////    dp.make_tetrahedron(dp0,dp1,dp2,dp3);
////    exportOFF(dp, path+"dp");

//    // examples from TLS
//    Tetrahedron dtet;
//    importOff("/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp1.off", dtet);
//    std::vector<Plane_cgal> planes;
//    importOff("/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp1.off", planes);

//    // intersection with own function
//    double cvol = 0.0;
//    int pc = 0;
//    int ib = tetIntersectionFun(dtet, planes, cvol, pc);

//    // NEF intersection
//    Polyhedron dp;
//    importOff("/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp1.off", dp);
//    Polyhedron sp;
//    importOff("/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp1.off", sp);

//    CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> sensor_modifier(sp);
//    Polyhedron_Exact sp_exact;
//    sp_exact.delegate(sensor_modifier);
//    Nef_polyhedron snef(sp_exact);

//    CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> dp_modifier(dp);
//    Polyhedron_Exact dp_exact;
//    dp_exact.delegate(dp_modifier);
//    Nef_polyhedron dnef(dp_exact);

//    Nef_polyhedron fullnef = dnef*snef;

//    Polyhedron_Exact intersection_tet_exact;
//    fullnef.convert_to_polyhedron(intersection_tet_exact);

//    CGAL::Polyhedron_copy_3<Polyhedron_Exact, Polyhedron::HalfedgeDS> modifier_rev(intersection_tet_exact);
//    Polyhedron full_poly;
//    full_poly.delegate(modifier_rev);
//    exportOFF(full_poly, path+"intersection");

//    double vol_full1 = CGAL::Polygon_mesh_processing::volume(full_poly);

//    // comparision
//    std::cout << "intersection: " << ib << " my vol: " << cvol << "    nef vol: " << vol_full1 << std::endl;

//}






//-----------------------------------------------------------------------------
//-------------------------------  TEST FUN  ----------------------------------
//-----------------------------------------------------------------------------
bool
tetIntersectionTest
(const Tetrahedron& tetrahedron_0,
 const Tetrahedron& tetrahedron_1)
{
  // This algorithm checks whether two tetrahedra intersect.

  // Algorithm and source code from Fabio Ganovelli, Federico Ponchio
  // and Claudio Rocchini: Fast Tetrahedron-Tetrahedron Overlap
  // Algorithm, Journal of Graphics Tools, 7(2), 2002. DOI:
  // 10.1080/10867651.2002.10487557. Source code available at
  // http://web.archive.org/web/20031130075955/http://www.acm.org/jgt/papers/GanovelliPonchioRocchini02/tet_a_tet.html


  // Get the vertices as points
  std::vector<Point> V1(4), V2(4);
  for (std::size_t i = 0; i < 4; ++i)
  {
    V1[i] = tetrahedron_0.vertex(i);
    V2[i] = tetrahedron_1.vertex(i);;
  }

  // Get the vectors between V2 and V1[0]
  std::vector<Point> P_V1(4);
  for (std::size_t i = 0; i < 4; ++i){
      Vector temp = V2[i]-V1[0];
      P_V1[i] = Point(temp.x(), temp.y(), temp.z());

  }

  // Data structure for edges of V1 and V2
  std::vector<Vector> e_v1(5), e_v2(5);
  e_v1[0] = V1[1] - V1[0];
  e_v1[1] = V1[2] - V1[0];
  e_v1[2] = V1[3] - V1[0];
  Point n = crossP(e_v1[1], e_v1[0]);

  // Maybe flip normal. Normal should be outward.
  if (dot(n, e_v1[2]) > 0)
    n = neg(n);
  std::vector<int> masks(4);
  std::vector<std::vector<double>> Coord_1(4, std::vector<double>(4));
  if (separating_plane_face_A_1(P_V1, n, Coord_1[0], masks[0]))
    return false;
  n = crossP(e_v1[0],e_v1[2]);

  // Maybe flip normal
  if (dot(n, e_v1[1]) > 0)
    n = neg(n);
  if (separating_plane_face_A_1(P_V1, n, Coord_1[1], masks[1]))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 0, 1))
    return false;
  n = crossP(e_v1[2], e_v1[1]);

  // Maybe flip normal
  if (dot(n, e_v1[0]) > 0)
      n = neg(n);
  if (separating_plane_face_A_1(P_V1, n, Coord_1[2], masks[2]))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 0, 2))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 1,2))
    return false;
  e_v1[4] = V1[3] - V1[1];
  e_v1[3] = V1[2] - V1[1];
  n = crossP(e_v1[3], e_v1[4]);

  // Maybe flip normal. Note the < since e_v1[0]=v1-v0.
  if (dot(n, e_v1[0]) < 0)
      n = neg(n);
  if (separating_plane_face_A_2(V1, V2, n, Coord_1[3], masks[3]))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 0, 3))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 1, 3))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 2, 3))
    return false;
  if ((masks[0] | masks[1] | masks[2] | masks[3] )!= 15)
    return true;

  // From now on, if there is a separating plane, it is parallel to a
  // face of b.
  std::vector<Point> P_V2(4);
  for (std::size_t i = 0; i < 4; ++i){
      Vector temp = V1[i] - V2[0];
      P_V2[i] = Point(temp.x(), temp.y(), temp.z());
  }

  e_v2[0] = V2[1] - V2[0];
  e_v2[1] = V2[2] - V2[0];
  e_v2[2] = V2[3] - V2[0];
  n = crossP(e_v2[1], e_v2[0]);

  // Maybe flip normal
  if (dot(n, e_v2[2])>0)
    n = neg(n);
  if (separating_plane_face_B_1(P_V2, n))
    return false;
  n=crossP(e_v2[0], e_v2[2]);

  // Maybe flip normal
  if (dot(n, e_v2[1]) > 0)
    n = neg(n);
  if (separating_plane_face_B_1(P_V2, n))
    return false;
  n = crossP(e_v2[2], e_v2[1]);

  // Maybe flip normal
  if (dot(n, e_v2[0]) > 0)
    n = neg(n);
  if (separating_plane_face_B_1(P_V2, n))
    return false;
  e_v2[4] = V2[3] - V2[1];
  e_v2[3] = V2[2] - V2[1];
  n = crossP(e_v2[3], e_v2[4]);

  // Maybe flip normal. Note the < since e_v2[0] = V2[1] - V2[0].
  if (dot(n, e_v2[0]) < 0)
    n = neg(n);
  if (separating_plane_face_B_2(V1, V2, n))
    return false;

  return true;
}
//-----------------------------------------------------------------------------
bool
separating_plane_face_A_1(const std::vector<Point>& pv1,
                          const Point& n,
                          std::vector<double>& coord,
                          int&  mask_edges)
{
  // Helper function for tetrahedron-tetrahedron collision test:
  // checks if plane pv1 is a separating plane. Stores local
  // coordinates and the mask bit mask_edges.

  mask_edges = 0;
  const int shifts[4] = {1, 2, 4, 8};

  for (std::size_t i = 0; i < 4; ++i)
  {
    coord[i] = dot(pv1[i],n);
    if (coord[i] > 0)
      mask_edges |= shifts[i];
  }

  return (mask_edges == 15);
}
//-----------------------------------------------------------------------------
bool
separating_plane_face_A_2(const std::vector<Point>& V1,
                          const std::vector<Point>& V2,
                          const Point& n,
                          std::vector<double>& coord,
                          int&  mask_edges)
{
  // Helper function for tetrahedron-tetrahedron collision test:
  // checks if plane v1,v2 is a separating plane. Stores local
  // coordinates and the mask bit mask_edges.

  mask_edges = 0;
  const int shifts[4] = {1, 2, 4, 8};

  for (std::size_t i = 0; i < 4; ++i)
  {
    coord[i] = dot((V2[i] - V1[1]),n);
    if (coord[i] > 0)
      mask_edges |= shifts[i];
  }

  return (mask_edges == 15);
}
bool separating_plane_face_B_1(const std::vector<Point>& P_V2,
                  const Point& n)
{
  return ((dot(P_V2[0],n) > 0) &&
      (dot(P_V2[1],n) > 0) &&
      (dot(P_V2[2],n) > 0) &&
      (dot(P_V2[3],n) > 0));
}

// Helper function for collides_tetrahedron_tetrahedron: checks if
// plane v1, v2 is a separating plane.
bool separating_plane_face_B_2(const std::vector<Point>& V1,
                  const std::vector<Point>& V2,
                  const Point& n)
{
  return ((dot((V1[0] - V2[1]),n) > 0) &&
      (dot((V1[1] - V2[1]),n) > 0) &&
      (dot((V1[2] - V2[1]),n) > 0) &&
      (dot((V1[3] - V2[1]),n) > 0));
}


//-----------------------------------------------------------------------------
bool
separating_plane_edge_A(
  const std::vector<std::vector<double>>& coord_1,
  const std::vector<int>& masks, int f0, int f1)
{
  // Helper function for tetrahedron-tetrahedron collision: checks if
  // edge is in the plane separating faces f0 and f1.

  const std::vector<double>& coord_f0 = coord_1[f0];
  const std::vector<double>& coord_f1 = coord_1[f1];

  int maskf0 = masks[f0];
  int maskf1 = masks[f1];

  if ((maskf0 | maskf1) != 15) // if there is a vertex of b
    return false; // included in (-,-) return false

  maskf0 &= (maskf0 ^ maskf1); // exclude the vertices in (+,+)
  maskf1 &= (maskf0 ^ maskf1);

  // edge 0: 0--1
  if ((maskf0 & 1) && // the vertex 0 of b is in (-,+)
      (maskf1 & 2)) // the vertex 1 of b is in (+,-)
    if ((coord_f0[1]*coord_f1[0] - coord_f0[0]*coord_f1[1]) > 0)
      // the edge of b (0,1) intersect (-,-) (see the paper)
      return false;

  if ((maskf0 & 2) &&
      (maskf1 & 1))
    if ((coord_f0[1]*coord_f1[0] - coord_f0[0]*coord_f1[1]) < 0)
      return false;

  // edge 1: 0--2
  if ((maskf0 & 1) &&
      (maskf1 & 4))
    if ((coord_f0[2]*coord_f1[0] - coord_f0[0]*coord_f1[2]) > 0)
      return false;

  if ((maskf0 & 4) &&
      (maskf1 & 1))
    if ((coord_f0[2]*coord_f1[0] - coord_f0[0]*coord_f1[2]) < 0)
      return false;

  // edge 2: 0--3
  if ((maskf0 & 1) &&
      (maskf1 & 8))
    if ((coord_f0[3]*coord_f1[0] - coord_f0[0]*coord_f1[3]) > 0)
      return false;

  if ((maskf0 & 8) &&
      (maskf1 & 1))
    if ((coord_f0[3]*coord_f1[0] - coord_f0[0]*coord_f1[3]) < 0)
      return false;

  // edge 3: 1--2
  if ((maskf0 & 2) &&
      (maskf1 & 4))
    if ((coord_f0[2]*coord_f1[1] - coord_f0[1]*coord_f1[2]) > 0)
      return false;

  if ((maskf0 & 4) &&
      (maskf1 & 2))
    if ((coord_f0[2]*coord_f1[1] - coord_f0[1]*coord_f1[2]) < 0)
      return false;

  // edge 4: 1--3
  if ((maskf0 & 2) &&
      (maskf1 & 8))
    if ((coord_f0[3]*coord_f1[1] - coord_f0[1]*coord_f1[3]) > 0)
      return false;

  if ((maskf0 & 8) &&
      (maskf1 & 2))
    if ((coord_f0[3]*coord_f1[1] - coord_f0[1]*coord_f1[3]) < 0)
      return false;

  // edge 5: 2--3
  if ((maskf0 & 4) &&
      (maskf1 & 8))
    if ((coord_f0[3]*coord_f1[2] - coord_f0[2]*coord_f1[3]) > 0)
      return false;

  if ((maskf0 & 8) &&
      (maskf1 & 4))
    if ((coord_f0[3]*coord_f1[2] - coord_f0[2]*coord_f1[3]) < 0)
      return false;

  // Now there exists a separating plane supported by the edge shared
  // by f0 and f1.
  return true;
}
//-----------------------------------------------------------------------------



