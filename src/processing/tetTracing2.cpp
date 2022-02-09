//#include <base/cgal_typedefs.h>
//#include <IO/fileIO.h>
//#include <processing/rayTracing.h>
//#include <processing/tetIntersection.h>

//namespace tetTracing2{

///////////////////////////////////////////////////////////////////
///////////////////////// Tetrahedron tracing /////////////////////
///////////////////////////////////////////////////////////////////
/////
/////
/////
/////
///// cannot simply look for the sensor tetrahedron from the same 3 points, because the connection is broken when combining different sensors. Meaning I will have
///// Delaunay surface triangles that are formed by points from different sensors.
/////
/////
//void idxCells(Delaunay& Dt){
//    Delaunay::Finite_cells_iterator fci;
//    int i=0;
//    for(fci=Dt.finite_cells_begin();fci!=Dt.finite_cells_end();fci++){
//        fci->info().global_idx = i++;
//    }
//}

//// TODO: measure the time of
//// do_intersect(tri, tri)
//// and tetIntersectionFun()

//int traverseCells(Delaunay& Dt,
//                  Cell_handle& first_cell,
//                  Cell_handle& current_cell, std::unordered_set<Cell_handle>& processed,
//                  std::vector<Plane>& planes,
//                  Triangle& sensor_tri,
//                  double noise,
//                  int outside_weight,
//                  double infinite_score){

//    // if current cell is not infinite then go on
//    // processed state is already checked before calling this function, so no need to check again
//    if(!Dt.is_infinite(current_cell))
//    {

////        int fidx = first_cell->info().global_idx;
////        int cidx = current_cell->info().global_idx;
////        std::cout << "first cell: " << fidx << "    second cell: " << cidx << std::endl;

//        // TODO: fix this!!
//        // this intersection test is actually not valid:
//        // a Dt.tetrahedron can lie completely inside a sensor tetrahedron and thus
//        // not be detected by this test!!!
//        bool dointersect = false;
//        int t=0;
//        while(!dointersect && t < 4)
//        {
//            Triangle tri = Dt.triangle(current_cell, t++);
//            dointersect = CGAL::do_intersect(sensor_tri,tri);
//        }
//        if(!dointersect){
////            processed.insert(current_cell);
//            return 0;
//        }

//        Tetrahedron current_tet = Dt.tetrahedron(current_cell);
//        double vol = 0.0;
//        tetIntersectionFun(current_tet, planes, vol);
//        if(!isnan(vol)){
//            double score = vol*outside_weight/abs(current_tet.volume());
////            double score = vol*outside_weight;
//            current_cell->info().outside_score+=score;
//        }
//        else{
//            std::cout << "NaN hit. Intersection? " << std::endl;
////            exportOFF(sp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp");
//            exportOFF(current_tet, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp");
//            processed.insert(current_cell);
//            return 0;
//        }
//        processed.insert(current_cell);

//        // get neighbouring cells
//        for(int ci = 0; ci < 4; ci++){
//            Facet fac = std::make_pair(current_cell, ci);
//            Facet mirror_fac = Dt.mirror_facet(fac);
//            Cell_handle newCell = mirror_fac.first;
//            if(processed.find(newCell) == processed.end()){
//                traverseCells(Dt, first_cell, newCell, processed, planes, sensor_tri,
//                              noise, outside_weight, infinite_score);
//            }
//        }
//    }
//    else{
//        // give infinite cell high value
//        current_cell->info().outside_score+=infinite_score;
//        processed.insert(current_cell);
//    }
//    return 1;
//}


////int whichSide(std::vector<Point> S, Point D, Point P){
////    // S v e r t i c e s a r e p r o j e c t e d t o t h e form P+t ∗D .
////    // R e t u r n v a l u e i s +1 i f a l l t > 0 ,
////    // −1 i f a l l t < 0 , 0 o t h e r w i s e , i n w h i c h c a s e t h e l i n e s p l i t s t h e p o l y g o n .
////    int positive = 0;
////    int negative = 0;
////    for(int i = 0; i < S.size(); i++){

////        Vector temp = S[i] - P;
////        double t = D.x() * temp.x() + D.y() * temp.y() + D.z() * temp.z();
////        if( t > 0 )
////            positive++;
////        else if( t < 0 )
////            negative++;
////        if(positive && negative )
////            return 0;
////    }
////    return( positive? +1 : -1 );
////}

////std::vector<Vector> getTetNormals(const Tetrahedron& tet){

////    std::vector<Vector> normals(4);

////    std::vector<Point> tet_pts;
////    tet_pts.push_back(tet.vertex(0));
////    tet_pts.push_back(tet.vertex(1));
////    tet_pts.push_back(tet.vertex(2));
////    tet_pts.push_back(tet.vertex(3));

////    Vector v1, v2;
////    // 0nd
////    v1 = tet.vertex(3) - tet.vertex(1);
////    v2 = tet.vertex(2) - tet.vertex(1);
////    normals[0] = cross(v1,v2);
////    // 1nd
////    v1 = tet.vertex(3) - tet.vertex(0);
////    v2 = tet.vertex(2) - tet.vertex(0);
////    normals[1] = cross(v1,v2);
////    // 2nd
////    v1 = tet.vertex(1) - tet.vertex(0);
////    v2 = tet.vertex(3) - tet.vertex(0);
////    normals[2] = cross(v1,v2);
////    // 2nd
////    v1 = tet.vertex(1) - tet.vertex(0);
////    v2 = tet.vertex(2) - tet.vertex(0);
////    normals[3] = cross(v1,v2);

////    return normals;
////}



////bool testTetIntersection(Tetrahedron C0, Tetrahedron C1){
////    // Test faces of C0 for seperation.
////    // Because of the counter-clockwise ordering, the projection interval for C0 is [m, 0] where m <= 0 .
////    // Only try to deterline if C1 is on the 'positive' side of the line

////    std::vector<Point> P0;
////    std::vector<Point> P1;
////    for(int i = 0; i < 4; i++){
////        P0.push_back(C0.vertex(i));
////        P1.push_back(C1.vertex(i));
////    }

////    for(int i = 0; i < 4; i++){


////        Point D =
////        D = C0 . F ( i ) . n o r m a l ;
////        // outward pointing
////        if( whichSide( C1 . V , D, C0 . F ( i ) . v e r t e x ) > 0 ){
////        // C1 i s e n t i r e l y on ‘ p o s i t i v e ’ s i d e o f l i n e C0 . F ( i ) . v e r t e x+t ∗D
////            return false ;
////        }
////    }
////    // T e s t f a c e s o f C1 f o r s e p a r a t i o n .
////    // Because of the counter-clockwise ordering
////    // t h e p r o j e c t i o n i n t e r v a l f o r C1 i s [m, 0 ] w h e r e m <= 0 .
////    // Only t r y t o d e t e r m i n e
////    // i f C0 i s on t h e ‘ p o s i t i v e ’ s i d e o f t h e l i n e .
////    for(int i = 0 ; i < C1 . L ; i ++){
////        D = C1 . F ( i ) . n o r m a l ;
////        // o u t w a r d p o i n t i n g
////        if( whichSide( C0.V, D, C1.F(i).vertex ) > 0 ){
////            // C0 i s e n t i r e l y on ‘ p o s i t i v e ’ s i d e o f l i n e C1 . F ( i ) . v e r t e x+t ∗D
////            return false;
////        }
////    }
////    // T e s t c r o s s p r o d u c t o f p a i r s o f e d g e s , one from e a c h p o l y h e d r o n .
////    for(int i = 0 ; i < C0.M; i ++){
////        for(int j = 0 ; j < C1.M; j ++){
////            D = Cross( C0.E ( i ) , C1 . E ( j ) ) ;
////            int side0 = WhichSide(C0.V, D, C0.E(i).vertex ) ;
////            if(side0 == 0 ){
////                continue;
////            }
////            int side1 = WhichSide( C1 . V , D, C0 . E ( i ) . v e r t e x ) ;
////            if( side1 == 0 ){
////            continue ;
////            }
////            if( side0 * side1 < 0 ){
////            // C0 and C1 a r e on ‘ o p p o s i t e ’ s i d e s o f l i n e C0 . E ( i ) . v e r t e x+t ∗D
////            return false ;
////            }
////        }
////    }
////    return true ;
////}






//void firstCell(Delaunay& Dt, std::vector<std::vector<int>>& sensor_polys, int outside_weight, int infinite_score=1){
//    std::cout << "Start marking crossed tets..." << std::endl;
//    auto start = std::chrono::high_resolution_clock::now();

//    // save Delaunay vertex handle in sensor_infos vector for each sensor poly
//    std::vector<std::vector<Vertex_handle>> sensor_infos(sensor_polys.size());
//    Delaunay::Finite_vertices_iterator vit;
//    for(vit = Dt.finite_vertices_begin(); vit != Dt.finite_vertices_end(); vit++){
//        std::vector<int> sensor_tets = vit->info().sensor_tet;
//        for(int i = 0; i < sensor_tets.size(); i++){
//            int current_sensor_tet = sensor_tets[i];
//            sensor_infos[current_sensor_tet].push_back(vit);
//            int a = 5;
//        }
//    }
//    idxCells(Dt);


//    for(;;){





//    }






//    // iterate over all sensor triangles/tetrahedrons
//    int hit = 0;
///*    for(int k = 0; k < sensor_polys.size(); k++){
//        std::unordered_set<Cell_handle> processed;
////        std::cout << "sensor poly " << k << std::endl;

//        if(sensor_infos[k].size()<3)
//            continue;

//        // iterate over all 3 points of the sensor triangle
//        for(int j = 0; j < 3; j++){

//            // get the corresponding Dt vertex handle for the current point
//            // so we are now considering the Dt vertex sensor_infos[k][j]
//            Vertex_handle current_vertex = sensor_infos[k][j];
//            Point a(current_vertex->point());
//            Point b(sensor_infos[k][(j+1)%3]->point());
//            Point c(current_vertex->info().sensor_pos);
//            Triangle tri(a,b,c);
//            Ray ray(current_vertex->point(), current_vertex->info().sensor_vec);
//            // vector of incident cells to the current vertex
//            std::vector<Cell_handle> inc_cells;
//            Dt.incident_cells(current_vertex, std::back_inserter(inc_cells));
//            for(std::size_t c=0; c < inc_cells.size(); c++){
//                Cell_handle current_cell = inc_cells[c];
////                int fcidx = current_cell->info().global_idx;
////                std::cout << "first cell: " << fcidx << std::endl;
//                // if current cell is not infinite and is not already processed (for this sensor poly) then go on
//                if(!Dt.is_infinite(current_cell)){
//                    if(processed.find(current_cell) == processed.end())
//                    {
//                        int cellBasedVertexIndex = current_cell->index(current_vertex);
//                        Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);

//                        bool dointersect = CGAL::do_intersect(ray,tri);
//                        if(dointersect){
//    //                        if(processed.find(current_cell) == processed.end()){
//                                // make a polyhedron of the sensor tet
//                                // not actually necessary anymore, just for outputting it, while debugging
//                                Point sp0 = sensor_infos[k][0]->point();
//                                Point sp1 = sensor_infos[k][1]->point();
//                                Point sp2 = sensor_infos[k][2]->point();
//                                Point sp3 = sensor_infos[k][j]->info().sensor_pos;
//                                Point spc = CGAL::centroid(sp0,sp1,sp2,sp3);
//                                Polyhedron sp;
//                                sp.make_tetrahedron(sp0,sp1,sp2,sp3);

//                                // The plane is oriented such that p, q and r are oriented in a positive sense (that is counterclockwise) when seen from the positive side of h.
//                                // from: https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html
//                                // and tetrahedron orientation can be found here: https://doc.cgal.org/latest/Triangulation_3/index.html
//                                // and the cell centroid has to be on the NEGATIVE side of the plane
//                                // sensor planes
//                                std::vector<Plane> planes(4);
//                                planes[0] = Plane(sp0,sp2,sp1);
//                                planes[1] = Plane(sp0,sp1,sp3);
//                                planes[2] = Plane(sp1,sp2,sp3);
//                                planes[3] = Plane(sp0,sp3,sp2);
//                                for(int i = 0; i<4; i++){
//                                    if(!planes[i].has_on_negative_side(spc)){
//                                        planes[i]=planes[i].opposite();
//                                    }

//                                }

//                                Tetrahedron current_tet = Dt.tetrahedron(current_cell);
//                                // get the noise dependant weight of the current cell
//                                double noise = current_cell->vertex(0)->info().sigma +
//                                                current_cell->vertex(1)->info().sigma +
//                                                current_cell->vertex(2)->info().sigma +
//                                                current_cell->vertex(3)->info().sigma;

//                                // intersection
//                                double vol = 0;
//                                tetIntersectionFun(current_tet, planes, vol);

//                                // TODO: investigate why there are sometimes NaNs!!
//                                if(!isnan(vol)){
////                                    double score = vol*outside_weight*noise/abs(current_tet.volume());
//                                    double score = vol*outside_weight/abs(current_tet.volume());
////                                    double score = vol*outside_weight;
//                                    current_cell->info().outside_score+=score;
//    //                                std::cout << "inside score: " << current_cell->info().inside_score <<
//    //                                             "  outside score " << current_cell->info().outside_score << std::endl;
//    //                                exportOFF(sp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp");
//    //                                exportOFF(current_tet, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp");
//                                }
//                                else{
//                                    std::cout << "NaN hit. Intersection? " << std::endl;
//    //                                exportOFF(sp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp");
//    //                                exportOFF(current_tet, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp");
//                                    processed.insert(current_cell);
//                                    continue;
//                                }
//                                processed.insert(current_cell);

//                                // everything works quite well, but one problem is also that sensor mesh has wholes!!
//                                // thats why I have the errors in the door example
//                                // in fact, bad cells have inside and outside scores
//                                // so change to probability somehow?
//                                // furthermore also try with theory that all can be carved out with outside score
//                                // and triangle behind it is inside

//                                // update:
//                                // working well with carving out outside triangles. However, is it possible that
//                                // I am still missing outside triangles, due to some weird circumstances.
//                                // could try something like TODO:
//                                // if(close && vol==0)
//                                //      still check neighbouring cells

//                                //// traverse neighbouring cells
//                                for(int ci = 0; ci < 4; ci++){
//                                    Facet fac = std::make_pair(current_cell, ci);
//                                    Facet mirror_fac = Dt.mirror_facet(fac);
//                                    Cell_handle newCell = mirror_fac.first;
//                                    if(processed.find(newCell) == processed.end()){
//                                        traverseCells(Dt, current_cell, newCell, processed, planes, tri,
//                                                      noise, outside_weight, infinite_score);
//                                    }
//                                }
//    //                        }//end of IF-already-processed
//                        }
//                    }
//                } // if cell is not infinite end
//                else{
//                    // give infinite cell high value
//                    current_cell->info().outside_score+=infinite_score;
//                    processed.insert(current_cell);
//                }
//            }
//        }//end of iteration over all three triangles of a sensor tetrahedron
//    }*///end of iteration over all sensor tetrahedrons


//    auto stop = std::chrono::high_resolution_clock::now();
//    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    std::cout << "Tet intersection done in " << full_duration.count() << "s" << std::endl;

//}// end of firstCell function

//}// end of namespace



