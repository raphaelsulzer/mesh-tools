//#include <exe/inputParser.h>

//#include <base/cgal_typedefs.h>
//#include <IO/fileIO.h>
//#include <IO/ttIO.h>
//#include <IO/ethIO.h>
//#ifdef OpenMVS
//#include <IO/mvIO.h>
//#endif
//#ifdef COLMAP
//#include <IO/colmapIO.h>
//#endif
//#include <util/helper.h>
//#include <util/vectorArithmetic.h>

//#include <processing/meshProcessing.h>
//#include <processing/edgeManifoldness.h>
//#include <processing/graphCut.h>

//#include <processing/pointSetProcessing.h>
//#include <processing/normalAndSensorProcessing.h>
//#include <processing/evaluation.h>
//#include <processing/rayTracingFacet.h>

//#include <learning/learning.h>
//#include <learning/learningMath.h>
//#include <learning/learningRayTracing.h>
//#include <learning/learningRayTracingGroundTruth.h>
//#include <learning/learningIO.h>

//#include <CGAL/refine_mesh_3.h>

//#ifdef Open3D
//#include "open3d/Open3D.h"
//#include "open3d/geometry/TetraMesh.h"
//#endif

//#include <CGAL/optimal_bounding_box.h>
//#include <boost/filesystem.hpp>
//using namespace boost::filesystem;

////#include <CGAL/Polygon_mesh_processing/clip.h>
////#include <CGAL/Surface_mesh_default_triangulation_3.h>
////#include <CGAL/Complex_2_in_triangulation_3.h>
////#include <CGAL/make_surface_mesh.h>
////#include <CGAL/Implicit_surface_3.h>
////#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
////#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

////Delaunay gDt;

////typedef Delaunay::Geom_traits GT;
////typedef GT::FT FT;

////// default triangulation for Surface_mesher
////typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
////// c2t3
////typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
//////typedef Tr::Geom_traits GT;
////typedef GT::Sphere_3 Sphere_3;
////typedef GT::FT FT;
////typedef FT (*Function)(Point);
////typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

////double smoothedOccupancy(Cell_handle& ch, Point& p, bool gc){

////    EPICK k;

////    // query point to triangle distances
////    double dp0 = sqrt(CGAL::internal::squared_distance_to_triangle(p,ch->vertex(0)->point(),ch->vertex(1)->point(),ch->vertex(2)->point(),k));
////    double dp1 = sqrt(CGAL::internal::squared_distance_to_triangle(p,ch->vertex(0)->point(),ch->vertex(1)->point(),ch->vertex(3)->point(),k));
////    double dp2 = sqrt(CGAL::internal::squared_distance_to_triangle(p,ch->vertex(1)->point(),ch->vertex(2)->point(),ch->vertex(3)->point(),k));
////    double dp3 = sqrt(CGAL::internal::squared_distance_to_triangle(p,ch->vertex(0)->point(),ch->vertex(2)->point(),ch->vertex(3)->point(),k));

////    Point bc = barycenter(ch);
////    // barycenter to triangle distances
////    double dc0 = sqrt(CGAL::internal::squared_distance_to_triangle(bc,ch->vertex(0)->point(),ch->vertex(1)->point(),ch->vertex(2)->point(),k));
////    double dc1 = sqrt(CGAL::internal::squared_distance_to_triangle(bc,ch->vertex(0)->point(),ch->vertex(1)->point(),ch->vertex(3)->point(),k));
////    double dc2 = sqrt(CGAL::internal::squared_distance_to_triangle(bc,ch->vertex(1)->point(),ch->vertex(2)->point(),ch->vertex(3)->point(),k));
////    double dc3 = sqrt(CGAL::internal::squared_distance_to_triangle(bc,ch->vertex(0)->point(),ch->vertex(2)->point(),ch->vertex(3)->point(),k));

////    // d(P,f_i)/d(c,f_i)
////    double d0 = dp0/dc0;
////    double d1 = dp1/dc1;
////    double d2 = dp2/dc2;
////    double d3 = dp3/dc3;

////    // min_i d(P,f_i)/d(c,f_i)
////    double mind = min(min(min(d0,d1),d2),d3);

////    // x
////    double x = (1.0+mind)/2.0;

////    if(0.5 > x || x > 1.0){
////        cout << x << endl;
////        cout << ch->vertex(0)->point() << endl;
////        cout << ch->vertex(1)->point() << endl;
////        cout << ch->vertex(2)->point() << endl;
////        cout << ch->vertex(3)->point() << endl;
////    }
////    assert(0.5<=x && x<=1.0);


////    // z_i
////    double z0 = max(0.0,1-d0)/2.0;
////    double z1 = max(0.0,1-d1)/2.0;
////    double z2 = max(0.0,1-d2)/2.0;
////    double z3 = max(0.0,1-d3)/2.0;
////    assert(0.0<=z0 && z0<=0.5);
////    assert(0.0<=z1 && z1<=0.5);
////    assert(0.0<=z2 && z2<=0.5);
////    assert(0.0<=z3 && z3<=0.5);

////    if(!gc){ // without graph cut
////        double occ0 = ch->neighbor(3)->info().outside_score - 0.5; // triangle 0,1,2
////        double occ1 = ch->neighbor(2)->info().outside_score - 0.5; // triangle 0,1,3
////        double occ2 = ch->neighbor(0)->info().outside_score - 0.5; // triangle 1,2,3
////        double occ3 = ch->neighbor(1)->info().outside_score - 0.5; // triangle 0,2,3
////        double occ = ch->info().outside_score - 0.5;
////        return (x*occ + z0*occ0 + z1*occ1 + z2*occ2 + z3*occ3) / (x+z0+z1+z2+z3);
////    }
////    else{ // with graph cut
////        bool stet = false;
////        // check if tet is surface tet
////        for(int i = 0; i < 4; i++){
////            if(ch->neighbor(i)->info().gc_label != ch->info().gc_label){
////                stet = true;
////                break;
////            }
////        }
////        if(stet){ // use continuous (smoothed) occupancy
////            double occ0 = ch->neighbor(3)->info().outside_score - 0.5; // triangle 0,1,2
////            double occ1 = ch->neighbor(2)->info().outside_score - 0.5; // triangle 0,1,3
////            double occ2 = ch->neighbor(0)->info().outside_score - 0.5; // triangle 1,2,3
////            double occ3 = ch->neighbor(1)->info().outside_score - 0.5; // triangle 0,2,3
////            double occ = ch->info().outside_score - 0.5;

//////            double occ0 = ch->neighbor(3)->info().gc_label - 0.5; // triangle 0,1,2
//////            double occ1 = ch->neighbor(2)->info().gc_label - 0.5; // triangle 0,1,3
//////            double occ2 = ch->neighbor(0)->info().gc_label - 0.5; // triangle 1,2,3
//////            double occ3 = ch->neighbor(1)->info().gc_label - 0.5; // triangle 0,2,3
//////            double occ = ch->info().gc_label - 0.5;
////            return (x*occ + z0*occ0 + z1*occ1 + z2*occ2 + z3*occ3) / (x+z0+z1+z2+z3);
////        }
////        else{ // use discrete occupancy (smoothing shouldn't change anything??)
////            return ch->info().gc_label - 0.5;
////        }
////    }
////}


////FT gc_occ_function (Point p){


////    // get the cell of the test point
////    Delaunay::Locate_type lt;
////    int li, lj;
////    auto ch = gDt.locate(p,lt,li,lj);
////    if(lt == Delaunay::OUTSIDE_CONVEX_HULL || lt == Delaunay::OUTSIDE_AFFINE_HULL){
////        const FT occ = 0.5;
////        return occ;
////    }
////    // check for neighboring labels
////    for(int i = 0; i < 4; i++){
////        if(ch->neighbor(i)->info().gc_label != ch->info().gc_label){
////            const FT occ = ch->info().outside_score - 0.5;
////            return occ;
////        }
////    }
////    const FT occ = ch->info().gc_label - 0.5;
////    return occ;


////}
////FT occ_function (Point p){

////    // get the cell of the test point and return its "signed occupancy"
////    Delaunay::Locate_type lt;
////    int li, lj;
////    auto ch = gDt.locate(p,lt,li,lj);
////    if(lt == Delaunay::OUTSIDE_CONVEX_HULL || lt == Delaunay::OUTSIDE_AFFINE_HULL){
////        const FT occ = 0.5;
////        return occ;
////    }
////    const FT occ = ch->info().outside_score - 0.5;
////    return occ;
////}
////FT gc_smoothed_occ_function (Point p){

////    // get the cell of the test point
////    Delaunay::Locate_type lt;
////    int li, lj;
////    auto ch = gDt.locate(p,lt,li,lj);
////    if(lt == Delaunay::OUTSIDE_CONVEX_HULL || lt == Delaunay::OUTSIDE_AFFINE_HULL){
////        const FT occ = 0.5;
////        return occ;
////    }
////    const FT occ = smoothedOccupancy(ch,p,1);
////    return occ;


////}
////FT smoothed_occ_function (Point p){
////    // get the cell of the test point and return its "signed occupancy"
////    Delaunay::Locate_type lt;
////    int li, lj;
////    auto ch = gDt.locate(p,lt,li,lj);
////    if(lt == Delaunay::OUTSIDE_CONVEX_HULL || lt == Delaunay::OUTSIDE_AFFINE_HULL){
////        const FT occ = 0.5;
////        return occ;
////    }
////    const FT occ = smoothedOccupancy(ch,p,0);
////    return occ;
////}



////int isoExtraction(dirHolder dir, dataHolder& data, runningOptions options){

////    auto start = std::chrono::high_resolution_clock::now();

////    std::cout << "\nExtract isosurface..." << endl;
////    cout << "\t-lower bound facet angle: " << options.iso_options[0] << endl;
////    cout << "\t-upper bound radius surface Delaunay balls: " << options.iso_options[1] << endl;
////    cout << "\t-upper bound for facet center-center distances: " << options.iso_options[2] << endl;



////    Tr tr;            // 3D-Delaunay triangulation
////    C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation

////    gDt = data.Dt;

////    int init = 20;

////    // defining the surface
////    if(options.optimization){
////        cout << "\t-use graph cut labels" << endl;
////        if(options.smooth_field){
////            cout << "\t-use smoothed field" << endl;
////            Surface_3 surface(gc_smoothed_occ_function,             // pointer to function
////                            Sphere_3(CGAL::ORIGIN, 100.0*100.0)); // bounding sphere with squared radius
////            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(options.iso_options[0],  // angular bound
////                                                           options.iso_options[1],  // radius bound
////                                                           options.iso_options[2]); // distance bound

////            cout << "\t-defined field" << endl;

////            // meshing surface
////            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(),init);
////        }
////        else{
////            cout << "\t-use unsmoothed field" << endl;
////            Surface_3 surface(gc_occ_function,             // pointer to function
////                            Sphere_3(CGAL::ORIGIN, 100.0*100.0)); // bounding sphere with squared radius
////            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(options.iso_options[0],  // angular bound
////                                                           options.iso_options[1],  // radius bound
////                                                           options.iso_options[2]); // distance bound

////            cout << "\t-defined field" << endl;

////            // meshing surface
////            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(),init);
////        }


////    }
////    else{
////        cout << "\t-use continuous labels" << endl;
////        if(options.smooth_field){
////            cout << "\t-use smoothed field" << endl;
////            Surface_3 surface(smoothed_occ_function,             // pointer to function
////                            Sphere_3(CGAL::ORIGIN, 100.0*100.0)); // bounding sphere with squared radius
////            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(options.iso_options[0],  // angular bound
////                                                           options.iso_options[1],  // radius bound
////                                                           options.iso_options[2]); // distance bound
////            cout << "\t-defined field" << endl;

////            // meshing surface
////            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(),init);
////        }
////        else{
////            cout << "\t-use unsmoothed field" << endl;
////            Surface_3 surface(occ_function,             // pointer to function
////                            Sphere_3(CGAL::ORIGIN, 100.0*100.0)); // bounding sphere with squared radius
////            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(options.iso_options[0],  // angular bound
////                                                           options.iso_options[1],  // radius bound
////                                                           options.iso_options[2]); // distance bound
////            cout << "\t-defined field" << endl;
////            // meshing surface
////            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(),init);
////        }
////    }

////    // criteria see here: https://doc.cgal.org/latest/Surface_mesher/classCGAL_1_1Surface__mesh__default__criteria__3.html

//////    a lower bound on the minimum angle in degrees of the surface mesh facets.
//////    an upper bound on the radius of surface Delaunay balls. A surface Delaunay ball is a ball circumscribing a facet, centered on the surface and empty of vertices. Such a ball exists for each facet of the current surface mesh. Indeed the current surface mesh is the Delaunay triangulation of the current sampling restricted to the surface which is just the set of facets in the three dimensional Delaunay triangulation of the sampling that have a Delaunay surface ball.
//////    an upper bound on the center-center distances of the surface mesh facets. The center-center distance of a surface mesh facet is the distance between the facet circumcenter and the center of its surface Delaunay ball.

////    //  topology	is the set of topological constraints which have to be verified by each surface facet. See section Delaunay Refinement for further details. Note that if one parameter is set to 0, then its corresponding criteria is ignored.


////    SurfaceMesh sm;
////    std::ofstream out(dir.path+dir.write_file+"_isosurface.off");


////    CGAL::output_surface_facets_to_off(out,c2t3);

////    cout << "\t-make surface mesh" << endl;
//////    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
//////    out << sm << std::endl;

//////    std::cout << "\t-number of points: " << tr.number_of_vertices() << "\n";
////    auto stop = chrono::high_resolution_clock::now();
////    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
////    cout << "\t-done after " << duration.count() << "s" << endl;
////}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////// CONTROL FUNCTIONS /////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//int surfaceReconstruction(dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& exportO){

//    // TOOD: add an assert for checking that there are no duplicate points in the input for the 3DT, because the ray tracing crashes
//    // if there are. Maybe it would however be best to change my data.points+data.infos vector to a CLGA::Point_set_3
//    // with property maps for each info.



//    ///////////////////////////////
//    ///////// FILE NAMING /////////
//    ///////////////////////////////
//    options.scoring="_"+options.scoring;
//    dir.rw_string = double2string(options.area_reg_weight+options.angle_reg_weight+options.cc_reg_weight+options.sv_reg_weight);
//    if(options.percentage_of_outliers > 0.0)
//        dir.ol_string = double2string(options.percentage_of_outliers*100);
//    if(dir.read_file.empty())
//        dir.read_file = dir.write_file;
//    if(dir.write_file.empty())
//        dir.write_file = dir.read_file+dir.ol_string+options.scoring+dir.rw_string;
//    else
//        dir.write_file+=options.scoring+dir.rw_string;

//    ///////////////////////////////
//    ///////// IMPORT DATA /////////
//    ///////////////////////////////
//    // ground truth
//    if(!dir.gt_poly_file.empty()){
//        if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "off"){
//            if(importOFFMesh(dir.path+dir.gt_poly_file, data.gt_poly))
//                return 1;
////            CGAL::Polygon_mesh_processing::keep_largest_connected_components(data.gt_poly, 1);
//        }
//        #ifdef RECONBENCH
//        else if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "mpu"){
//            if(importImplicit(dir, data))
//                return 1;
//        }
//        #endif
//        else{
//            cerr << "\nnot a valid ending for a ground truth file. choose either .mpu or .off" << endl;
//            return 1;
//        }
//        options.ground_truth = 1;
//    }
//    // tt input
////    if(options.data_source == "tt"){
////        // bounding box
////        cout << "\nRead bounding box from " << dir.path+dir.read_file+"_obb.ply" << endl;
////        std::ifstream in(dir.path+dir.read_file+"_obb.ply");
////        CGAL::read_ply(in,data.sm_obb);
////        if(data.sm_obb.is_empty()){
////            cout << "\nERROR: empty bounding box" << endl;
////            return 1;
////        }
//////        SurfaceMesh gt_mesh;
//////        CGAL::copy_face_graph(data.gt_poly, gt_mesh);
//////        CGAL::Polygon_mesh_processing::triangulate_faces(data.sm_obb);
//////        assert(CGAL::Polygon_mesh_processing::does_bound_a_volume(data.sm_obb));
//////        CGAL::Polygon_mesh_processing::clip(gt_mesh,data.sm_obb);
//////        data.gt_poly.clear();
//////        CGAL::copy_face_graph(gt_mesh,data.gt_poly);
////        // sensor locations
////        if(importTanksAndTemplesScannerLocations(dir, data, options))
////            return 1;
////    }
////    if(!options.gt_scan_source.empty()){
////        if(options.gt_scan_source == "eth"){
////            #ifdef PCL
////            if(importETH3D(dir,data)){
////                cout << "\nERROR: in importETH3D()" << endl;
////                return 1;
////            }
////            #endif
////            data.gt_points=data.points, data.points.clear();
////            data.gt_infos=data.infos, data.infos.clear();
////        }
////        else if(options.gt_scan_source == "lidar"){
////            dir.gt_scan_file;
////            // load the file from there;
////            data.gt_points=data.points, data.points.clear();
////            data.gt_infos=data.infos, data.infos.clear();
////        }
//////        else if(options.gt_scan_source == "tt"){
//////            if(importTanksAndTemples(dir,data,options)){
//////                cout << "\nERROR: in importTanksAndTemples()" << endl;
//////                return 1;
//////            }
//////            if(options.subsample_grid_spacing > 0)
//////                gridSimplify(data, options.subsample_grid_spacing);
//////            data.gt_points=data.points, data.points.clear();
//////            data.gt_infos=data.infos, data.infos.clear();
//////        }
////        else{
////            cout << "ERROR: not a valid ground truth scan source [--gts]" << endl;
////            return 1;
////        }
////    }
//    // reconstruction input
//    if(options.data_source == "scan" || options.data_source == "tt"){
//        if(!data.gt_poly.size_of_vertices() > 0){
//            cout << "\nERROR: you need to load a ground truth polygon (with -g) for scanning it" << endl;
//            return 1;
//        };
//        if(scanObjectClosed(data, options))
//            return 1;
//    }
//    else if(options.data_source == "lidar"){
//        importPLYPoints(dir, data);
//    }
//    else if(options.data_source == "npz"){
//        importNPZ(dir, data);
//    }
//    #ifdef COLMAP
//    else if(options.data_source == "colmap"){
//        readColmapFiles(dir, data);
//    }
//    #endif
//    #ifdef OpenMVS
//    else if(options.data_source == "omvs"){
//        if(loadOMVSScene(dir, data))
//            return 1;
//    }
//    #endif
////    else if(options.data_source == "tt"){
////        if(importTanksAndTemples(dir,data,options))
////            return 1;
////        options.ground_truth = 0;
////    }
////    #ifdef PCL
////    else if(options.data_source == "eth"){
////        if(importETH3D(dir,data)){
////            cout << "\nERROR: in importETH3D()" << endl;
////            return 1;
////        }
////    }
////    #endif
//    else{
//        cout << "ERROR: not a valid reconstruction input" << endl;
//        return 1;
//    }


//    ///////////////////////////////////
//    ///////// PREPROCESS DATA /////////
//    ///////////////////////////////////
//    if(!dir.transformation_file.empty()){
//        // import translation matrix
////        if(importTransformationMatrix(dir,data))
////            return 1;
//        if(applyTransformationMatrix(data))
//            return 1;
//    }
//    // calc gt obb
////    if(!options.gt_isclosed){
////        cout << "\nMake bounding box of open GT for cropping input..." << endl;
////        assert(options.ground_truth);
////        // get centroid of ground truth
////        std::array<Point, 8> obb_points;
////        CGAL::oriented_bounding_box(data.gt_poly,obb_points);
////        double mx = 0, my = 0, mz = 0;
////        for(const auto p : obb_points){
////            mx+=p.x();
////            my+=p.y();
////            mz+=p.z();
////        }
////        data.gt_centroid = Point(mx/8,my/8,mz/8+1000);
////        cout << "\t-ground truth centroid for ray tracing: " << data.gt_centroid << "-1000" << endl;
////        CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
////                              obb_points[4], obb_points[5], obb_points[6], obb_points[7], data.sm_obb);
////        CGAL::Polygon_mesh_processing::triangulate_faces(data.sm_obb);

////        clipWithBoundingBox(data.points,data.infos,data.sm_obb);

////    }

////    if(options.scale > 0.0)
////        scalePointSet(dir, data, options.scale);

//    if(options.scale > 0.0)
//        standardizePointSet(dir, data, options.scale);

////    if(options.scoring == "_rt")
////        orderSensors(data);
////    else
////        cout << "\nSensors not ordered for ray-tracing" << endl;

//    if(exportO.cameras){
//        dir.suffix = "_cameras";
//        exportCameraCenter(dir, data);
//    }
//    // export the scanned points
//    if(exportO.scan){
//        dir.suffix = "_scan";
//        exportPLY(dir, data.points, data.infos, exportO);
//    }


//    //////////////////////////////////////
//    /////// DELAUNAY TRIANGULATION ///////
//    //////////////////////////////////////
//    if(options.Dt_epsilon > 0.0)
//        makeAdaptiveDelaunayWithInfo(data, options.Dt_epsilon);
//    else
//        makeDelaunayWithInfo(data);

//    if(options.labatut_sigma == -1.0 && (options.scoring == "_rt" || options.scoring == "_clrt"))
//        options.labatut_sigma = pcaKNN(data.Dt);
//    else{
//        cout << "\nd parameter of Labatu set to " << options.labatut_sigma << endl;
//    }

////    cout << "Mean edge length after scaling: " << calcMeanEdgeLength(data) << endl;





////    // add camera points to delaunay
////    typedef CGAL::Labeled_mesh_domain_3<EPICK> Mesh_domain;
////    // Triangulation
////    typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
////    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
////    // Criteria
////    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

////    Mesh_domain domain(data.Dt);
////    Mesh_criteria criteria();
////    CGAL::refine_mesh_3(data.Dt, domain, criteria);




//    ///////////////////////////////////
//    /////// TETRAHEDRON SCORING ///////
//    ///////////////////////////////////
//    // make an index, necessary for graphExport e.g.
//    options.make_global_cell_idx=1;
//    indexDelaunay(data, options);
//    // put score on cells
//    if(options.scoring == "_cs"){
//        if(options.gt_isclosed){
//            if(labelObjectWithClosedGroundTruth(dir,data,options.number_of_points_per_cell))
//                return 1;
//        }
//        else{
//            if(labelObjectWithOpenGroundTruth(dir,data, options.number_of_points_per_cell))
//                return 1;
//        }
//        #ifdef RECONBENCH
//        else
//            labelObjectWithImplicit(data, options.number_of_points_per_cell);
//        #endif
//    }
//    else if(options.scoring == "_csrt"){
//        auto rc = processing::RayCaster(data, options);
//        if(options.gt_isclosed){
//            if(labelObjectWithClosedGroundTruth(dir,data, options.number_of_points_per_cell))
//                return 1;
//        }
//        else{
//            if(labelObjectWithOpenGroundTruth(dir,data, options.number_of_points_per_cell))
//                return 1;
//        }
//        rc.run(1);
//    }
//    else if(options.scoring == "_rt"){
//        // NOTE: careful! the processing rayTracing does not have the cell_distance_pairs function
//        // resulting in the fact that the ray cell edges will not be exported at all
//        if(options.score_type == "")
//            processing::rayTracing(data.Dt, options);
//        else{
//            auto rc = processing::RayCaster(data, options);
//            rc.run(1);
//        }
//    }
//    else if(options.scoring == "_lrt"){
//        // this constructs the features which are exported to the cell ray graph,
//        // important to note the difference between score, and features, where score is not used in any way in the learning.
//        if(options.scale == 0.0)
//            cout << "\nConsider turning on scaling with --sn, if your learning data is not yet scaled to a unit cube!" << endl;
//        learning::rayTracing(data.Dt, options);
//        // this is only for making a reconstruction of a lrt, but has no influence on the features or the learning
////        learning::aggregateScoreAndLabel(data.Dt);
//    }
//    else if(options.scoring == "_lrtcs"){
//        learning::rayTracing(data.Dt, options);
//        if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "off"){
//            if(options.gt_isclosed){
//                if(labelObjectWithClosedGroundTruth(dir,data, options.number_of_points_per_cell))
//                    return 1;
//            }
//            else{
//                if(labelObjectWithOpenGroundTruth(dir,data, options.number_of_points_per_cell))
//                    return 1;
//            }
//        }
//        #ifdef RECONBENCH
//        if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "mpu"){
//            labelObjectWithImplicit(data, options.number_of_points_per_cell);
//        }
//        #endif

//    }
//    else if(options.scoring == "_cl"){
//        if(loadPrediction(dir, data, options))
//            return 1;
//    }
//    else if(options.scoring == "_clrt"){
//        auto rc = processing::RayCaster(data, options);
//        if(loadPrediction(dir, data, options))
//            return 1;
//        rc.run(1);
//    }
//    else if(options.scoring == "_clcs"){
//        if(options.gt_isclosed){
//            if(labelObjectWithClosedGroundTruth(dir,data, options.number_of_points_per_cell))
//                return 1;
//        }
//        else{
//            if(labelObjectWithOpenGroundTruth(dir,data, options.number_of_points_per_cell))
//                return 1;
//        }
//        if(loadPrediction(dir, data, options))
//            return 1;
//        // TODO: tweak the graphCut function to use unaries from cell.prediction_inside/outside but sv - binaries from cell.inside/outside_score

//    }



//    ////////////////////////////
//    /////// OPTIMIZATION ///////
//    ////////////////////////////
//    if(options.optimization == 1){
//        options.gco_iterations = -1;
//        if(options.score_type == "mine" || options.scoring == "_cl" || options.scoring == "_cs" || options.scoring == "_clcs" || options.scoring == "_lrtcs")
//            graphCutTet(data, options);
//        else
//            graphCutFacet(data,options);
//    }
//    else{
//        cout << "\nLabel cells without optimization by taking max score" << endl;
//        cout << "\t-infinite cells will be labelled outside" << endl;
//        for(auto cit = data.Dt.all_cells_begin(); cit != data.Dt.all_cells_end(); cit++){
//            if(data.Dt.is_infinite(cit)){
//                cit->info().gc_label = 1;
//                continue;
//            }
//            // this means untraversed cells and 50/50 cells will be labelled as inside
//            cit->info().gc_label = cit->info().outside_score > cit->info().inside_score ? 1 : 0;
//        }
//    }

//    ////////////////////////////
//    ////////// EXPORT //////////
//    ////////////////////////////


//    //// export features and labels
////    if(options.scoring == "_lrtcs" || options.scoring == "_lrt"){
////        // check if labels directory exists, if not create it
////        path lpath(dir.path);
////        // this shitty problem here is not present on the laptop
////        lpath /= string("gt");
////        if(!is_directory(lpath))
////            create_directory(lpath);
////        // export the features and labels
////        exportGraph(dir, options, data.Dt);
////    }



//    // make new 3DT with old 3DT cell centers, attribute cell score to vertices


////    #ifdef Open3D
//////    Delaunay Dt;
//////    Vertex_handle vh;
//////    vertex_info vi;
//////    Point p1,p2,p3,p4,centroid;
//////    for(auto cit = data.Dt.finite_cells_begin(); cit != data.Dt.finite_cells_end(); cit++){
//////        p1 = cit->vertex(0)->point();
//////        p2 = cit->vertex(1)->point();
//////        p3 = cit->vertex(2)->point();
//////        p4 = cit->vertex(3)->point();
//////        centroid = CGAL::centroid(p1,p2,p3,p4);
//////        vh = Dt.insert(centroid);
//////        vi.occupancy = cit->info().outside_score;
//////        vh->info() = vi;
//////    }
//////    exportIsoSurface(dir, Dt, 0);
//////    exportConvexHull(dir, Dt);

////    if(exportO.isosurface)
////        exportMITSurface(dir, data.Dt, 1, 0.5);

////    #endif



//    //// export interface
//    if(exportO.interface)
//        exportInterface(dir, data, options, exportO);

//    #ifdef OpenMVS
//    if(options.clean_mesh)
//        omvsCleanMesh(dir,data);
//    #endif

////    if(exportO.coloredFacets)
////        exportColoredFacets(dir, data.Dt, options.optimization);
//    if(exportO.cellScore){
//        exportCellCenter(dir, data.Dt);
//        exportCellScore(dir, data.Dt);
//    }
//    if(exportO.convexHull){
//        exportConvexHull(dir, data.Dt);
//    }
//    ////////// create surface mesh ////////////
//    // set export options
//    meshProcessingOptions mpOptions;
//    mpOptions.try_to_close = options.try_to_close;
//    mpOptions.try_to_make_manifold = options.try_to_make_manifold;
//    mpOptions.number_of_components_to_keep = options.number_of_components_to_keep;
//    mpOptions.factor_for_removing_large_faces = options.factor_for_removing_large_faces;
//    if(exportO.mesh){
//        dir.suffix = "_mesh";
//        createSurfaceMesh(data, mpOptions);
//        exportPLY(dir, data.smesh);
//    }


//    ////////// sample surface mesh //////////
//    if(exportO.sampling){
//        data.points.clear();
//        data.infos.clear();
//        if(data.smesh.is_empty()){
//            createSurfaceMesh(data, mpOptions);
//        }
//        sampleMesh(data, exportO);
//        dir.suffix = "_sampled";
//        // turn of color and sensor_vec export, because mesh sampling obviously does not have that
//        exportO.sensor_vec = false;
//        exportO.sensor_position = false;
//        exportO.color = false;
//        exportO.normals = false;
//        exportPLY(dir, data.points, data.infos, exportO);
//    }

//    if(options.evaluate_mesh){
//        if(data.smesh.is_empty()){
//            createSurfaceMesh(data, mpOptions);
//        }
//        if(options.ground_truth){
//            double iou;
//            if(calcIOU(data, exportO.sampling_method_option, iou))
//                return 1;
//            printMeshEvaluation(dir,data,iou);
//        }
//        else
//            printMeshEvaluation(dir,data,-1.0);

//    }

//    return 0;
//}



//int main(int argc, char const *argv[]){


//    cliParser ip("clf");
//    if(ip.parse(argc, argv))
//        return 1;
//    if(ip.getInput())
//        return 1;
//    if(ip.getOutput())
//        return 1;
//    if(ip.getLabatut())
//        return 1;

//    auto start = std::chrono::high_resolution_clock::now();
//    cout << "\n-----DELAUNAY-GRAPH-CUT-BASED SURFACE RECONSTRUCTION-----" << endl;
//    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

//    dataHolder data;
//    if(surfaceReconstruction(ip.dh, data, ip.ro, ip.eo))
//        return 1;

//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    cout << "\n-----DELAUNAY-GRAPH-CUT-BASED SURFACE RECONSTRUCTION FINISHED in "<< duration.count() << "s -----\n" << endl;

//    return 0;

//}



