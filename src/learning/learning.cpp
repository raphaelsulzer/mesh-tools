#include <random>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
//#include <util/helper.h>
//#include <util/geometricOperations.h>

#include <processing/meshProcessing.h>
#include <processing/pointSetProcessing.h>

#include <learning/learning.h>
#include <learning/learningMath.h>
#include <learning/learningRayTracing.h>
#include <learning/learningIO.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <boost/foreach.hpp>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Side_of_triangle_mesh.h>

using namespace std;

int labelReal(dirHolder dir, dataHolder& data, int sampling_points){

    auto start = std::chrono::high_resolution_clock::now();

    if(data.gt_poly.size_of_vertices() == 0){
        cout << "\nERROR: you need to load a ground truth polygon (with -g) for scanning it" << endl;
        return 1;
    }

    assert(sampling_points > 0);

    cout << "\nMake bounding box of open GT for cropping input..." << endl;
    // get centroid of ground truth
    getOrientedBoundingBox(data.gt_poly,data.bb_array,data.bb_surface_mesh);
    data.gt_centroid = getBoundingBoxCentroid(data.bb_array,1000.0);
    cout << "\t-ground truth centroid for ray tracing: " << data.gt_centroid << "-1000" << endl;
    clipWithBoundingBox(data.points,data.infos,data.bb_surface_mesh);


    cout << "\nLabelling real with open ground truth and ray tracing..." << endl;
    cout << "\t-sample " << sampling_points << endl;

    // Initialize the point-in-polyhedron tester
    // Construct AABB tree with a KdTree

    SurfaceMesh smesh;
    CGAL::copy_face_graph(data.gt_poly, smesh);
    CGAL::Polygon_mesh_processing::remove_connected_components_of_negligible_size(smesh);
//    dir.suffix = "_repaired_mesh";
//    exportOFF(dir,smesh);
//    cout << "gt centroid: " << data.gt_centroid << endl;

//    Tree tree(faces(smesh).first, faces(smesh).second, smesh);

    vector<face_descriptor> fds;
    for(face_descriptor fd : faces(smesh)){
        if(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(fd,smesh))
            fds.push_back(fd);

    }
    Tree tree(fds.begin(), fds.end(), smesh);
    tree.accelerate_distance_queries();

//    Polyhedron poly;
//    CGAL::copy_face_graph(data.sm_obb, poly);
//    Tree bb_tree(faces(data.sm_obb).first, faces(data.sm_obb).second, data.sm_obb);
//    const Point_inside inside_bb(data.sm_obb);

    int icount = 0;
    int ncount = 0;
    Delaunay::All_cells_iterator aci;
    // random sampler
    CGAL::Random random(42);
    vector<Point> sampled_points;
    double iperc;
    double operc;
    int newCellIndex = 0;
//    vector<Point> all_sampled_points;
//    vector<vertex_info> all_sampled_infos;
//    int progress = 0;
//    int total = data.Dt.number_of_cells();
    for(aci = data.Dt.all_cells_begin(); aci != data.Dt.all_cells_end(); aci++){
//        cout << progress++ << endl;
//        if((progress*100/total) % 10 == 0 && (++progress*100/total) % 10 != 0)
//            cout << progress << "/" << total << endl;

        aci->info().global_idx=newCellIndex++;

        if(data.Dt.is_infinite(aci)){
            aci->info().gc_label = 1; // more likely that it is outside
            aci->info().inside_score = 0.0;
            aci->info().outside_score = 1.0;
            continue;
        }

        // sample points inside the tet and check if they are inside or outside
        double cinside=0;
        double coutside=0;
        Tetrahedron current_tet = data.Dt.tetrahedron(aci);
        CGAL::Random_points_in_tetrahedron_3<Point> tet_point_sampler(current_tet, random);
        sampled_points.clear();
        CGAL::cpp11::copy_n(tet_point_sampler, sampling_points, std::back_inserter(sampled_points));
        for(auto const& sampled_point : sampled_points){


//            all_sampled_points.push_back(sampled_point);
//            vertex_info vi;
//            if(inside_bb(sampled_point) == CGAL::ON_UNBOUNDED_SIDE){
//                vi.color = CGAL::violet();
//                all_sampled_infos.push_back(vi);
//                cinside+=1.0;
//            }
//            else{
            Segment seg(sampled_point,data.gt_centroid);
            int n = tree.number_of_intersected_primitives(seg);
//            cout << "number of intersected primitives: " << n << endl;
            if(n % 2){
                cinside+=1.0;
//                vi.color = CGAL::red();
//                all_sampled_infos.push_back(vi);
            }
            else{
                coutside+=1.0;
//                vi.color = CGAL::blue();
//                all_sampled_infos.push_back(vi);
            }
//            if(n == 0)
//                vi.color = CGAL::yellow();
//            else if(n == 1)
//                vi.color = CGAL::green();
//            else if(n == 2)
//                vi.color = CGAL::blue();
//            else if(n == 3)
//                vi.color = CGAL::red();
//            else
//                vi.color = CGAL::black();
//            all_sampled_infos.push_back(vi);


//            }
        }
        iperc = cinside/(cinside+coutside);
        operc = coutside/(cinside+coutside);
        aci->info().outside_score = operc;
        aci->info().inside_score = iperc;
//        aci->info().gt_outside = operc;
//        aci->info().gt_inside = iperc;
        if(operc>iperc){
            ncount++;
//            aci->info().gc_label = 1;
        }
        else{
            icount++;
//            aci->info().gc_label = 0;
        }

    }
//    dir.suffix = "_sampled";
//    exportOptions eo;
//    eo.color = true;
//    exportPLY(dir,all_sampled_points,all_sampled_infos,eo);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Labelled object with " << icount << "/" << ncount << "(inside/outside) = " << icount+ncount << " cells" << endl;
    cout << "\t-in "<<duration.count() << "s" << endl;
    if(icount == 0 || ncount == 0){
        cout << "\nERROR: cells could not be labelled" << endl;
        return 1;
    }
    return 0;
}



int labelSRD(dirHolder dir, dataHolder& data, int sampling_points){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nLabelling SRD by sampling points inside cells..." << endl;
    cout << "\t-sample " << sampling_points << endl;

    if(!data.gt_poly.is_closed()){
        cout << "\nERROR: ground truth is not closed" << endl;
        return 1;
    }
    assert(sampling_points > 0);

    getOrientedBoundingBox(data.gt_poly,data.bb_array,data.bb_surface_mesh);

    // Initialize the point-in-polyhedron tester
    // Construct AABB tree with a KdTree
    SurfaceMesh smesh;
    CGAL::copy_face_graph(data.gt_poly, smesh);
    Tree tree(faces(smesh).first, faces(smesh).second, smesh);
    tree.accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    const Point_inside inside_mesh(tree);
    const Point_inside inside_bbox(data.bb_surface_mesh);


    int icount = 0;
    int ncount = 0;
    Delaunay::All_cells_iterator aci;
    // random sampler
    CGAL::Random random(42);
    vector<Point> sampled_points;
    double iperc;
    double operc;
    int newCellIndex = 0;
    vector<Point> all_sampled_points;
    vector<vertex_info> all_sampled_infos;
    for(aci = data.Dt.all_cells_begin(); aci != data.Dt.all_cells_end(); aci++){

        aci->info().global_idx=newCellIndex++;

        if(data.Dt.is_infinite(aci)){
            aci->info().gc_label = 1; // more likely that it is outside
            aci->info().inside_score = 0.0;
            aci->info().outside_score = 1.0;
            continue;
        }

        // sample points inside the tet and check if they are inside or outside
        double cinside=0;
        double coutside=0;
        Tetrahedron current_tet = data.Dt.tetrahedron(aci);
        CGAL::Random_points_in_tetrahedron_3<Point> tet_point_sampler(current_tet, random);
        sampled_points.clear();
        CGAL::cpp11::copy_n(tet_point_sampler, sampling_points, std::back_inserter(sampled_points));
        for(auto const& sampled_point : sampled_points){
            all_sampled_points.push_back(sampled_point);
            vertex_info vi;

            if(sampled_point.y() > data.bb_array[0].y()){
                coutside+=1.0;
//                vi.color = CGAL::blue();
            }
            else if(inside_bbox(sampled_point) == CGAL::ON_UNBOUNDED_SIDE){
                cinside+=1.0;
//                vi.color = CGAL::red();
            }
            else if(inside_mesh(sampled_point) == CGAL::ON_BOUNDED_SIDE){
                cinside+=1.0;
//                vi.color = CGAL::red();
            }
            else{
                coutside+=1.0;
//                vi.color = CGAL::blue();
            }
//            all_sampled_infos.push_back(vi);

        }
        iperc = cinside/(cinside+coutside);
        operc = coutside/(cinside+coutside);
        aci->info().outside_score = operc;
        aci->info().inside_score = iperc;

        if(operc>iperc){
            ncount++;
//            aci->info().gc_label = 1;
        }
        else{
            icount++;
//            aci->info().gc_label = 0;
        }
    }
//    dir.suffix = "_sampled";
//    exportOptions eo;
//    eo.color = true;
//    exportPLY(dir,all_sampled_points,all_sampled_infos,eo);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Labelled object with " << icount << "/" << ncount << "(inside/outside) = " << icount+ncount << " cells" << endl;
    cout << "\t-in "<<duration.count() << "s" << endl;
    if(icount == 0 || ncount == 0){
        cout << "\nERROR: cells could not be labelled" << endl;
        return 1;
    }
    return 0;
}


int labelObject(dirHolder dir, dataHolder& data, int sampling_points){

    auto start = std::chrono::high_resolution_clock::now();

    if(!data.gt_poly.is_closed()){
        cout << "\nERROR: ground truth is not closed" << endl;
        return 1;
    }

    assert(sampling_points > 0);

    cout << "\nLabelling object by sampling points inside cells..." << endl;
    cout << "\t-sample " << sampling_points << endl;

    // Initialize the point-in-polyhedron tester
    // Construct AABB tree with a KdTree
    SurfaceMesh smesh;
    CGAL::copy_face_graph(data.gt_poly, smesh);
    Tree tree(faces(smesh).first, faces(smesh).second, smesh);
//    Tree tree(faces(data.gt_poly).first, faces(data.gt_poly).second, data.gt_poly);
    tree.accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    const Point_inside inside_tester(tree);
//    Polyhedron poly;
//    CGAL::copy_face_graph(data.sm_obb, poly);
////    Tree bb_tree(faces(data.sm_obb).first, faces(data.sm_obb).second, data.sm_obb);
//    const Point_inside inside_bb(poly);
    
    int icount = 0;
    int ncount = 0;
    Delaunay::All_cells_iterator aci;
    // random sampler
    CGAL::Random random(42);
    vector<Point> sampled_points;
    double iperc;
    double operc;
    int newCellIndex = 0;
    vector<Point> all_sampled_points;
    vector<vertex_info> all_sampled_infos;
    for(aci = data.Dt.all_cells_begin(); aci != data.Dt.all_cells_end(); aci++){

        aci->info().global_idx=newCellIndex++;

        if(data.Dt.is_infinite(aci)){
            aci->info().gc_label = 1; // more likely that it is outside
            aci->info().inside_score = 0.0;
            aci->info().outside_score = 1.0;
            continue;
        }

        // sample points inside the tet and check if they are inside or outside
        double cinside=0;
        double coutside=0;
        Tetrahedron current_tet = data.Dt.tetrahedron(aci);
        CGAL::Random_points_in_tetrahedron_3<Point> tet_point_sampler(current_tet, random);
        sampled_points.clear();
        CGAL::cpp11::copy_n(tet_point_sampler, sampling_points, std::back_inserter(sampled_points));
        for(auto const& sampled_point : sampled_points){
//            all_sampled_points.push_back(sampled_point);
//            vertex_info vi;

////            cout<<"\t check point "<< pcount++ <<"/"<<sampling_points<<endl;
//            if(inside_bb(sampled_point) != CGAL::ON_BOUNDED_SIDE){
//                vi.color = CGAL::green();
//                all_sampled_infos.push_back(vi);
//                coutside+=1.0;
//                continue;
//            }

            if(inside_tester(sampled_point) == CGAL::ON_BOUNDED_SIDE){
                cinside+=1.0;
//                vi.color = CGAL::red();
//                all_sampled_infos.push_back(vi);
            }
            else{
                coutside+=1.0;
//                vi.color = CGAL::blue();
//                all_sampled_infos.push_back(vi);
            }
        }
        iperc = cinside/(cinside+coutside);
        operc = coutside/(cinside+coutside);
        aci->info().outside_score = operc;
        aci->info().inside_score = iperc;

        if(operc>iperc){
            ncount++;
//            aci->info().gc_label = 1;
        }
        else{
            icount++;
//            aci->info().gc_label = 0;
        }

    }
//    dir.suffix = "_sampled";
//    exportOptions eo;
//    eo.color = true;
//    exportPLY(dir,all_sampled_points,all_sampled_infos,eo);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Labelled object with " << icount << "/" << ncount << "(inside/outside) = " << icount+ncount << " cells" << endl;
    cout << "\t-in "<<duration.count() << "s" << endl;
    if(icount == 0 || ncount == 0){
        cout << "\nERROR: cells could not be labelled" << endl;
        return 1;
    }
    return 0;
}


#ifdef RECONBENCH
int labelObjectWithImplicit(dataHolder& data, int sampling_points){

    auto start = std::chrono::high_resolution_clock::now();

    assert(sampling_points > 0);

    cout << "\nLabelling cells by sampling points inside cells: " << endl;
    cout << "\t-sample " << sampling_points << endl;


    int icount = 0;
    int ncount = 0;
    Delaunay::All_cells_iterator aci;
    // random sampler
    CGAL::Random random(42);
    vector<Point> sampled_points;
    double iperc;
    double operc;
    int newCellIndex = 0;
    for(aci = data.Dt.all_cells_begin(); aci != data.Dt.all_cells_end(); aci++){

        if(data.Dt.is_infinite(aci)){
            aci->info().gc_label = 1; // more likely that it is outside
            aci->info().inside_score = 0.2;
            aci->info().outside_score = 0.8;
            aci->info().global_idx=newCellIndex++;
            continue;
        }
        aci->info().global_idx=newCellIndex++;

        // sample points inside the tet and check if they are inside or outside
        double cinside=0;
        double coutside=0;
        Tetrahedron current_tet = data.Dt.tetrahedron(aci);
        CGAL::Random_points_in_tetrahedron_3<Point> tet_point_sampler(current_tet, random);
        sampled_points.clear();
        CGAL::cpp11::copy_n(tet_point_sampler, sampling_points, std::back_inserter(sampled_points));
        for(auto const& sampled_point : sampled_points){


            auto sampled_point_rb = point2Vector3RB(sampled_point);

            if(data.implicit->function(sampled_point_rb) <= 0.0){
                cinside+=1.0;
//                Color col = CGAL::make_array((unsigned char)255, (unsigned char)0, (unsigned char)0);
//                vertex_info vi;
//                vi.color = col;
//                all_sampled_infos.push_back(vi);
            }
            else{
                coutside+=1.0;
//                Color col = CGAL::make_array((unsigned char)0, (unsigned char)0, (unsigned char)255);
//                vertex_info vi;
//                vi.color = col;
//                all_sampled_infos.push_back(vi);
            }
        }
        iperc = cinside/(cinside+coutside);
        operc = coutside/(cinside+coutside);
        aci->info().outside_score = operc;
        aci->info().inside_score = iperc;
        if(operc>iperc){
            ncount++;
            aci->info().gc_label = 1;
        }
        else{
            icount++;
            aci->info().gc_label = 0;
        }

    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Labelled object with " << icount << "/" << ncount << "(inside/outside) = " << icount+ncount << " cells" << endl;
    cout << "\t-in "<<duration.count() << "s" << endl;
    if(icount == 0 || ncount == 0)
        return 0;
    return 1;
}
#endif
