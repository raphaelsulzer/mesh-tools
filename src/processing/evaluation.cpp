#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

using namespace std;

#include <CGAL/Polygon_mesh_processing/distance.h>
int sampleMesh(dataHolder& data, exportOptions exportO){

    cout << "\nSample surface mesh..." << endl;


    if(exportO.sampling_method == "gs"){
        CGAL::Polygon_mesh_processing::sample_triangle_mesh(data.smesh, std::back_inserter(data.points),
                                                            CGAL::parameters::use_grid_sampling(true)
                                                            .grid_spacing(exportO.sampling_method_option)
                                                            .do_sample_vertices(false)
                                                            .do_sample_edges(false)
                                                            );
    }
    else if(exportO.sampling_method == "ps"){
        CGAL::Polygon_mesh_processing::sample_triangle_mesh(data.smesh, std::back_inserter(data.points),
                                                            CGAL::parameters::use_random_uniform_sampling(true)
                                                            .number_of_points_on_faces(exportO.sampling_method_option)
                                                            .do_sample_vertices(false)
                                                            .do_sample_edges(false)
                                                            );
    }
    else if(exportO.sampling_method == "as"){
        double points_to_sample = (int)(CGAL::Polygon_mesh_processing::area(data.smesh)*exportO.sampling_method_option);
        CGAL::Polygon_mesh_processing::sample_triangle_mesh(data.smesh, std::back_inserter(data.points),
                                                            CGAL::parameters::use_random_uniform_sampling(true)
                                                            .number_of_points_on_faces(points_to_sample)
                                                            .do_sample_vertices(false)
                                                            .do_sample_edges(false)
                                                            );
    }
    else{
        cout << "ERROR: not a valid sampling options" << endl;
        return 1;
    }

    cout << "\t-sampled " << data.points.size() << endl;
    if(exportO.sampling_method == "gs")
        cout << "\t-with grid spacing " << exportO.sampling_method_option << endl;
    if(exportO.sampling_method == "ps")
        cout << "\t-with random uniform sampling " << endl;
    if(exportO.sampling_method == "as")
        cout << "\t-with random uniform sampling of " << exportO.sampling_method_option << " per area unit" << endl;



};


//// how many cells per edge
//vector<int> n_edges;
//int count;
//for(auto fei = data.Dt.finite_edges_begin(); fei != data.Dt.finite_edges_end(); fei++){
//    auto inc_cells = data.Dt.incident_cells(*fei);
//    auto first_cell = inc_cells;
//    count = 1;
//    while(first_cell != ++inc_cells)
//        count++;
//    n_edges.push_back((double)count);
//}

//double sum = 0.0;
//for(auto n : n_edges)
//    sum+=n;
//double mean = sum / (double)n_edges.size();
//sum = 0.0;
//for(auto n : n_edges)
//    sum+=sqrt((n-mean)*(n-mean));
//double var = sum / (double)n_edges.size();
//cout << "n_edes" << n_edges.size() << endl;
//cout << "mean " << mean << endl;
//cout << "var " << sqrt(var) << endl;

//// edge length
double calcMeanEdgeLength(dataHolder& data){

    double sum = 0;
    for(auto fei = data.Dt.finite_edges_begin(); fei != data.Dt.finite_edges_end(); fei++){
        int i = fei->second;
        int j = fei->third;
        auto cell = fei->first;

        auto p1 = cell->vertex(i)->point();
        auto p2 = cell->vertex(j)->point();

        sum+=sqrt(CGAL::squared_distance(p1,p2));
    }
    return sum / data.Dt.number_of_edges();
}





#include <CGAL/Polygon_mesh_processing/distance.h>
void checkMeshQuality(vector<Point>& ori_points, Polyhedron& surface_mesh){

    double max_dist_to_point_set =
          CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(surface_mesh, ori_points, 4000);
    cout << "Max distance to point set (precision): " << max_dist_to_point_set << endl;
    double max_dist_to_triangle_mesh =
          CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(ori_points, surface_mesh);
    cout << "Max distance to tiangle mesh (recall): " << max_dist_to_triangle_mesh << endl;
}

#include <CGAL/Side_of_triangle_mesh.h>
typedef CGAL::Side_of_triangle_mesh<Polyhedron, EPICK> Point_inside;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
double calcIOU(Polyhedron& gt, SurfaceMesh& reconm, int test_points, double& iou){

    cout << "\nCalculate intersection over union..." << endl;

    Polyhedron recon;
    CGAL::copy_face_graph(reconm, recon);

    if(test_points < 1){
        cerr << "Specify number of test points with '--output_sampling ps,X'\n" << endl;
        return 1;
    }
    else{
        cout << "\t-sample " << test_points << " inside the bounding box" << endl;
    }

    if(!gt.is_closed()){
        cerr << "Ground truth mesh is not closed" << endl;
        return 1;
    }
    if(!recon.is_closed()){
        cout << "\nWARNING: Reconstructed mesh is not closed (, probably due to dublication of non-manifold edges)." << endl;
        cout << "Point-in-polyhedron test used here still allows to give a reasonable number for IoU." << endl;
        cout << "TODO: would be to leave the graph-cut surface non-manifold (, which however the SurfaceMesh or " << endl;
        cout << "Polyhedron CGAL datastructures do not allow). Then a completeley valid IoU could be calculated" << endl;
//        iou = 0.0;
//        return 0;
    }




    // get all points of gt and reconstuction
    vector<Point> all_points;
    Polyhedron::Vertex_iterator vit;
    for(vit = recon.vertices_begin(); vit != recon.vertices_end(); vit++){
        all_points.push_back(vit->point());
    }
    for(vit = gt.vertices_begin(); vit != gt.vertices_end(); vit++){
        all_points.push_back(vit->point());
    }

    Iso_cuboid bb = CGAL::bounding_box(all_points.begin(), all_points.end());
    Point centroid = CGAL::centroid(all_points.begin(), all_points.end());
    double max_edge_length = max(max(bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin()), bb.zmax()-bb.zmin());
    CGAL::Random_points_in_cube_3<Point>gen(max_edge_length);

    // Initialize the point-in-polyhedron tester
    Tree gt_tree(faces(gt).first, faces(gt).second, gt);
    gt_tree.accelerate_distance_queries();
    const Point_inside gt_inside_tester(gt_tree);

    // Initialize the point-in-polyhedron tester
    Tree recon_tree(faces(recon).first, faces(recon).second, recon);
    recon_tree.accelerate_distance_queries();
    const Point_inside recon_inside_tester(recon_tree);

    int successfully_tested_points = 0;
    int intersection_points = 0;
    int union_points = 0;
    while(successfully_tested_points < test_points){

        Point rand_point = *gen++;
        rand_point = Point(rand_point.x()+centroid.x(), rand_point.y()+centroid.y(), rand_point.z()+centroid.z());

        if(gt_inside_tester(rand_point) == CGAL::ON_BOUNDED_SIDE &&
                recon_inside_tester(rand_point) == CGAL::ON_BOUNDED_SIDE){
            intersection_points++;
            union_points++;
            successfully_tested_points++;
        }
        else if(gt_inside_tester(rand_point) == CGAL::ON_BOUNDED_SIDE &&
             recon_inside_tester(rand_point) != CGAL::ON_BOUNDED_SIDE){
            union_points++;
            successfully_tested_points++;
        }
        else if(gt_inside_tester(rand_point) != CGAL::ON_BOUNDED_SIDE &&
             recon_inside_tester(rand_point) == CGAL::ON_BOUNDED_SIDE){
            union_points++;
            successfully_tested_points++;
        }

    }

    iou = double(intersection_points)*100.0/double(test_points);
    return 0;


}


#include <CGAL/Side_of_triangle_mesh.h>
typedef CGAL::Side_of_triangle_mesh<Polyhedron, EPICK> Point_inside;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> PolyhedronTree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
double calcIOU(dataHolder& data, int test_points, double& iou){

    cout << "\nCalculate intersection over union..." << endl;

    if(test_points < 1){
        cerr << "Specify number of test points with '--output_sampling ps,X'\n" << endl;
        return 1;
    }
    else{
        cout << "\t-sample " << test_points << " inside the bounding box" << endl;
    }

    CGAL::copy_face_graph(data.smesh, data.reconstructed_poly);

    if(!data.gt_poly.is_closed()){
        cerr << "Ground truth mesh is not closed" << endl;
        return 1;
    }


    // get all points of gt and reconstuction
    vector<Point> all_points;
    for(auto pit = data.remaining_points.begin(); pit != data.remaining_points.end(); pit++){
        all_points.push_back(*pit);
    }
    for(auto vit = data.gt_poly.vertices_begin(); vit != data.gt_poly.vertices_end(); vit++){
        all_points.push_back(vit->point());
    }

    Iso_cuboid bb = CGAL::bounding_box(all_points.begin(), all_points.end());
    Point centroid = CGAL::centroid(all_points.begin(), all_points.end());
    double max_edge_length = max(max(bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin()), bb.zmax()-bb.zmin());
    CGAL::Random_points_in_cube_3<Point>gen(max_edge_length);

    // Initialize the point-in-polyhedron tester
//    SurfaceMesh smesh;
//    CGAL::copy_face_graph(data.gt_poly, smesh);
////    vector<face_descriptor> fds;
////    for(face_descriptor fd : faces(smesh)){
////        if(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(fd,smesh))
////            fds.push_back(fd);

////    }
//    Tree gt_tree(faces(smesh).first, faces(smesh).second, smesh);
    PolyhedronTree gt_tree(faces(data.gt_poly).first, faces(data.gt_poly).second, data.gt_poly);
    gt_tree.accelerate_distance_queries();
    const Point_inside gt_inside_tester(gt_tree);

    int successfully_tested_points = 0;
    int intersection_points = 0;
    int union_points = 0;
    while(successfully_tested_points < test_points){

        Point rand_point = *gen++;
        rand_point = Point(rand_point.x()+centroid.x(), rand_point.y()+centroid.y(), rand_point.z()+centroid.z());

        bool gt_inside = false;
        if(gt_inside_tester(rand_point) == CGAL::ON_BOUNDED_SIDE)
            gt_inside = true;

        bool recon_inside = false;
        if(data.Dt.locate(rand_point)->info().gc_label == 0.0)
            recon_inside = true;

        if(gt_inside && recon_inside){
            intersection_points++;
            union_points++;
            successfully_tested_points++;
        }
        else if(gt_inside && !recon_inside){
            union_points++;
            successfully_tested_points++;
        }
        else if(!gt_inside && recon_inside){
            union_points++;
            successfully_tested_points++;
        }

    }

//    cout << "\t-(intersection / union) = (" << double(intersection_points)/double(test_points) << " / " << double(union_points)/double(test_points) << ")" << endl;

    iou = double(intersection_points)*100.0/double(test_points);
    return 0;

}





int printMeshEvaluation(dirHolder dir, dataHolder& data, double iou)
{

    typedef boost::graph_traits<SurfaceMesh>::face_descriptor face_descriptor;
    SurfaceMesh::Property_map<face_descriptor, std::size_t> fccmap =  data.smesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
    size_t components = CGAL::Polygon_mesh_processing::connected_components(data.smesh, fccmap);

    double area = CGAL::Polygon_mesh_processing::area(data.smesh);

    cout << "\nMesh evaluation..." << endl;
    cout << "\t-of file " <<  dir.write_file+dir.suffix << endl;
    cout << "\t-components: " << components << endl;
    cout << "\t-vertices: " << data.smesh.number_of_vertices() << endl;
    cout << "\t-edges: " << data.smesh.number_of_edges() << endl;
    cout << "\t-faces: " << data.smesh.number_of_faces() << endl;
    cout << "\t-area: " << area << endl;
    cout << "\t-closed: " << CGAL::is_closed(data.smesh) << endl;
    cout << "\t-sintersecting: " << PMP::does_self_intersect(data.smesh) << endl;
    if(iou > -1)
        cout << "\t-iou: " << iou << endl;

    return 0;
}



