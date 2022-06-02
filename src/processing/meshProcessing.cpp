#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <processing/pointSetProcessing.h>
#include <processing/meshProcessing.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/polygon_mesh_processing.h>


namespace VSA = CGAL::Surface_mesh_approximation;
namespace PMP = CGAL::Polygon_mesh_processing;
using namespace std;


////////////////////////////////////////////////////////////
////////////////////// Triangulation ///////////////////////
////////////////////////////////////////////////////////////
void makeDelaunayWithInfo(dataHolder& data)
{
    auto start = chrono::high_resolution_clock::now();

    cout << "\nDelaunay triangulation..." << endl;
    cout << "\t-" <<  data.points.size()  << " input points" << endl;

    assert(data.points.size() == data.infos.size());

    // make the triangulation
    data.Dt.insert( boost::make_zip_iterator(boost::make_tuple(data.points.begin(), data.infos.begin() )),
             boost::make_zip_iterator(boost::make_tuple(data.points.end(), data.infos.end() ) )  );
//    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple(data.points.begin(), data.infos.begin() )),
//    boost::make_zip_iterator(boost::make_tuple(data.points.end(), data.infos.end() ) )  );
//    data.Dt = Dt;


//    if(options.insert_sensor){
//        for(int i = 0; i < data.points.size(); i++)
//            data.Dt.insert(data.infos[i].sensor_positions[0]);
//    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-" << data.Dt.number_of_vertices() << " output points" << endl;
    cout << "\t-" << data.Dt.number_of_finite_facets() << " finite facets" << endl;
    cout << "\t-" << data.Dt.number_of_finite_cells() << "(+" << data.Dt.number_of_cells()-data.Dt.number_of_finite_cells() << ")" << " cells" << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;

}
// "Multi-view reconstruction preserving weakly-supported surfaces"
// and https://github.com/CGAL/cgal/tree/master/Triangulation_3/examples/Triangulation_3
// meaning, only insert point in triangulation, if no point in d-distance is already in there
// by looking for nearest neighbour in triangulation, which can be done efficiently with CGAL
// and checking if this nearest neighbour is further away then d-distance. if not increase a_vis of this neighbour by one.
void makeAdaptiveDelaunayWithInfo(dataHolder& data, double epsilon)
{
    auto start = chrono::high_resolution_clock::now();
    cout << "\nAdaptive Delaunay triangulation..." << endl;
    cout << "\t-" <<  data.points.size()  << " input points" << endl;


    if(data.scale_factor > 0.0){
        cout << "\t-scaling to unit cube is activate, so epislon will be scaled with scaling factor " << data.scale_factor << endl;
        epsilon *= data.scale_factor;
    }
    cout << "\t-with " << epsilon << " spacing" << endl;

    epsilon*=epsilon;


    // make the triangulation
//    Delaunay Dt;
    Vertex_handle v;
    v = data.Dt.insert(data.points[0]);
    v->info() = data.infos[0];
    v = data.Dt.insert(data.points[1]);
    v->info() = data.infos[1];
    v = data.Dt.insert(data.points[2]);
    v->info() = data.infos[2];
    assert( data.Dt.is_valid() );
    Cell_handle c;
    Vertex_handle current_vertex;
    Vertex_handle closest_vertex;
    Vertex_handle new_vertex;
    double dist;
    double new_dist;
    for(int i = 4; i < data.points.size(); i++){
        c=data.Dt.locate(data.points[i]);
        // TODO: define the distance to the four triangles of this cell, not to the four points!!
        dist=10e10;
        for(int ci = 0; ci < 4; ci++){
            current_vertex = c->vertex(ci);
            if(!data.Dt.is_vertex(current_vertex))
                continue;
            new_dist=CGAL::squared_distance(current_vertex->point(), data.points[i]);
            if(new_dist<dist){
                dist=new_dist;
                closest_vertex = current_vertex;
            }
        }
        if(dist>epsilon){
            new_vertex = data.Dt.insert(data.points[i]);
            new_vertex->info() = data.infos[i];
        }
        else{
            // TODO: this could be activate again, once sensor_positions is a set (controlled by the sensor index),
            // to avoid inserting the same sensor again and again for the same vertex.
            closest_vertex->info().sensor_positions.insert(closest_vertex->info().sensor_positions.end(),
                        data.infos[i].sensor_positions.begin(), data.infos[i].sensor_positions.end());

//            closest_vertex->info().alpha+=1; // increase the number of points this vertex represents
//            closest_vertex->info().number_of_rays+=1;
        }
    }

    // remove duplicate sensors

//    auto end = v.end();
//    for (auto it = v.begin(); it != end; ++it) {
//        end = std::remove(it + 1, end, *it);
//    }

//    v.erase(end, v.end());



    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);

    cout << "\t-" << data.Dt.number_of_vertices() << " output points" << endl;
    cout << "\t-" << data.Dt.number_of_finite_facets() << " finite facets" << endl;
    cout << "\t-" << data.Dt.number_of_finite_cells() << "(+" << data.Dt.number_of_cells()-data.Dt.number_of_finite_cells() << ")" << " cells" << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;
}

void indexDelaunay(dataHolder& data, runningOptions options){    
    // cells
    if(options.make_global_cell_idx){
        int gcidx = 0;
        Delaunay::All_cells_iterator aci;
        for(aci = data.Dt.all_cells_begin(); aci != data.Dt.all_cells_end(); aci++)
            aci->info().global_idx = gcidx++;
    }
    if(options.make_finite_cell_idx){
        int fcidx = 0;
        Delaunay::Finite_cells_iterator fci;
        for(fci = data.Dt.finite_cells_begin(); fci != data.Dt.finite_cells_end(); fci++)
            fci->info().finite_idx = fcidx++;
    }
    // vertex
    if(options.make_global_vertex_idx){
        int gvidx = 0;
        Delaunay::All_vertices_iterator avi;
        for(avi = data.Dt.all_vertices_begin(); avi != data.Dt.all_vertices_end(); avi++)
            avi->info().global_idx = gvidx++;
    }
    if(options.make_finite_vertex_idx){
        int fvidx = 0;
        Delaunay::Finite_vertices_iterator fvi;
        for(fvi = data.Dt.finite_vertices_begin(); fvi != data.Dt.finite_vertices_end(); fvi++)
            fvi->info().finite_idx = fvidx++;
    }
}


////////////////////////////////////////////////////////////
//////////////////////// scan object ////////////////
////////////////////////////////////////////////////////////
#include <CGAL/Random.h>
#include <random>
typedef CGAL::Side_of_triangle_mesh<SurfaceMesh, EPICK> Point_inside;
typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
bool scanObjectClosed(dataHolder& data, runningOptions& options){

    auto start = std::chrono::high_resolution_clock::now();


    cout << "\nScan closed ground truth surface: " << endl;
    cout << "\t-with " << options.number_of_points_to_scan << " points " << endl;
    cout << "\t-with Gaussian noise with mean 0.0 and std " << options.noise_std << endl;
    cout << "\t-with " << options.percentage_of_outliers * 100 << "% outliers" << endl;
    cout << "\t-and " << options.number_of_cameras << " cameras" << endl;


    // get all the triangles from the ground truth polyhedron model
    if(!data.gt_poly.is_valid()){
        cout << "GROUND TRUTH IS NOT VALID!\n" << endl;
        return 1;
    }
    if(!data.gt_poly.is_closed()){
        cout << "GROUND TRUTH IS NOT CLOSED!\n" << endl;
        return 1;
    }

    if(!data.gt_poly.is_pure_triangle())
        CGAL::Polygon_mesh_processing::triangulate_faces(data.gt_poly);


    SurfaceMesh smesh;
    CGAL::copy_face_graph(data.gt_poly, smesh);

    Tree tree(faces(smesh).first, faces(smesh).second, smesh);
//    Tree tree(faces(data.gt_poly).first, faces(data.gt_poly).second, data.gt_poly);
    tree.accelerate_distance_queries();

    // make the synthetic rays
    // this gets an axis aligned bounding box
    Bbox bb = CGAL::Polygon_mesh_processing::bbox(data.gt_poly);
    Point bb_center = Point((bb.xmin()+bb.xmax())/2,
                            (bb.ymin()+bb.ymax())/2,
                            (bb.zmin()+bb.zmax())/2);

    double radius = sqrt(pow(bb.xmin()-bb_center.x(),2)+
                     pow(bb.ymin()-bb_center.y(),2)+
                     pow(bb.zmin()-bb_center.z(),2));

    CGAL::Random random(42);

    CGAL::Random_points_on_sphere_3<Point> points_on_sphere3(3.5*radius, random); // for cameras
    CGAL::Random_points_on_sphere_3<Point> points_on_sphere2(2.5*radius, random); // for cameras
    CGAL::Random_points_on_sphere_3<Point> points_on_sphere1(radius, random); // for aims
    vector<Point> cameras, aims;

    int k = options.number_of_cameras % 2;
    CGAL::cpp11::copy_n(points_on_sphere3, options.number_of_cameras/2 + k, std::back_inserter(cameras));
    CGAL::cpp11::copy_n(points_on_sphere2, options.number_of_cameras/2, std::back_inserter(cameras));
    for(int i = 0; i < cameras.size(); i++){
        cameras[i] = Point(cameras[i].x() + bb_center.x(),
                        cameras[i].y() + bb_center.y(),
                        cameras[i].z() + bb_center.z());
        data.sensor_map.insert(pair<int,Point>(i,cameras[i]));
    }

    CGAL::cpp11::copy_n(points_on_sphere1, options.number_of_points_to_scan, std::back_inserter(aims));

    // for noise/outliers
    // TODO: add more noise to the sensors (e.g. mimic wrong ray). without adding noise to the surface
    // this could simply be done by choosing a different camera in the very last step
    default_random_engine rand_gen;
    rand_gen.seed(42);
//    uniform_real_distribution<double> dn(-radius/100,radius/100); // 1 cm noise
    normal_distribution<double> rnoise(0.0, options.noise_std); // Gaussian noise with mean 0.0 and std noise_std

    // the following two randoms are for choosing a random camera and a random aim to shoot at, or simply put, for producing a random ray
    std::uniform_int_distribution<int> rcamera(0, options.number_of_cameras-1);
    std::uniform_int_distribution<int> rpoint(0, options.number_of_points_to_scan-1);
    // scan the polyhedron by finding "first" intersections of all triangles and rays
    int n_outliers = options.number_of_points_to_scan * options.percentage_of_outliers;
    double rn1, rn2, rn3;
    while(data.points.size() < (options.number_of_points_to_scan - n_outliers)){

        // randomly choose a camera and a point from the generated cameras and points
        int sensor_idx = rcamera(rand_gen);
        Point c = cameras[sensor_idx];
        Point a = aims[rpoint(rand_gen)];
        Ray ray(c,
                Point(bb_center.x()+a.x(), bb_center.y()+a.y(), bb_center.z()+a.z()));

        Point closest_intersection_point;

        Ray_intersection intersection = tree.first_intersection(ray);
        if(intersection){
          if(boost::get<Point>(&(intersection->first))){
            const Point* p = boost::get<Point>(&(intersection->first) );
            auto fd = intersection->second;
            Vector gt_normal = PMP::compute_face_normal(fd,smesh);

            rn1 = rnoise(rand_gen);
            rn2 = rnoise(rand_gen);
            rn3 = rnoise(rand_gen);
            closest_intersection_point=Point(p->x()+rn1,
                                             p->y()+rn2,
                                             p->z()+rn3);

            data.points.push_back(closest_intersection_point);
            vertex_info vinf;
            vinf.gt_normal = gt_normal;
            vinf.normal = gt_normal;
            vinf.sensor_positions.push_back(ray.source());
            vinf.sensor_vec = ray.source() - closest_intersection_point;
            vinf.sensor_idx = sensor_idx;
            // TODO: so, e.g. here just say v_inf.sensor_idx = sensor_idx + 1 for an outlier
            vinf.outlier = 0;
            data.infos.push_back(vinf);
          }
        }
    }

//    vector<Point> outliers;
//    if(data.points.size() < options.number_of_points_to_scan){
//        CGAL::Random_points_in_sphere_3<Point> points_in_sphere(radius*3, random);
//        CGAL::cpp11::copy_n(points_in_sphere, options.number_of_points_to_scan, std::back_inserter(outliers));
//    }


    double xo,yo,zo;

    // add outliers
    while(data.points.size() < options.number_of_points_to_scan){

        // randomly choose a camera and an outlier point from the generated cameras and outlier points
        int sensor_idx = rcamera(rand_gen);
        Point c = cameras[sensor_idx];

        xo = ((double) rand() / (RAND_MAX))-0.5;
        yo = ((double) rand() / (RAND_MAX))-0.5;
        zo = ((double) rand() / (RAND_MAX))-0.5;

        Point o = Point(xo,yo,zo);
//        Point o = outliers[rpoint(rand_gen)];
        Point out = Point(bb_center.x()+o.x(), bb_center.y()+o.y(), bb_center.z()+o.z());
        Ray ray(Point(bb_center.x()+c.x(), bb_center.y()+c.y(), bb_center.z()+c.z()),
                out);

        data.points.push_back(out);
        vertex_info vinf;
        vinf.gt_normal = CGAL::NULL_VECTOR;
        vinf.normal = CGAL::NULL_VECTOR;
        vinf.sensor_positions.push_back(ray.source());
        vinf.sensor_vec = ray.source() - out;
        vinf.sensor_idx = sensor_idx;
        vinf.outlier = 1;
        data.infos.push_back(vinf);
    }
    data.has_sensor = true;
    data.has_normal = true;
    data.has_gt_normal = true;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done after " << duration.count() << "s" << endl;

    return 0;

}


bool scanObjectOpen(dataHolder& data, runningOptions& options){

    auto start = std::chrono::high_resolution_clock::now();


    cout << "\nScan closed ground truth surface: " << endl;
    cout << "\t-with " << options.number_of_points_to_scan << " points " << endl;
    cout << "\t-with Gaussian noise with mean 0.0 and std " << options.noise_std << endl;
    cout << "\t-with " << options.percentage_of_outliers * 100 << "% outliers" << endl;
    cout << "\t-and " << options.number_of_cameras << " cameras" << endl;


    // get all the triangles from the ground truth polyhedron model
    if(!data.gt_poly.is_valid()){
        cout << "GROUND TRUTH IS NOT VALID!\n" << endl;
        return 1;
    }
    if(!data.gt_poly.is_closed()){
        cout << "GROUND TRUTH IS NOT CLOSED!\n" << endl;
        return 1;
    }

    if(!data.gt_poly.is_pure_triangle())
        CGAL::Polygon_mesh_processing::triangulate_faces(data.gt_poly);


    SurfaceMesh smesh;
    CGAL::copy_face_graph(data.gt_poly, smesh);

    Tree tree(faces(smesh).first, faces(smesh).second, smesh);
    tree.accelerate_distance_queries();

    // make the synthetic rays
    // this gets an axis aligned bounding box
    Bbox bb = CGAL::Polygon_mesh_processing::bbox(data.gt_poly);
    Point bb_center = Point((bb.xmin()+bb.xmax())/2,
                            (bb.ymin()+bb.ymax())/2,
                            (bb.zmin()+bb.zmax())/2);

    double radius = sqrt(pow(bb.xmin()-bb_center.x(),2)+
                     pow(bb.ymin()-bb_center.y(),2)+
                     pow(bb.zmin()-bb_center.z(),2));

    CGAL::Random random(42);

    CGAL::Random_points_in_sphere_3<Point> points_in_sphere3(1.0*radius, random); // for cameras
    CGAL::Random_points_in_triangle_mesh_3<Polyhedron> points_in_gt(data.gt_poly,random); // for aims
    vector<Point> temp_cameras, cameras, aims;

    CGAL::cpp11::copy_n(points_in_sphere3, 10000, std::back_inserter(temp_cameras));
    // Initialize the point-in-polyhedron tester
    const Point_inside inside_tester(tree);
    int i = 0;
    while(cameras.size() < options.number_of_cameras){
        if(inside_tester(temp_cameras[i]) == CGAL::ON_UNBOUNDED_SIDE){
            Point cam = Point(temp_cameras[i].x() + bb_center.x(),
                              temp_cameras[i].y() + bb_center.y(),
                              temp_cameras[i].z() + bb_center.z());
            cameras.push_back(cam);
            data.sensor_map.insert(pair<int,Point>(i,cam));
        }
        i++;
    }

    CGAL::cpp11::copy_n(points_in_gt, options.number_of_points_to_scan, std::back_inserter(aims));

    // for noise/outliers
    // TODO: add more noise to the sensors (e.g. mimic wrong ray). without adding noise to the surface
    // this could simply be done by choosing a different camera in the very last step
    default_random_engine rand_gen;
    rand_gen.seed(42);
//    uniform_real_distribution<double> dn(-radius/100,radius/100); // 1 cm noise
    normal_distribution<double> rnoise(0.0, options.noise_std); // Gaussian noise with mean 0.0 and std noise_std

    // the following two randoms are for choosing a random camera and a random aim to shoot at, or simply put, for producing a random ray
    std::uniform_int_distribution<int> rcamera(0, options.number_of_cameras-1);
    std::uniform_int_distribution<int> rpoint(0, options.number_of_points_to_scan-1);
    // scan the polyhedron by finding "first" intersections of all triangles and rays
    int n_outliers = options.number_of_points_to_scan * options.percentage_of_outliers;
    double rn1, rn2, rn3;
    while(data.points.size() < (options.number_of_points_to_scan - n_outliers)){

        // randomly choose a camera and a point from the generated cameras and points
        int sensor_idx = rcamera(rand_gen);
        Point c = cameras[sensor_idx];
        Point a = aims[rpoint(rand_gen)];
        Ray ray(c,
                Point(bb_center.x()+a.x(), bb_center.y()+a.y(), bb_center.z()+a.z()));

        Point closest_intersection_point;

        Ray_intersection intersection = tree.first_intersection(ray);
        if(intersection){
          if(boost::get<Point>(&(intersection->first))){
            const Point* p = boost::get<Point>(&(intersection->first) );
            auto fd = intersection->second;
            Vector gt_normal = PMP::compute_face_normal(fd,smesh);

            rn1 = rnoise(rand_gen);
            rn2 = rnoise(rand_gen);
            rn3 = rnoise(rand_gen);
            closest_intersection_point=Point(p->x()+rn1,
                                             p->y()+rn2,
                                             p->z()+rn3);

            data.points.push_back(closest_intersection_point);
            vertex_info vinf;
            vinf.gt_normal = gt_normal;
            vinf.normal = gt_normal;
            vinf.sensor_positions.push_back(ray.source());
            vinf.sensor_vec = ray.source() - closest_intersection_point;
            vinf.sensor_idx = sensor_idx;
            // TODO: so, e.g. here just say v_inf.sensor_idx = sensor_idx + 1 for an outlier
            vinf.outlier = 0;
            data.infos.push_back(vinf);
          }
        }
    }

    vector<Point> outliers;
    if(data.points.size() < options.number_of_points_to_scan){
        CGAL::Random_points_in_sphere_3<Point> points_in_sphere(radius*3, random);
        CGAL::cpp11::copy_n(points_in_sphere, options.number_of_points_to_scan, std::back_inserter(outliers));
    }
    // add outliers
    while(data.points.size() < options.number_of_points_to_scan){

        // randomly choose a camera and an outlier point from the generated cameras and outlier points
        int sensor_idx = rcamera(rand_gen);
        Point c = cameras[sensor_idx];
        Point o = outliers[rpoint(rand_gen)];
        Point out = Point(bb_center.x()+o.x(), bb_center.y()+o.y(), bb_center.z()+o.z());
        Ray ray(Point(bb_center.x()+c.x(), bb_center.y()+c.y(), bb_center.z()+c.z()),
                out);

        data.points.push_back(out);
        vertex_info vinf;
        vinf.gt_normal = CGAL::NULL_VECTOR;
        vinf.normal = CGAL::NULL_VECTOR;
        vinf.sensor_positions.push_back(ray.source());
        vinf.sensor_vec = ray.source() - out;
        vinf.sensor_idx = sensor_idx;
        vinf.outlier = 1;
        data.infos.push_back(vinf);
    }
    data.has_sensor = true;
    data.has_normal = true;
    data.has_gt_normal = true;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done after " << duration.count() << "s" << endl;

    return 0;

}




//bool scanObjectOpen(dataHolder& data, runningOptions& options){

//    auto start = std::chrono::high_resolution_clock::now();


//    cout << "\nScan open ground truth surface: " << endl;
//    cout << "\t-with " << options.number_of_points_to_scan << " points " << endl;
//    cout << "\t-with Gaussian noise with mean 0.0 and std " << options.noise_std << endl;
//    cout << "\t-with " << options.percentage_of_outliers * 100 << "% outliers" << endl;
//    cout << "\t-and " << options.number_of_cameras << " cameras" << endl;

//    // get all the triangles from the ground truth polyhedron model
//    if(!data.gt_poly.is_valid()){
//        cout << "GROUND TRUTH IS NOT VALID!\n" << endl;
//        return 1;
//    }
//    if(!data.gt_poly.is_closed()){
//        cout << "GROUND TRUTH IS NOT CLOSED!\n" << endl;
//        return 1;
//    }

//    CGAL::Polygon_mesh_processing::remove_degenerate_faces(data.gt_poly);


//    if(!data.gt_poly.is_pure_triangle())
//        CGAL::Polygon_mesh_processing::triangulate_faces(data.gt_poly);


//    SurfaceMesh smesh;
//    CGAL::copy_face_graph(data.gt_poly, smesh);

//    Tree tree(faces(smesh).first, faces(smesh).second, smesh);
//    tree.accelerate_distance_queries();

//    // make the synthetic rays


//    CGAL::Random random(42);


//    // make the synthetic rays
//    // this gets an axis aligned bounding box
//    Bbox bb = CGAL::Polygon_mesh_processing::bbox(data.gt_poly);
//    Point bb_center = Point((bb.xmin()+bb.xmax())/2,
//                            (bb.ymin()+bb.ymax())/2,
//                            (bb.zmin()+bb.zmax())/2);

//    double radius = sqrt(pow(bb.xmin()-bb_center.x(),2)+
//                     pow(bb.ymin()-bb_center.y(),2)+
//                     pow(bb.zmin()-bb_center.z(),2));


//    CGAL::Random_points_on_sphere_3<Point> points_on_sphere1(1.5*radius, random); // for cameras
//    CGAL::Random_points_on_sphere_3<Point> points_on_sphere2(2.5*radius, random); // for cameras
//    CGAL::Random_points_on_sphere_3<Point> points_on_sphere3(3.5*radius, random); // for cameras
//    CGAL::Random_points_in_triangle_mesh_3<Polyhedron> points_in_gt(data.gt_poly,random); // for aims

//    vector<Point> temp_cameras, cameras, aims;
//    CGAL::cpp11::copy_n(points_in_gt, options.number_of_points_to_scan, std::back_inserter(aims));

////    if(data.sampled_sensors.empty()){
//    CGAL::cpp11::copy_n(points_on_sphere3, 200, std::back_inserter(temp_cameras));
//    CGAL::cpp11::copy_n(points_on_sphere2, 400, std::back_inserter(temp_cameras));
//    CGAL::cpp11::copy_n(points_on_sphere1, 400, std::back_inserter(temp_cameras));
//    std::uniform_int_distribution<int> tcamera(0, 999);

//    default_random_engine rand_gen;
//    rand_gen.seed(42);

//    int i = 0;
//    const Point_inside inside_tester(tree);
//    while(cameras.size()< options.number_of_cameras){
//        int rindx=tcamera(rand_gen);
//        Point c = Point(temp_cameras[rindx].x() + bb_center.x(),
//                        temp_cameras[rindx].y() + bb_center.y(),
//                        temp_cameras[rindx].z() + bb_center.z());
//        if(c.y() > bb_center.y() && inside_tester(c) == CGAL::ON_UNBOUNDED_SIDE){
//            cameras.push_back(c);
//            data.sensor_map.insert(pair<int,Point>(i++,c));
//        }
//    }

//    // for noise/outliers
//    // TODO: add more noise to the sensors (e.g. mimic wrong ray). without adding noise to the surface
//    // this could simply be done by choosing a different camera in the very last step

////    uniform_real_distribution<double> dn(-radius/100,radius/100); // 1 cm noise
//    normal_distribution<double> rnoise(0.0, options.noise_std); // Gaussian noise with mean 0.0 and std noise_std

//    // the following two randoms are for choosing a random camera and a random aim to shoot at, or simply put, for producing a random ray
//    std::uniform_int_distribution<int> rcamera(0, options.number_of_cameras-1);
//    std::uniform_int_distribution<int> rpoint(0, options.number_of_points_to_scan-1);
//    // scan the polyhedron by finding "first" intersections of all triangles and rays
//    int n_outliers = options.number_of_points_to_scan * options.percentage_of_outliers;
//    double rn1, rn2, rn3;
//    while(data.points.size() < (options.number_of_points_to_scan - n_outliers)){

//        // randomly choose a camera and a point from the generated cameras and points
//        int sensor_idx = rcamera(rand_gen);
//        Point c = cameras[sensor_idx];
//        Point a = aims[rpoint(rand_gen)];
//        Ray ray(c,a);

//        Point closest_intersection_point;

//        Ray_intersection intersection = tree.first_intersection(ray);
//        if(intersection){
//          if(boost::get<Point>(&(intersection->first))){
//            const Point* p = boost::get<Point>(&(intersection->first) );

//            auto fd = intersection->second;
//            Vector gt_normal = PMP::compute_face_normal(fd,smesh);

//            rn1 = rnoise(rand_gen);
//            rn2 = rnoise(rand_gen);
//            rn3 = rnoise(rand_gen);
//            closest_intersection_point=Point(p->x()+rn1,p->y()+rn2,p->z()+rn3);
//            data.points.push_back(closest_intersection_point);
//            vertex_info vinf;
//            vinf.gt_normal = gt_normal;
//            vinf.normal = gt_normal;
//            vinf.sensor_positions.push_back(ray.source());
//            vinf.sensor_vec = ray.source() - closest_intersection_point;
//            vinf.sensor_idx = sensor_idx;
//            // TODO: so, e.g. here just say v_inf.sensor_idx = sensor_idx + 1 for an outlier
//            vinf.outlier = 0;
//            data.infos.push_back(vinf);
//          }
//        }
//    }

//    vector<Point> outliers;
////    if(data.points.size() < options.number_of_points_to_scan){
////        CGAL::Random_points_in_sphere_3<Point> points_in_sphere(radius*3, random);
////        CGAL::cpp11::copy_n(points_in_sphere, options.number_of_points_to_scan, std::back_inserter(outliers));
////    }
////    if(data.points.size() < options.number_of_points_to_scan){
////        getOrientedBoundingBox(data.gt_poly,data.bb_array,data.bb_surface_mesh);
////        CGAL::Random_points_in_triangle_mesh_3<SurfaceMesh> points_in_bb(data.bb_surface_mesh,random); // for aims
////        CGAL::cpp11::copy_n(points_in_bb, options.number_of_points_to_scan, std::back_inserter(outliers));
////    }
//    if(data.points.size() < options.number_of_points_to_scan){
//        std::uniform_real_distribution<double> oxr(bb.xmin(),bb.xmax());
//        std::uniform_real_distribution<double> oyr(bb.ymin(),bb.ymax());
//        std::uniform_real_distribution<double> ozr(bb.zmin(),bb.zmax());
//        for(int i = 0; i < int(options.number_of_points_to_scan/7); i++){
//            outliers.push_back(Point(oxr(rand_gen)+rn1,oyr(rand_gen),ozr(rand_gen)+rn3));
//            outliers.push_back(Point(oxr(rand_gen)+rn2,oyr(rand_gen),ozr(rand_gen)));
//            outliers.push_back(Point(oxr(rand_gen)+rn3,oyr(rand_gen),ozr(rand_gen)+rn1));
//            outliers.push_back(Point(oxr(rand_gen)+rn3,oyr(rand_gen)+rn1,ozr(rand_gen)));
//            outliers.push_back(Point(oxr(rand_gen)+rn3,oyr(rand_gen)+rn2,ozr(rand_gen)));
//        }
//    }
//    std::uniform_int_distribution<int> routlier(0, outliers.size()-1);
//    // add outliers
//    double vx,vy,vz, vl;
//    while(data.points.size() < options.number_of_points_to_scan){

//        // randomly choose a camera and an outlier point from the generated cameras and outlier points
//        int sensor_idx = rcamera(rand_gen);
//        Point c = cameras[sensor_idx];
//        Point o = outliers[routlier(rand_gen)];
//        Ray ray(c,o);

//        vx = rnoise(rand_gen);
//        vy = rnoise(rand_gen);
//        vz = rnoise(rand_gen);
//        vl = sqrt(vx*vx+vy*vy+vz*vz);

//        Vector rand_normal = Vector(vx/vl,vy/vl,vz/vl);


//        data.points.push_back(o);
//        vertex_info vinf;
//        vinf.gt_normal = rand_normal;
//        vinf.normal = rand_normal;
//        vinf.sensor_positions.push_back(ray.source());
//        vinf.sensor_vec = ray.source() - o;
//        vinf.sensor_idx = sensor_idx;
//        vinf.outlier = 1;
//        data.infos.push_back(vinf);
//    }
//    data.has_sensor = true;
//    data.has_normal = true;
//    data.has_gt_normal = true;

//    auto stop = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
//    cout << "\t-done after " << duration.count() << "s" << endl;

//    return 0;

//}



////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////
#include <CGAL/Polygon_mesh_processing/repair.h>
void preprocessSensorMesh(vector<Point>& points,
                          vector<vertex_info>& infos,
                          vector<std::vector<int>>& polys){


        auto start = chrono::high_resolution_clock::now();


        assert(points.size() == infos.size());
        map<Point, vertex_info> pim;
        for(int i = 0; i < points.size(); i++){
            pim[points[i]] = infos[i];
        }

        Polyhedron sensor_mesh;
        PMP::orient_polygon_soup(points, polys);
        PMP::polygon_soup_to_polygon_mesh(points, polys, sensor_mesh);
        PMP::remove_degenerate_faces(sensor_mesh);


        vector<Point> new_points;
        vector<vertex_info> new_infos;
        Polyhedron::Vertex_iterator svi;
        int id = 0;
        for(svi = sensor_mesh.vertices_begin(); svi != sensor_mesh.vertices_end(); svi++){
            new_points.push_back(svi->point());
            vertex_info info = pim.find(svi->point())->second;
            new_infos.push_back(info);
            svi->id() = id++;
        }
        Polyhedron::Facet_iterator sfi;
        id = 0;
        for(sfi = sensor_mesh.facets_begin(); sfi != sensor_mesh.facets_end(); sfi++){
            sfi->id() = id++;
        }
//        Polyhedron::Edge_iterator sei;
//        id = 0;
//        for(sei = sensor_mesh.edges_begin(); sei != sensor_mesh.edges_end(); sei++){
//            sei->id() = id++;
//        }
        vector<std::vector<int>> new_polys;
        for(sfi = sensor_mesh.facets_begin(); sfi != sensor_mesh.facets_end(); sfi++){
            Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
            std::vector<int> ids;
            do{ids.push_back(circ->vertex()->id());}
            while ( ++circ != sfi->facet_begin());
            new_polys.push_back(ids);
        }

        points = new_points;
        infos = new_infos;
        polys = new_polys;

        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
        cout << "Mesh preprocessing done in " << duration.count() << "s" << endl;
}

////////////////////////////////////////////////////////////
/////////////////// postprocessing functions ////////////////
////////////////////////////////////////////////////////////
void edgeCollapse(dataHolder& data, meshProcessingOptions mo)
{

    auto start = chrono::high_resolution_clock::now();

    cout << "\nEdge collapse..." << endl;
    cout << "\t-keep " << mo.keep_edges*100 << "\% of edges" << endl;
    cout << "\t-before: " << data.smesh.number_of_halfedges()/2 << endl;

    // This is a stop predicate (defines when the algorithm terminates).
    SMS::Count_ratio_stop_predicate<Polyhedron> stop_criterion(mo.keep_edges);

    // This the actual call to the simplification algorithm.
    // The surface mesh and stop conditions are mandatory arguments.
    // The index maps are needed because the vertices and edges
    // of this surface mesh lack an "id()" field.
    int r = SMS::edge_collapse(data.smesh,stop_criterion);



    cout << "\t-after " << (data.smesh.number_of_halfedges()/2) << endl;
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;

}






void variationalShapeApproximation(dataHolder& data, meshProcessingOptions mo){

    auto start = chrono::high_resolution_clock::now();

    cout << "\nVariational Shape Approximation..." << endl;

    // output indexed triangles
    vector<Point> anchors;
    vector<array<size_t,3>> triangles; // triplets of indices
    // free function interface with named parameters
    bool is_manifold = VSA::approximate_triangle_mesh(data.smesh,
      CGAL::parameters::seeding_method(VSA::HIERARCHICAL). // hierarchical seeding
      max_number_of_proxies(mo.max_number_of_proxies). // seeding with maximum number of proxies
      number_of_iterations(30). // number of clustering iterations after seeding
      anchors(back_inserter(anchors)). // anchor vertices
      triangles(back_inserter(triangles))); // indexed triangles
//    cout << "\t-anchor vertices: " << anchors.size() << endl;
//    cout << "\t-triangles: " << triangles.size() << endl;

    // convert from soup to surface mesh
//    PMP::orient_polygon_soup(anchors, triangles);
//    SurfaceMesh output;

//    data.smesh.collect_garbage();
//    data.smesh.clear();
//    data.smesh = output;
    data.remaining_points.clear();
    data.remaining_facets.clear();
    data.remaining_points = anchors;
    data.remaining_facets = triangles;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;

    createSurfaceMesh(data, mo);
}

void removeLargeFaces(dataHolder& data, int factor){
    // get the area of the full mesh:
    double full_area = PMP::area(data.poly);
    int n_facets = data.poly.size_of_facets();
    cout << "\t-area of the mesh: " << full_area << endl;
    cout << "\t-number of facets of the mesh: " << n_facets << endl;
    set<Polyhedron::Facet_iterator> fs_to_remove;
    for (Polyhedron::Facet_iterator f_it = data.poly.facets_begin(); f_it != data.poly.facets_end(); ++f_it)
    {
        // bool remove_this_facet = (test if f_it should be removed)
        double face_area = PMP::face_area(f_it, data.poly);
        if(face_area > full_area * 10/n_facets)
        {
            fs_to_remove.insert(f_it);
        }
    }
    cout << "\t-number of facets that will be removed due to their size: " << fs_to_remove.size() << endl;
    for (auto fi = fs_to_remove.begin(); fi != fs_to_remove.end();  ++fi)
    {
        Polyhedron::Halfedge_around_facet_circulator f_circ = (*fi)->facet_begin();
        data.poly.erase_facet(f_circ);
    }

}


void fillHoles(SurfaceMesh& mesh){

    typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   halfedge_descriptor;
    typedef boost::graph_traits<SurfaceMesh>::face_descriptor       face_descriptor;
    typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     vertex_descriptor;

    unsigned int nb_holes = 0;
    unsigned int unfilled_holes = 0;
    std::vector<halfedge_descriptor> border_cycles;
    // collect one halfedge per boundary cycle
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
    for(halfedge_descriptor h : border_cycles)
    {

      std::vector<face_descriptor>  patch_facets;
      std::vector<vertex_descriptor> patch_vertices;
      bool success = std::get<0>(
        CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                  mesh,
                  h,
                  std::back_inserter(patch_facets),
                  std::back_inserter(patch_vertices)) );
      if(!success)
          unfilled_holes++;
      else
          nb_holes++;
//      std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
//      std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
//      std::cout << "  Is fairing successful: " << success << std::endl;
    }
    std::cout << "\t-number of filled holes: " << nb_holes << std::endl;
    if(unfilled_holes)
        std::cout << "\t-number of unfilled holes: " << nb_holes << std::endl;


}
void createSurfaceMesh(dataHolder& data, meshProcessingOptions& options)
{
    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nCreate a manifold and closed mesh..." << endl;

    data.smesh.collect_garbage();
    data.smesh.clear();
    data.smesh.collect_garbage();

    if(data.remaining_points.size() < 1){
        data.remaining_points = data.points;
        data.remaining_facets = data.facets;
    }

    cout << "\t-is manifold: " << PMP::is_polygon_soup_a_polygon_mesh(data.remaining_facets) << endl;

    if(options.try_to_make_manifold){
        cout << "\t-try to make manifold" << endl;
        bool manifold = PMP::orient_polygon_soup(data.remaining_points, data.remaining_facets);
        cout << "\t-has been made manifold: " << manifold << endl;
    }
    PMP::polygon_soup_to_polygon_mesh(data.remaining_points, data.remaining_facets, data.smesh);

    SurfaceMesh::Property_map<face_descriptor, std::size_t> fccmap =  data.smesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
    size_t components = CGAL::Polygon_mesh_processing::connected_components(data.smesh, fccmap);

    cout << "\t-number of components: " << components << endl;
    int erased_components = 0;
    if(options.number_of_components_to_keep > 0){
        erased_components = CGAL::Polygon_mesh_processing::keep_largest_connected_components(data.smesh, options.number_of_components_to_keep);
        cout << "\t-erased components: " << erased_components << endl;
    }

    bool closed = CGAL::is_closed(data.smesh);
    cout << "\t-is closed: " << closed << endl;

    if(options.try_to_close && !closed){
        cout << "\t-try to close holes" << endl;
        fillHoles(data.smesh);
        cout << "\t-has been closed: " << CGAL::is_closed(data.smesh) << endl;
    }
    // reorient
//    PMP::orient(data.smesh);
    bool oriented = PMP::is_outward_oriented(data.smesh);
    cout << "\t-is outward oriented: " << oriented << endl;



    cout << "\t-number of vertices: " << data.smesh.number_of_vertices() << endl;
    cout << "\t-number of edges: " << data.smesh.number_of_edges() << endl;
    cout << "\t-number of faces: " << data.smesh.number_of_faces() << endl;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
//    cout << "\t-Created duplicate vertices to generated polygon surface? " << oriented << endl;
//    cout << "\t-0 = additional vertices added for orientation and for ensuring manifoldness." << endl;
//    cout << "\t-Number of components of surface mesh: " << components << endl;
//    if(options.number_of_components_to_keep > 0)
//        cout << "\t-Number of erased components from surface mesh: " << erased_components << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;

    // remove the CC property map again
    data.smesh.remove_property_map(fccmap);


}


#ifdef OpenMVS
// fix with PI from here: https://github.com/cdcseacave/openMVS/issues/643
// otherwise nameclash with some CGAL headers
#pragma push_macro("PI")
#undef PI
#include "MVS.h"
#pragma pop_macro("PI")
#include <boost/filesystem.hpp>
using namespace MVS;

int omvsCleanMesh(dirHolder& dir, dataHolder& data){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nOMVS clean mesh..." << endl;
    cout << "\t-from densify file " << dir.path+"openMVS/densify_file.mvs" << endl;
    cout << "\t-from mesh " << dir.path+dir.write_file+dir.suffix+".ply" << endl;

    Scene scene(0);
    scene.Load(dir.path+"openMVS/densify_file.mvs");

    if(!scene.mesh.Load(dir.path+dir.write_file+dir.suffix+".ply")){
            cout << "ERROR in loading the mesh file" << endl;
            return 1;
    };


    int nItersFixNonManifold = 4; // also hardcoded to 4 in omvs code
    for (unsigned i=0; i<nItersFixNonManifold; ++i)
        if (!scene.mesh.FixNonManifold())
            break;

//    dir.suffix = "_nm";
//    scene.mesh.Save(dir.path+dir.write_file+dir.suffix+".ply");

    // clean the mesh
    scene.mesh.Clean(1.f, 20.f, true, 30.f, 2.f, false);
    scene.mesh.Clean(1.f, 0.f, true, 30.f, 0, false); // extra cleaning trying to close more holes
    scene.mesh.Clean(1.f, 0.f, false, 0, 0, true); // extra cleaning to remove non-manifold problems created by closing holes

    dir.suffix = "_cleaned";
    scene.mesh.Save(dir.path+dir.write_file+dir.suffix+".ply");


    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-saved to " << dir.path+dir.write_file+dir.suffix+".ply" << endl;
    cout << "\t-done after " << duration.count() << "s" << endl;

    return 0;
}

#endif


////////////////////////////////////////////////////////////
/////////////////// info functions ////////////////
////////////////////////////////////////////////////////////

Point getBoundingBoxCentroid(std::array<Point, 8>& obb_points, double z_offset){

    double mx = 0, my = 0, mz = 0;
    for(const auto p : obb_points){
        mx+=p.x();
        my+=p.y();
        mz+=p.z();
    }
    return Point(mx/8,my/8,mz/8+z_offset);

}

#include <CGAL/optimal_bounding_box.h>
void getOrientedBoundingBox(Polyhedron& poly, std::array<Point, 8>& obb_points, SurfaceMesh& obb){


    CGAL::oriented_bounding_box(poly,obb_points);
    CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                          obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb);
    CGAL::Polygon_mesh_processing::triangulate_faces(obb);

}



void calculateMeanTriangleArea(dataHolder& data){

    int n_facets = data.Dt.number_of_finite_facets();
    cout << "\nCalculate mean triangle area..." << endl;

    Delaunay::Finite_facets_iterator ffi;
    for (ffi = data.Dt.finite_facets_begin(); ffi != data.Dt.finite_facets_end(); ffi++)
        data.mean_triangle_area+=sqrt(data.Dt.triangle(*ffi).squared_area());
    data.mean_triangle_area = data.mean_triangle_area / n_facets;
    cout << "\t-mean triangle area: " << data.mean_triangle_area << endl;

}

void calculateMeanCellVolume(dataHolder& data){

    int n_cells = data.Dt.number_of_finite_cells();
    cout << "\nCalculate mean cell volume..." << endl;

    Delaunay::Finite_cells_iterator ffi;
    for(ffi = data.Dt.finite_cells_begin(); ffi != data.Dt.finite_cells_end(); ffi++){
//        Tetrahedron tet = ;
        data.mean_cell_volume+=data.Dt.tetrahedron(ffi).volume();
    }

    data.mean_cell_volume = data.mean_cell_volume / n_cells;
    cout << "\t-mean cell volume: " << data.mean_cell_volume << endl;

};


void createVertexSet(dataHolder& data){

    auto start = chrono::high_resolution_clock::now();
    cout << "\nCreate vertex set with properties..." << endl;
    // init point set
    data.point_set.clear();
    data.point_set.reserve(data.smesh.number_of_vertices());
    Point_set::Property_map<Vector> normals = data.point_set.add_normal_map().first;
    Point_set::Property_map<SurfaceMesh::Vertex_index> vindex = data.point_set.add_property_map<SurfaceMesh::Vertex_index>("v:index").first;

    SurfaceMesh::Property_map<vertex_descriptor,Vector> vnormals;
    vnormals = data.smesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_vertex_normals(data.smesh, vnormals);

    // fill point_set with the points from mesh vertices
    for(auto it = data.smesh.vertices_begin (); it < data.smesh.vertices_end (); ++it){
        // point-normal-vertex_handle map
        Point_set::iterator pit = data.point_set.insert(data.smesh.point(*it));
        normals[*pit] = vnormals[*it];
        vindex[*pit] = *it;
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;

}


