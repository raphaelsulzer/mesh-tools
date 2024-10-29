#include <base/cgal_typedefs.h>
#include <CGAL/estimate_scale.h>
#include <util/helper.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/convex_hull_3_to_face_graph.h>
#include <CGAL/grid_simplify_point_set.h>
#include <processing/evaluation.h>
#include <processing/meshProcessing.h>
#include <CGAL/optimal_bounding_box.h>

#include <IO/fileIO.h>

using namespace std;

typedef CGAL::Side_of_triangle_mesh<SurfaceMesh, EPICK> PointInsideSurfaceMesh;

////////////////////////////////////////////////////////////
//////////////////////////// INFOS /////////////////////////
////////////////////////////////////////////////////////////
int clipWithBoundingBox(vector<Point>& input_points, vector<vertex_info>& input_infos, SurfaceMesh& sm_obb){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nClip point cloud with bounding box..." << endl;
    cout << "\t-" << input_points.size() << " input points" << endl;

//    Polyhedron poly;
//    CGAL::copy_face_graph(sm_obb, poly);
    const PointInsideSurfaceMesh inside_tester(sm_obb);


//    Tree tree(faces(data.sm_obb).first, faces(data.sm_obb).second, data.sm_obb);
//    tree.accelerate_distance_queries();
//    // Initialize the point-in-polyhedron tester
//    const Point_inside inside_tester(tree);

    vector<Point> clipped_points;
    vector<vertex_info> clipped_infos;

    assert(input_points.size() == input_infos.size());

    for(int i = 0; i < input_points.size(); i++){
        if(inside_tester(input_points[i]) == CGAL::ON_BOUNDED_SIDE){
            clipped_points.push_back(input_points[i]);
            clipped_infos.push_back(input_infos[i]);
        }
    }

    input_points.clear();
    input_infos.clear();
    input_points = clipped_points;
    input_infos = clipped_infos;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-" << input_points.size() << " output points" << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;

    return 0;
}

void getOrientedBoundingBox(vector<Point>& points, SurfaceMesh& obb, runningOptions ro){

    cout << "\nMake oriented bounding box..." << endl;


    std::array<Point, 8> obb_points;
    CGAL::oriented_bounding_box(points,obb_points);
    if(ro.scale > 0.0){
        double mx = 0, my = 0, mz = 0;
        for(const auto p : obb_points){
            mx+=p.x();
            my+=p.y();
            mz+=p.z();
        }
        Point centroid(mx/8,my/8,mz/8);
        std::array<Point, 8> obb_scaled;
        int i = 0;
        for(const auto p : obb_points){
            Vector mover(ro.scale*(p.x() - centroid.x()),
                         ro.scale*(p.y() - centroid.y()),
                         ro.scale*(p.z() - centroid.z()));
            obb_scaled[i++] = p + mover;
        }
        CGAL::make_hexahedron(obb_scaled[0], obb_scaled[1], obb_scaled[2], obb_scaled[3],
                              obb_scaled[4], obb_scaled[5], obb_scaled[6], obb_scaled[7], obb);

    }
    else{
        CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                              obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb);
    }
    CGAL::Polygon_mesh_processing::triangulate_faces(obb);

}



double pcaKNN(Delaunay& Dt){

    typedef CGAL::Search_traits_3<EPICK> SearchTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nPCA to calculate noise per point... " << endl;


    // calculate the size of the smallest eigenvalue, which sould serve as a good measurement of noise. and use that as the sigma for the score computation
    // problem: this might have a good effect in noise areas, but it also weakens the few important votes in missing data areas

//    // get the point for every vertex
    Delaunay::Finite_vertices_iterator vft;
    std::vector<Point> all_points;
    for(vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        Point p = vft->point();
        all_points.push_back(p);
    }

    std::size_t NN =  CGAL::estimate_global_k_neighbor_scale(all_points);

    // build a kd-tree with all the points
    Tree tree(all_points.begin(), all_points.end());

    std::vector<double> all_eig3s;
    for(vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){

        // set up the neighborhood search function
        Neighbor_search search(tree, vft->point(), NN);

        std::vector<Vertex_handle> av;
        Dt.adjacent_vertices(vft, std::back_inserter(av));

        std::vector<Point> adj_points;
        int nn = 0;
        for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
            adj_points.push_back(it->first);
            nn++;
        }
        Point centroid = CGAL::centroid(adj_points.begin(), adj_points.end(), CGAL::Dimension_tag<0>());
        Eigen::MatrixXd m(3,nn);
        Vector p;
        for(int i = 0; i < nn; i++){
            p = adj_points[i]-centroid;
            m(0,i) = p.x();
            m(1,i) = p.y();
            m(2,i) = p.z();
        }
        if(nn<1)
            cout<<"no neighbours" << endl;
        Eigen::MatrixXd A(3,3);
        A = (m*(m.transpose()))/nn;

        // compute eigenvalues of the covariance matrix A:
        // implemented from here: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
        double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
        double eig1;
        double eig2;
        double eig3;
        Eigen::MatrixXd I(3,3);
        I = I.setIdentity();
        if (p1 == 0.0){
           eig1 = A(0,0);
           eig2 = A(1,1);
           eig3 = A(2,2);}
        else{
           double q = A.trace()/3.0;
           double p2 = pow((A(0,0) - q),2.0) + pow((A(1,1) - q),2.0) + pow((A(2,2) - q),2.0) + 2.0 * p1;
           double p = sqrt(p2 / 6.0);
           Eigen::MatrixXd B(3,3);
           B = (1.0 / p) * (A - q * I);
           double r = B.determinant() / 2.0;
           double phi;
           if (r <= -1.01)
              phi = M_PI / 3.0;
           else if (r >= 0.99)
              phi = 0.0;
           else
              phi = acos(r) / 3.0;
           // the eigenvalues satisfy eig3 <= eig2 <= eig1
           eig1 = q + 2.0 * p * cos(phi);
           eig3 = q + 2.0 * p * cos(phi + (2.0*M_PI/3.0));
           eig2 = 3.0 * q - eig1 - eig3;
        }

        // compute eigenvector of the third (smallest) eigenvalue:
//        Eigen::MatrixXd EV(3,3);
//        EV = (A-eig1*I)*(A-eig2*I);

        // save result as a map of Delaunay vertice and eigenvalue
        // save the result directly in the vertex_base of the Delaunay
        vft->info().sigma = eig3;
        all_eig3s.push_back(eig3);
    }

    size_t n = all_eig3s.size() / 2;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Estimated noise on neighbourhood size = " << NN << " neighbors." << endl;
    cout << "\t-Median noise is " << all_eig3s[n] << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;

    return all_eig3s[n];


}


double pcaKNN(Point_set& point_set){

    // this function calculates the size of the smallest eigenvalue, which sould serve as a good measurement of noise.


    typedef CGAL::Search_traits_3<EPICK> SearchTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;

    Point_set::Property_map<double> noise = point_set.add_property_map<double>("noise").first;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nPCA to calculate noise per point... " << endl;

    std::size_t NN =  CGAL::estimate_global_k_neighbor_scale(point_set);

    // build a kd-tree with all the points
    Tree tree(point_set.points().begin(), point_set.points().end());

    std::vector<double> all_eig3s;
    for(size_t idx = 0; idx < point_set.size(); idx++){

        auto p = point_set.point(idx);

        // set up the neighborhood search function
        Neighbor_search search(tree, p, NN);

        std::vector<Point> adj_points;
        int nn = 0;
        for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
            adj_points.push_back(it->first);
            nn++;
        }
        Point centroid = CGAL::centroid(adj_points.begin(), adj_points.end(), CGAL::Dimension_tag<0>());
        Eigen::MatrixXd m(3,nn);
        Vector v;
        for(int i = 0; i < nn; i++){
            v = adj_points[i]-centroid;
            m(0,i) = p.x();
            m(1,i) = p.y();
            m(2,i) = p.z();
        }
        if(nn<1)
            cout<<"no neighbours" << endl;
        Eigen::MatrixXd A(3,3);
        A = (m*(m.transpose()))/nn;

        // compute eigenvalues of the covariance matrix A:
        // implemented from here: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
        double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
        double eig1;
        double eig2;
        double eig3;
        Eigen::MatrixXd I(3,3);
        I = I.setIdentity();
        if (p1 == 0.0){
           eig1 = A(0,0);
           eig2 = A(1,1);
           eig3 = A(2,2);}
        else{
           double q = A.trace()/3.0;
           double p2 = pow((A(0,0) - q),2.0) + pow((A(1,1) - q),2.0) + pow((A(2,2) - q),2.0) + 2.0 * p1;
           double p = sqrt(p2 / 6.0);
           Eigen::MatrixXd B(3,3);
           B = (1.0 / p) * (A - q * I);
           double r = B.determinant() / 2.0;
           double phi;
           if (r <= -1.01)
              phi = M_PI / 3.0;
           else if (r >= 0.99)
              phi = 0.0;
           else
              phi = acos(r) / 3.0;
           // the eigenvalues satisfy eig3 <= eig2 <= eig1
           eig1 = q + 2.0 * p * cos(phi);
           eig3 = q + 2.0 * p * cos(phi + (2.0*M_PI/3.0));
           eig2 = 3.0 * q - eig1 - eig3;
        }

        // compute eigenvector of the third (smallest) eigenvalue:
//        Eigen::MatrixXd EV(3,3);
//        EV = (A-eig1*I)*(A-eig2*I);

        // save result as a map of Delaunay vertice and eigenvalue
        // save the result directly in the vertex_base of the Delaunay
        noise[idx]=eig3;
        all_eig3s.push_back(eig3);
    }

    size_t n = all_eig3s.size() / 2;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Estimated noise on neighbourhood size = " << NN << " neighbors." << endl;
    cout << "\t-Median noise is " << all_eig3s[n] << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;

    return all_eig3s[n];


}



// PCA with Delaunay neighborhood
void pcaDt(Delaunay& Dt){

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "PCA..." << std::endl;

    Delaunay::Finite_vertices_iterator vft;
    for(vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){

        std::vector<Vertex_handle> av;
        Dt.adjacent_vertices(vft, std::back_inserter(av));

        std::vector<Point> adj_points;
        std::vector<Vertex_handle>::iterator avt;
        int nn = 0;
        for(avt = av.begin(); avt < av.end(); avt++){
            if(!Dt.is_infinite(*avt)){
                adj_points.push_back(Dt.point(*avt));
                nn++;
            }
        }


        Point centroid = CGAL::centroid(adj_points.begin(), adj_points.end(), CGAL::Dimension_tag<0>());
        Eigen::MatrixXd m(3,nn);
        Vector v;
        for(int i = 0; i < nn; i++){
            v = adj_points[i]-centroid;
            m(0,i) = v.x();
            m(1,i) = v.y();
            m(2,i) = v.z();
        }
        if(nn<1)
            std::cout<<"no neighbours" << std::endl;
        Eigen::MatrixXd A(3,3);
        A = (m*(m.transpose()))/nn;

        // compute eigenvalues of the covariance matrix A:
        // implemented from here: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
        double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
        double eig1;
        double eig2;
        double eig3;
        double q,p2,p,r,phi;
        Eigen::MatrixXd B(3,3);
        Eigen::MatrixXd I(3,3);
        I = I.setIdentity();
        if (p1 <= 1e-8){
           eig1 = A(0,0);
           eig2 = A(1,1);
           eig3 = A(2,2);}
        else{
           q = A.trace()/3.0;
           p2 = pow((A(0,0) - q),2.0) + pow((A(1,1) - q),2.0) + pow((A(2,2) - q),2.0) + 2.0 * p1;
           p = sqrt(p2 / 6.0);
           B = (1.0 / p) * (A - q * I);
           r = B.determinant() / 2.0;
           if (r <= -1.01)
              phi = M_PI / 3.0;
           else if (r >= 0.99)
              phi = 0.0;
           else
              phi = acos(r) / 3.0;
           // the eigenvalues satisfy eig3 <= eig2 <= eig1
           eig1 = q + 2.0 * p * cos(phi);
           eig3 = q + 2.0 * p * cos(phi + (2.0*M_PI/3.0));
           eig2 = 3.0 * q - eig1 - eig3;
        }

        // compute eigenvector of the third (smallest) eigenvalue:
        Eigen::MatrixXd EV(3,3);
        EV = (A-eig1*I)*(A-eig2*I);
        vft->info().sigma = eig3;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "\t-Calculated noise per point with PCA on Delaunay neighborhood in " << duration.count() << "s" << std::endl;
}


////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////
//#ifdef Open3D
//bool pointInsideBbox(Point& p, const open3d::geometry::OrientedBoundingBox& bb){

//    const auto min = bb.GetMinBound();
//    const auto max = bb.GetMaxBound();

//    if(         p.x() > min.x() && p.x() < max.x()
//           &&   p.y() > min.y() && p.y() < max.y()
//           &&   p.z() > min.z() && p.z() < max.z())
//        return true;
//    else
//        return false;
//}
//bool pointInsideBbox(Point& p, const open3d::geometry::AxisAlignedBoundingBox& bb){

//    const auto min = bb.GetMinBound();
//    const auto max = bb.GetMaxBound();

//    if(         p.x() > min.x() && p.x() < max.x()
//           &&   p.y() > min.y() && p.y() < max.y()
//           &&   p.z() > min.z() && p.z() < max.z())
//        return true;
//    else
//        return false;
//}

////int applyCropping(dirHolder dir, dataHolder& data){

////    cout << "\nCrop point cloud with crop file..." << endl;

////    const string cropfile = dir.path + dir.crop_file + ".json";
////    cout << "\t-read cropfile " << dir.crop_file + ".json" << endl;
////    open3d::visualization::SelectionPolygonVolume bounding_poly;
////    open3d::io::ReadIJsonConvertible(cropfile, bounding_poly);

////    vector<Eigen::Vector3d> points_eigen;
////    for(int i =  0; i < data.points.size(); i++){
////        points_eigen.push_back(Eigen::Vector3d(data.points[i].x(),
////                                               data.points[i].y(),
////                                               data.points[i].z()));

////    }
////    vector<Point> cropped_points;
////    vector<vertex_info> cropped_infos;
////    auto indices = bounding_poly.CropInPolygon(points_eigen);
////    for(int i = 0; i < indices.size(); i++){
////        cropped_points.push_back(data.points[indices[i]]);
////        cropped_infos.push_back(data.infos[indices[i]]);

////    }
////    cout << "\t-decimated cloud from " << data.points.size() << " to " << cropped_points.size() << endl;

////    data.points = cropped_points;
////    data.infos = cropped_infos;

////}
//#endif


bool pointInsideBbox(Point& p, Bbox& bb){

    if(p.x() > bb.xmin() && p.x() < bb.xmax()
           && p.y() > bb.ymin() && p.y() < bb.ymax()
           && p.z() > bb.zmin() && p.z() < bb.zmax())
        return true;
    else
        return false;
}


bool pointInsideBbox(Point& p, SurfaceMesh& sm_obb){

//    Polyhedron poly;
//    CGAL::copy_face_graph(sm_obb, poly);
    const PointInsideSurfaceMesh inside_tester(sm_obb);

    return (inside_tester(p) == CGAL::ON_BOUNDED_SIDE);
}

typedef CGAL::Side_of_triangle_mesh<Polyhedron, EPICK> PointInsidePolyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> PolyhedronPrimitive;
typedef CGAL::AABB_traits_3<EPICK, PolyhedronPrimitive> PolyhedronTraits;
typedef CGAL::AABB_tree<PolyhedronTraits> PolyhedronTree;
typedef boost::optional<PolyhedronTree::Intersection_and_primitive_id<Ray>::Type> PolyhedronRayIntersection;
int applyCropping(dataHolder& data){

    cout << "\nCrop reconstruction to extend of convex hull of ground truth..." << endl;

    assert(data.gt_poly.size_of_vertices() > 0);

//    Bbox bbox = CGAL::Polygon_mesh_processing::bbox(data.gt_poly);


    // use convex hull of the polyhedron instead of bbox
    Delaunay T;
    T.insert(data.gt_poly.points_begin(), data.gt_poly.points_end());
    Polyhedron chull;
    CGAL::convex_hull_3_to_face_graph(T, chull);
    // Initialize the point-in-polyhedron tester
    // Construct AABB tree with a KdTree
    // typedef of tree is coming from learning.h
    PolyhedronTree tree(faces(chull).first, faces(chull).second, chull);
    tree.accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    const PointInsidePolyhedron inside_tester(tree);

    vector<vertex_info> infos;
    vector<Point> points;
    for(int i =  0; i < data.points.size(); i++){
//        if(pointInsideBbox(data.points[i], bbox)){
        if(inside_tester(data.points[i]) == CGAL::ON_BOUNDED_SIDE){
            points.push_back(data.points[i]);
            infos.push_back(data.infos[i]);
        }
    }

    cout << "\t-decimated cloud from " << data.points.size() << " to " << points.size() << endl;

    data.points = points;
    data.infos = infos;

}

void applyTransformationToPoint(Point& p, Eigen::Matrix4d& rt){

    Eigen::Vector4d slh(p.x(), p.y(), p.z(), 1.0);

    Eigen::Vector4d RTslh = rt*slh;
    p = Point(RTslh[0] / RTslh[3],
                RTslh[1] / RTslh[3],
                RTslh[2] / RTslh[3]);
}

int applyTransformationMatrix(dataHolder& data){

    cout << "\nApply allignment matrix to " << data.points.size()  << " points" << std::endl;

    assert(data.points.size()==data.infos.size());

    for(int i = 0; i < data.points.size(); i++){

        applyTransformationToPoint(data.points[i], data.transformation_matrix);

        assert(data.infos[i].sensor_positions.size() > 0);
        for(int j = 0; j < data.infos[i].sensor_positions.size(); j++)
            applyTransformationToPoint(data.infos[i].sensor_positions[j], data.transformation_matrix);

        data.infos[i].sensor_vec = Vector(data.infos[i].sensor_positions[0] - data.points[i]);
    }
    return 0;
}

void standardizePointSet(dirHolder dir, dataHolder& data, double scale_factor){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nScale " << data.points.size()  << " points..." << std::endl;

    Point centroid = CGAL::centroid(data.points.begin(), data.points.end());
    Iso_cuboid bb = CGAL::bounding_box(data.points.begin(), data.points.end());
    double scale = scale_factor/max(max(bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin()), bb.zmax()-bb.zmin());

    for(int i=0; i < data.points.size(); i++){

        // standardize the point coordinates
        data.points[i] = Point((data.points[i].x() - centroid.x()) * scale,
                               (data.points[i].y() - centroid.y()) * scale,
                               (data.points[i].z() - centroid.z()) * scale);
        // standardize the sensor_positions
        for(int j = 0; j < data.infos[i].sensor_positions.size(); j++){

            data.infos[i].sensor_positions[j] = Point((data.infos[i].sensor_positions[j].x() - centroid.x()) * scale,
                                                      (data.infos[i].sensor_positions[j].y() - centroid.y()) * scale,
                                                      (data.infos[i].sensor_positions[j].z() - centroid.z()) * scale);
        }
        data.infos[i].sensor_vec = Vector(data.infos[i].sensor_positions[0] - data.points[i]);
    }
    map<int, Point> standardized_sensor_map;
    for(int i=0; i < data.sensor_map.size(); i++){
        auto current = data.sensor_map.find(i);
        standardized_sensor_map[current->first]=Point((current->second.x() - centroid.x()) * scale,
                                           (current->second.y() - centroid.y()) * scale,
                                           (current->second.z() - centroid.z()) * scale);
    }
    data.sensor_map = standardized_sensor_map;

//    data.translation_vector = centroid;
//    data.scale_factor = scale;

//    // export to file
//    ofstream scale_file;
//    scale_file.open(dir.path+dir.read_file+"_std.txt");
//    scale_file << "# translation vector // scale factor" << endl;
//    scale_file << data.translation_vector << endl;
//    scale_file << data.scale_factor << endl;
//    scale_file.close();

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-scaled points to unit cube with " << endl;
    cout << "\t-meaning shifted centroid to Point(0,0,0) and scaled points so max. edge length of bounding box is 1" << endl;
    cout << "\t-Parameters applied:" << endl;
    cout << "\t\t-Translation vector: " << centroid << endl;
    cout << "\t\t-Scale: " << scale << endl;
    cout << "\t-exported data to " << dir.write_file+"_std.txt" << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;
}

//#include <CGAL/optimal_bounding_box.h>
//void scalePointSet(dirHolder dir, dataHolder& data, double scale){

//    auto start = std::chrono::high_resolution_clock::now();
//    cout << "\nScale " << data.points.size()  << " points" << std::endl;

//    std::array<Point, 8> obb_points;
//    CGAL::oriented_bounding_box(data.points,obb_points);
//    double mx = 0, my = 0, mz = 0;
//    for(const auto p : obb_points){
//        mx+=p.x();
//        my+=p.y();
//        mz+=p.z();
//    }
//    Point centroid(mx/8,my/8,mz/8);


//    makeDelaunayWithInfo(data);
//    double thisScale = calcMeanEdgeLength(data);
//    data.Dt.clear();

//    scale = scale / thisScale;


//    //    std::array<Point, 8> obb_scaled;
//    //    int i = 0;
//    //    for(const auto p : obb_points){
//    //        Vector mover(0.1*(p.x() - centroid.x()),
//    //                     0.1*(p.y() - centroid.y()),
//    //                     0.1*(p.z() - centroid.z()));
//    //        obb_scaled[i++] = p + mover;
//    //    }
//    //    SurfaceMesh obb_sm;
//    //    CGAL::make_hexahedron(obb_scaled[0], obb_scaled[1], obb_scaled[2], obb_scaled[3],
//    //                          obb_scaled[4], obb_scaled[5], obb_scaled[6], obb_scaled[7], obb_sm);

//    for(int i=0; i < data.points.size(); i++){

//        // standardize the point coordinates
//        data.points[i] = Point((data.points[i].x() - centroid.x()) * scale,
//                               (data.points[i].y() - centroid.y()) * scale,
//                               (data.points[i].z() - centroid.z()) * scale);
//        // standardize the sensor_positions
//        for(int j = 0; j < data.infos[i].sensor_positions.size(); j++){

//            data.infos[i].sensor_positions[j] = Point((data.infos[i].sensor_positions[j].x() - centroid.x()) * scale,
//                                                      (data.infos[i].sensor_positions[j].y() - centroid.y()) * scale,
//                                                      (data.infos[i].sensor_positions[j].z() - centroid.z()) * scale);
//        }
//        data.infos[i].sensor_vec = Vector(data.infos[i].sensor_positions[0] - data.points[i]);
//    }
//    map<int, Point> standardized_sensor_map;
//    for(int i=0; i < data.sensor_map.size(); i++){
//        auto current = data.sensor_map.find(i);
//        standardized_sensor_map[current->first]=Point((current->second.x() - centroid.x()) * scale,
//                                           (current->second.y() - centroid.y()) * scale,
//                                           (current->second.z() - centroid.z()) * scale);
//    }
//    data.sensor_map = standardized_sensor_map;

//    data.translation_vector = centroid;
//    data.scale_factor = scale;

////    // export to file
////    ofstream scale_file;
////    scale_file.open(dir.path+dir.read_file+"_std.txt");
////    scale_file << "# translation vector // scale factor" << endl;
////    scale_file << data.translation_vector << endl;
////    scale_file << data.scale_factor << endl;
////    scale_file.close();

//    auto stop = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
////    cout << "\t-scaled points to unit cube with " << endl;
////    cout << "\t-meaning shifted centroid to Point(0,0,0) and scaled points so max. edge length of bounding box is 1" << endl;
////    cout << "\t-Parameters applied:" << endl;
////    cout << "\t\t-Translation vector: " << centroid << endl;
////    cout << "\t\t-Scale: " << scale << endl;
////    cout << "\t-exported data to " << dir.write_file+"_std.txt" << endl;
////    cout << "\t-done in " << duration.count() << "s" << endl;
//}




int reStandardizePointSet(dirHolder dir, dataHolder& data){

    // this function is most likely unnecessary, because I never load geometry from a file, only topology

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nRe-Scale " << data.points.size()  << " points" << std::endl;

    ifstream file(dir.path+dir.read_file+"_std.txt");

    cout << "\t-Load standardization data from " << dir.read_file+"_std.txt" << endl;

    if(!file){
        cout << "\nFILE DOES NOT EXIST OR IS EMPTY!" << endl;
        return 1;
    }

    string line;
    getline(file, line); // skipping the header
    getline(file, line); // reading the point
    vector<string> lineAsVector;
    splitString(line, lineAsVector, ' ');
    Point centroid = Point(atof(lineAsVector[0].c_str()),
            atof(lineAsVector[1].c_str()),
            atof(lineAsVector[2].c_str()));
    getline(file, line); // reading the scale
    double scale = atof(line.c_str());

    for(int i=0; i < data.points.size(); i++){

        // standardize the point coordinates
        data.points[i] = Point((data.points[i].x() / scale) + centroid.x(),
                               (data.points[i].y() / scale) + centroid.y(),
                               (data.points[i].z() / scale) + centroid.z());
        // standardize the sensor_positions
        for(int j = 0; j < data.infos[i].sensor_positions.size(); j++){

            data.infos[i].sensor_positions[j] = Point((data.infos[i].sensor_positions[j].x() / scale) + centroid.x(),
                                                      (data.infos[i].sensor_positions[j].y() / scale) + centroid.y(),
                                                      (data.infos[i].sensor_positions[j].z() / scale) + centroid.z());
        }
        data.infos[i].sensor_vec = Vector(data.infos[i].sensor_positions[0] - data.points[i]);
    }
    map<int, Point> standardized_sensor_map;
    for(int i=0; i < data.sensor_map.size(); i++){
        auto current = data.sensor_map.find(i);
        standardized_sensor_map[current->first]=Point(  (current->second.x() / scale) + centroid.x(),
                                                        (current->second.y() / scale) + centroid.y(),
                                                        (current->second.z() / scale) + centroid.z());
    }
    data.sensor_map = standardized_sensor_map;

    // do not set these values, here, to prevent them to be applied again in the adaptive Delaunay triangulation
//    data.translation_vector = centroid;
//    data.scale_factor = scale;


    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);

    cout << "\t-Parameters read and re-applied:" << endl;
    cout << "\t\t-Translation vector: " << centroid << endl;
    cout << "\t\t-Scale: " << scale << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;

    return 0;
}


void gridSimplify(vector<Point>& points, double cell_size){

    int before = points.size();

    cout << "\nGrid simplify "<< before << " points with cell size " << cell_size << endl;

    points.erase(CGAL::grid_simplify_point_set(points, cell_size),
                 points.end());
    // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
    std::vector<Point>(points).swap(points);


    cout << "\t-from " << before << " to " << points.size() << std::endl;
}


void gridSimplify(dataHolder& data, double cell_size){

    int before = data.points.size();

    cout << "\nGrid simplify "<< before << " points with cell size " << cell_size << endl;

    std::vector<std::size_t> indices(data.points.size());
    for(std::size_t i = 0; i < data.points.size(); ++i){
      indices[i] = i;
    }
    // simplification by clustering using erase-remove idiom
    std::vector<std::size_t>::iterator end;
    end = CGAL::grid_simplify_point_set(indices,
                                        cell_size,
                                        CGAL::parameters::point_map (CGAL::make_property_map(data.points)));
    std::size_t k = end - indices.begin();
    std::vector<Point> tmp_points(k);
    std::vector<vertex_info> tmp_infos(k);
    for(std::size_t i=0; i<k; ++i){
        tmp_points[i] = data.points[indices[i]];
        tmp_infos[i] = data.infos[indices[i]];
    }
    data.points.clear();
    data.infos.clear();
    for(std::size_t i=0; i<k; ++i){
        data.points.push_back(tmp_points[i]);
        data.infos.push_back(tmp_infos[i]);
    }
    cout << "\t-from " << before << " to " << data.points.size() << std::endl;
}


void sampleFun(vector<Point>& points_in, vector<vertex_info>& infos_in,
               const int percentage){

    int before = points_in.size();
    std::vector<Point> points_out;
    std::vector<vertex_info> infos_out;

    int mod = static_cast<int>(100/percentage);
    for(int i = 0; i < points_in.size(); i++){

        if(!(i%mod)){
            points_out.push_back(points_in[i]);
            infos_out.push_back(infos_in[i]);
        }
    }
    points_in.clear();
    infos_in.clear();
    points_in = points_out;
    infos_in = infos_out;
    points_in.shrink_to_fit();
    infos_in.shrink_to_fit();

    std::cout << "Points subsampled from " << before << " to " << points_in.size() << std::endl;
}
void sampleFun(std::vector<Point>& points_in,
               const int percentage){

    int before = points_in.size();
    std::vector<Point> points_out;

    int mod = static_cast<int>(100/percentage);
    for(int i = 0; i < points_in.size(); i++){
        if(!(i%mod)){
            points_out.push_back(points_in[i]);
        }
    }
    points_in.clear();
    points_in = points_out;
    points_in.shrink_to_fit();

    std::cout << "Points subsampled from " << before << " to " << points_in.size() << std::endl;
}




void createPointSet(dataHolder& data){

    cout << "\nCreate point set..." << endl;
    cout << "\t-with ";

    // This function converts a vector<Point> and vector<vertex_info> to a Point_set_3 with corresponding property maps

    Point_set::Property_map<unsigned char> reds, greens, blues;
    Point_set::Property_map<Vector> normals;
    Point_set::Property_map<Vector> sensor_vecs;
    Point_set::Property_map<vector<Point>> sensor_positions;

    if(data.has_color){
        cout << "colors, ";
        // unfortunately, for allowing to directly export colors with CGAL::write_ply_point_set(), one has to make them as simple types
        // and not e.g. as an object Color
        reds = data.point_set.add_property_map<unsigned char>("red").first;
        greens = data.point_set.add_property_map<unsigned char>("green").first;
        blues = data.point_set.add_property_map<unsigned char>("blue").first;
    }
    if(data.has_normal)
        cout << "normals, ";
        normals = data.point_set.add_normal_map().first;
    if(data.has_sensor){
        cout << "sensor positions";
        // for allowing to directly export sensor_vector with CGAL::write_ply_point_set(),
        // one has to call data.point_set.remove_normal_map() for doing that
        // what would be more practical is simply put the sensor vectors on the normal property map before exporting
//        data.point_set.add_property_map<double>("nx", 0.0).first;
//        data.point_set.add_property_map<double>("ny", 0.0).first;
//        data.point_set.add_property_map<double>("nz", 0.0).first;

        sensor_vecs = data.point_set.add_property_map<Vector>("sensor_vec", CGAL::NULL_VECTOR).first;
        sensor_positions = data.point_set.add_property_map<vector<Point>>("sensor_positions").first;
    }

    cout << endl;

    for(size_t idx = 0; idx < data.points.size(); idx++){

//        const Point_set::Index i = idx;
        Point_set::iterator it = data.point_set.insert(data.points[idx]);

        if(data.has_color){
            reds[*it] = (unsigned char)(data.infos[idx].color.r());
            greens[*it] = (unsigned char)(data.infos[idx].color.g());
            blues[*it] = (unsigned char)(data.infos[idx].color.b());
        }
        if(data.has_normal){
            normals[*it] = data.infos[idx].normal;
        }
        if(data.has_sensor){
            sensor_vecs[*it] = data.infos[idx].sensor_vec;
            for(const auto v : data.infos[idx].sensor_positions){
                sensor_positions[*it].push_back(v);
            }
        }
    }
}







