#include <base/cgal_typedefs.h>
#include <CGAL/estimate_scale.h>
#include <util/helper.h>
#include <processing/normalAndSensorProcessing.h>

#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>

#include <CGAL/mst_orient_normals.h>

#include <IO/fileIO.h>

typedef CGAL::Search_traits_3<EPICK> SearchTraits;

typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

using namespace std;

////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////
//// estimate normals of a point set
///
///
//TODO: rewrite this with a Point_set_3, (see structuring.cpp on how to do it), so I can get rid of the PointVectorPair in cgal_typedefs.h
//int estimateNormals(dataHolder& data, int nb_neighbors, int orient)
//{
//    // WARNING: This function may fail when Boost version 1.54 is used, because of the following bug: https://svn.boost.org/trac/boost/ticket/9012

//    auto start = std::chrono::high_resolution_clock::now();
//    cout << "\nEstimate normals of " << data.point_set.size()  << " points" << endl;





//    // add the points as the first element of the point vector pair
//    for(std::size_t i=0; i < data.points.size(); ++i){

//        auto it = data.point_set.insert(data.points[i]);
//        if(data.has_normal){
//            auto normals = data.point_set.add_normal_map().first;
//            normals[*it] = data.infos[i].normal;
//        }
//    }
//    if(data.has_normal){
//            cout << "\t-existing normals copied to point_set. Set data.has_normal to false before" << endl;
//            cout << "\t-calling this function if you want to recalculate the normals." << endl;
//            return 0;
//    }

//    if(nb_neighbors == 0){
//        nb_neighbors = CGAL::estimate_global_k_neighbor_scale(data.point_set);
//        cout << "\t-estimated neighborhood size to " << nb_neighbors << endl;
//    }
//    cout << "\t-using neighborhood size " << nb_neighbors << endl;


//    // from: https://doc.cgal.org/latest/Point_set_processing_3/index.html#Point_set_processing_3Scale
//    CGAL::pca_estimate_normals<Concurrency_tag>
//      (data.point_set, nb_neighbors,
//       data.point_set.parameters());

//    // Orients normals either towards sensor if its there, or with function.
//    if(!orient){
//        cout << "\t-normals will not be oriented" << endl;
//    }
//    else if(orient==1){
//        if(data.has_sensor){
//            cout << "\t-orient normals towards sensor" << endl;
//            for(int i = 0; i < pointVectorPairs.size(); i++){
//                // if dot product is smaller than zero, they are facing in opposite directions
//                if(pointVectorPairs[i].second * data.infos[i].sensor_vec <= 0)
//                    pointVectorPairs[i].second = -pointVectorPairs[i].second;
//            }
//        }
//        else{
//            cout << "ERROR: orient was set to sensor orientation, but data has no sensor." << endl;
//            return 1;
//        }
//    }
//    else if(orient==2){
//    // Note: mst_orient_normals() requires a range of points
//    // as well as property maps to access each point's position and normal.
//        cout << "\t-orient normals with mst" << endl;

//        std::vector<PointVectorPair>::iterator unoriented_points_it =    // this returns unoriented points that could be removed in the next step
//                CGAL::mst_orient_normals(pointVectorPairs,
//                                         nb_neighbors,
//                                         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
//                                                     normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
//    }
//    else{
//        cout << "ERROR: not a valid normal orientation method, set either to 0,1 or 2." << endl;
//        return 1;
//    }

//    data.pVP = pointVectorPairs;
//    assert(data.infos.size() == data.points.size());
//    for(std::size_t i=0; i < data.points.size(); i++){
//        data.infos[i].normal = pointVectorPairs[i].second;
//    }
//    data.has_normal = true;

//    auto stop = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
////    cout << "\t-" << unoriented_points_it->end() - unoriented_points_it->begin()  << " points could not be oriented" << endl;
//    cout << "\t-done in " << duration.count() << "s" << endl;
//}







int estimateNormals(dataHolder& data, runningOptions ro)
{
    // WARNING: This function may fail when Boost version 1.54 is used, because of the following bug: https://svn.boost.org/trac/boost/ticket/9012

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nEstimating normals..." << endl;
    cout << "\t-of " << data.points.size()  << " points" << endl;
    cout << "\t-using " << ro.normal_method << endl;
    if(ro.normal_orientation == 1)
        cout << "\t-orient normals using sensor positions" << endl;
    else if(ro.normal_orientation == 2)
        cout << "\t-orient normals using mst" << endl;
    else
        cout << "\t-without orientation" << endl;

    // initialise a vector of point-vector-pairs
    std::vector<PointVectorPair> pointVectorPairs(data.points.size());


    // add the points as the first element of the point vector pair
    for(std::size_t i=0; i < data.points.size(); ++i){
        pointVectorPairs[i].first = data.points[i];
        if(data.has_normal){
            pointVectorPairs[i].second = data.infos[i].normal;
        }
    }
    if(data.has_normal){
            cout << "\t-existing normals copied to data.pVP. Set data.has_normal to false before" << endl;
            cout << "\t-calling this function if you want to recalculate the normals." << endl;
            data.pVP = pointVectorPairs;
            return 0;
    }

    int nb = ro.normal_neighborhood;
    if(ro.normal_neighborhood == 0){
        nb = CGAL::estimate_global_k_neighbor_scale(data.points);
        cout << "\t-estimated neighborhood size to " << nb << endl;
    }
    cout << "\t-using neighborhood size " << nb << endl;


    // from: https://doc.cgal.org/latest/Point_set_processing_3/index.html#Point_set_processing_3Scale
    if(ro.normal_method == "pca"){
        CGAL::pca_estimate_normals<Concurrency_tag>
          (pointVectorPairs, nb,
           CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
           normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
    }
    else if(ro.normal_method == "jet"){
        CGAL::jet_estimate_normals<Concurrency_tag>
          (pointVectorPairs, nb,
           CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
           normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
    }
//    else if(ro.normal_method == "vcm"){
//        CGAL::vcm_estimate_normals<Concurrency_tag>
//          (pointVectorPairs, nb,
//           CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
//           normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
//    }
    else{
        cout << "\ERROR: not a valid method for estimating normals. Available methods are: pca, jet, vcm." << endl;
        return 1;
    }

    int flipped = 0;

    // Orients normals either towards sensor if its there, or with function.
    if(!ro.normal_orientation){
        cout << "\t-normals will not be oriented" << endl;
    }
    else if(ro.normal_orientation==1){
        assert(pointVectorPairs.size() == data.infos.size());
        if(data.has_sensor){
            cout << "\t-orient normals towards sensor" << endl;

            for(int i = 0; i < pointVectorPairs.size(); i++){
//                cout << "before: " << pointVectorPairs[i].second << endl;
//                cout << "sensor vec: " << data.infos[i].sensor_vec << endl;
//                cout << "scalar product: " << CGAL::scalar_product(pointVectorPairs[i].second, data.infos[i].sensor_vec) << endl;
//                double s = data.infos[i].sensor_vec.x() * pointVectorPairs[i].second.x()
//                        + data.infos[i].sensor_vec.y() * pointVectorPairs[i].second.y()
//                        + data.infos[i].sensor_vec.z() * pointVectorPairs[i].second.z();
//                cout << "s: " << s << endl;
                // if dot product is smaller than zero, they are facing in opposite directions
                if(CGAL::scalar_product(pointVectorPairs[i].second, data.infos[i].sensor_vec) < 0.0){
                    flipped++;
                    pointVectorPairs[i].second = -pointVectorPairs[i].second;
                }
//                cout << "after: " << pointVectorPairs[i].second << endl;
//                cout << endl;

            }
            cout << "\t\t-flipped " << flipped << "/" << pointVectorPairs.size() << endl;

        }
        else{
            cout << "ERROR: orient was set to sensor orientation, but data has no sensor." << endl;
            return 1;
        }
    }
    else if(ro.normal_orientation==2){
    // Note: mst_orient_normals() requires a range of points
    // as well as property maps to access each point's position and normal.
        cout << "\t-orient normals with mst" << endl;

        std::vector<PointVectorPair>::iterator unoriented_points_it =    // this returns unoriented points that could be removed in the next step
                CGAL::mst_orient_normals(pointVectorPairs,
                                         nb,
                                         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                                                     normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
    }
    else{
        cout << "ERROR: not a valid normal orientation method, set either to 0,1 or 2." << endl;
        return 1;
    }


    data.pVP = pointVectorPairs;
    assert(data.infos.size() == data.points.size());
    for(std::size_t i=0; i < data.points.size(); i++){
        data.infos[i].normal = pointVectorPairs[i].second;
    }
    data.has_normal = true;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
//    cout << "\t-" << unoriented_points_it->end() - unoriented_points_it->begin()  << " points could not be oriented" << endl;
    cout << "\t-done in " << duration.count() << "s" << endl;

    return 0;
}




#include <CGAL/Orthogonal_incremental_neighbor_search.h>
typedef boost::tuple<Point, vertex_info>                            point_info;
typedef CGAL::Search_traits_3<EPICK> SearchTraits;
typedef CGAL::Search_traits_adapter<point_info,
  CGAL::Nth_of_tuple_property_map<0, point_info>,
  SearchTraits>                                               SearchTraitsAdapted;

typedef CGAL::Orthogonal_incremental_neighbor_search<SearchTraitsAdapted> Incremental_neighbor_search;
typedef Incremental_neighbor_search::Tree Incremental_Tree;
int orientNormalsWithSecondCloud(dataHolder& data1, dataHolder& data2){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nOrient normals with nearest neighbor in second cloud..." << endl;

    assert(data1.points.size()==data1.infos.size());
    assert(data2.points.size()==data2.infos.size());

    Incremental_Tree tree(boost::make_zip_iterator(boost::make_tuple( data1.points.begin(),data1.infos.begin() )),
              boost::make_zip_iterator(boost::make_tuple( data1.points.end(), data1.infos.end() ) ));

    for(int i = 0; i < data2.points.size(); i++)
    {
        Incremental_neighbor_search search(tree, data2.points[i]);
        Incremental_neighbor_search::iterator it = search.begin();

        // index of nearest neighbor in AP
        int idx = boost::get<1>(it->first).global_idx;

        // if dot product is smaller than zero, they are facing in opposite directions
        if(data2.infos[i].normal * data1.infos[idx].sensor_vec <= 0)
            data2.infos[i].normal = -data2.infos[i].normal;
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;
}



bool sortSensors(const std::pair<double, Point> &a, const std::pair<double, Point> &b)
{return (a.first < b.first);}

void orderSensors(dataHolder& data){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nOrder sensors according to surface normal..." << endl;

    if(!data.has_normal){
        runningOptions ro;
        ro.normal_method = "pca";
        ro.normal_neighborhood = 0;
        ro.normal_orientation = 1;
        estimateNormals(data,ro);
    }



    for(int i=0; i < data.pVP.size(); i++){

        std::vector<std::pair<double,Point>> sensor_points;
        for(int j = 0; j<data.infos[i].sensor_positions.size(); j++){

            Vector v1 = data.infos[i].sensor_positions[j] - data.pVP[i].first;
            Vector v2 = data.infos[i].normal;
            double angle = std::acos( v1 * v2 / std::sqrt(v1.squared_length()*v2.squared_length()));
            if(CGAL::cross_product(v1,v2).z()<0)
                angle = 2*M_PI-angle;
            sensor_points.push_back(std::make_pair(angle, data.infos[i].sensor_positions[j]));
        }
        sort(sensor_points.begin(), sensor_points.end(), sortSensors);
        data.infos[i].sensor_vec = sensor_points[0].second-data.pVP[i].first;
        // clear the unordered positions
        data.infos[i].sensor_positions.clear();
        for(auto& sensor_position : sensor_points)
            data.infos[i].sensor_positions.push_back(sensor_position.second);
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;
}



