#include <exe/inputParser.h>
#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <IO/simplificationIO.h>
#include <util/helper.h>
#include <util/vectorArithmetic.h>
#include <processing/meshProcessing.h>
#include <processing/evaluation.h>
#include <processing/pointSetProcessing.h>
#include <exe/simplification.h>

#include "boost/tuple/tuple.hpp"
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/structure_point_set.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Optimal_bounding_box/oriented_bounding_box.h>
#include <CGAL/estimate_scale.h>

bool sort_by_area(pair<Polygon_2, double>& a,
                    pair<Polygon_2, double>& b)
{
    return a.second > b.second;
}

bool sort_by_dist(pair<int, double>& a,
                    pair<int, double>& b)
{
    return a.second < b.second;
}


void exportSegmentation(dirHolder& dir, dataHolder& data){


    ofstream out(dir.path+dir.write_file+dir.suffix);

    Point_set::Property_map<SurfaceMesh::Vertex_index> vpoint = data.point_set.property_map<SurfaceMesh::Vertex_index>("v:index").first;
    SurfaceMesh::Property_map<vertex_descriptor,int> vindex = data.smesh.add_property_map<vertex_descriptor, int>("v:index").first;

    // vertices
    out << "!vertices: ";
    int i = 0;
    for(auto v : data.smesh.vertices()){
        out << i << " ";
        vindex[v] = i;
        i++;
    }
    out << endl;

    // planes
    RansacPlaneRange planes = data.efficient_ransac.planes();
    auto it = planes.begin();
    while(it != planes.end()){
        EPICK::Plane_3 current_plane(*it->get());
        out << "#plane: " << current_plane.a() << " " << current_plane.b() << " " << current_plane.c() << " " << current_plane.d() << " ";
        out << "#area: " << (*it)->indices_of_assigned_points().size() << " #vertices: ";
        vector<size_t>::const_iterator index_it = (*it)->indices_of_assigned_points().begin();
        while (index_it != (*it)->indices_of_assigned_points().end()){
//            // Retrieve point.
            const auto pt = *(data.point_set.begin() + (*index_it));
            const auto vhandle = vpoint[pt];



//            auto handle = delaunay.insert(current_plane.to_2d(data.point_set.point(pt)));
//            // keep the 3D point for exporting
//            handle->info() = current_plane.projection(data.point_set.point(pt));
//            // Proceed with the next point.
//            index_it++;
            out << vindex[vhandle] << " ";
            index_it++;
        }
        out << endl;
        // Proceed with the next detected shape.
        it++;
    }

    data.smesh.remove_property_map(vindex);

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// RANSAC //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ransac(dirHolder& dir, dataHolder& data,  meshProcessingOptions& options){

    cout << "\nPlane detection..." << endl;
    auto start = std::chrono::high_resolution_clock::now();

    // RANSAC
//    Efficient_ransac efficient_ransac;
    Efficient_ransac::Parameters parameters;
//    efficient_ransac.set_input(data.ransac_points); // Call built-in property map);
    data.efficient_ransac.set_input(data.point_set,
                               data.point_set.point_map(),
                               data.point_set.normal_map()); // Call built-in property map);
    data.efficient_ransac.add_shape_factory<RansacPlane>();

    // get optimal bounding box
    parameters.epsilon = options.ransac_epsilon;
    if(!parameters.epsilon > 0.0){
        CGAL::oriented_bounding_box(data.points, data.bounding_box);
        parameters.epsilon = 0.01*sqrt(max(CGAL::squared_distance(data.bounding_box[0],data.bounding_box[2]),CGAL::squared_distance(data.bounding_box[0],data.bounding_box[4])));
        options.ransac_epsilon = parameters.epsilon;
    }
    parameters.normal_threshold = options.ransac_normal;
    cout << "\t-distance threshold: " << parameters.epsilon << endl;
    cout << "\t-normal threshold: " << parameters.normal_threshold << endl;
    // run
    data.efficient_ransac.detect(parameters);
    // get planes
    auto planes = data.efficient_ransac.planes();
    auto index_map = CGAL::Shape_detection::Point_to_shape_index_map<RansacTraits>(data.point_set, planes);

    for(const auto p : planes){
        data.ransac_planes.push_back(*p);
    }

    cout << "\t-number of detected planes: " << planes.size() << endl;

    if(options.ransac_regularization){
        cout << "\t-regularize planes" << endl;
        CGAL::regularize_planes(data.point_set,
//                                Point_map(),
                                data.point_set.point_map(),
                                planes,
                                CGAL::Shape_detection::Plane_map<RansacTraits>(),
                                index_map,
                                true,  // regularize parallelism
                                true,  // regularize orthogonality
                                false, // do not regularize coplanarity
                                true,  // regularize Z-symmetry (default)
                                10);   // 10 degrees of tolerance for parallelism / orthogonality
        cout << "\t-number of planes after regularization: " << planes.size() << endl;
    }

    const int nc=8;
    array<CGAL::Color,nc> color_palette{CGAL::red(), CGAL::blue(), CGAL::green(), CGAL::orange(), CGAL::gray(), CGAL::purple(), CGAL::violet(), CGAL::yellow()};

    Point_set::Property_map<unsigned char> red, green, blue;
    unsigned char white = 255;
    red = data.point_set.add_property_map<unsigned char>("red", white).first;
    green = data.point_set.add_property_map<unsigned char>("green", white).first;
    blue = data.point_set.add_property_map<unsigned char>("blue", white).first;

    int i = 0;
    for(Point_set::const_iterator pit = data.point_set.begin(); pit !=data.point_set.end(); pit++){
        int c = get(index_map, i++);
        if(c == -1)
            continue; // leave unassigned points white
        red[*pit]=(unsigned char)color_palette[c%nc].r();
        green[*pit]=(unsigned char)color_palette[c%nc].g();
        blue[*pit]=(unsigned char)color_palette[c%nc].b();
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done after "<< duration.count() << "s" << endl;
}


int pointInAlphaShape(Polygon_2 as_points, EPICK::Point_2 pt){

    // taken from here: https://doc.cgal.org/latest/Polygon/Polygon_2polygon_algorithms_8cpp-example.html
    switch(CGAL::bounded_side_2(as_points.begin (), as_points.end(), pt, EPICK())) {
      case CGAL::ON_BOUNDED_SIDE:
//        cout << " is inside the polygon.\n";
        return 1;
      case CGAL::ON_BOUNDARY:
//        cout << " is on the polygon boundary.\n";
        return 1;
      case CGAL::ON_UNBOUNDED_SIDE:
//        cout << " is outside the polygon.\n";
        return 0;
    }
}


#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Alpha_shape_vertex_base_2<EPICK>                              Vba;
typedef CGAL::Alpha_shape_face_base_2<EPICK>                                Fba;
typedef CGAL::Triangulation_data_structure_2<Vba,Fba>                       Tdsa;
typedef CGAL::Alpha_shape_2<CGAL::Delaunay_triangulation_2<EPICK,Tdsa>>     Alpha_shape_2;
//Polygon_2 alphaShapeToPolygon(dirHolder& dir, Alpha_shape_2& as, EPICK::Plane_3& plane){

//    vector<Segment> segments;

//    vector<deque<Alpha_shape_2::Vertex_handle>> chains;
//    deque<Alpha_shape_2::Vertex_handle> pf;

//    Alpha_shape_2::Alpha_shape_edges_iterator eit;
//    Alpha_shape_2::Face_handle face;
//    Alpha_shape_2::Vertex_handle vfirst;
//    Alpha_shape_2::Vertex_handle vsecond;
//    // find first border edge
//    for(eit = as.alpha_shape_edges_begin(); eit != as.alpha_shape_edges_end(); eit++){

//        // see here how to access the two vertex handles of the edge: https://stackoverflow.com/questions/4837179/getting-a-vertex-handle-from-an-edge-iterator
//        face = eit->first;
//        int i = eit->second;
//        vfirst = face->vertex(face->cw(i));
//        vsecond = face->vertex(face->ccw(i));

//        if(as.classify(vfirst) == Alpha_shape_2::REGULAR && as.classify(vsecond) == Alpha_shape_2::REGULAR)
//            break;
//    }
//    pf.push_back(vfirst);
//    pf.push_back(vsecond);
//    chains.push_back(pf);

//    int edge_added;
//    while(++eit != as.alpha_shape_edges_end()){

//        edge_added = 0;

//        face = eit->first;
//        int i = eit->second;
//        vfirst = face->vertex(face->cw(i));
//        vsecond = face->vertex(face->ccw(i));

//        if(as.classify(vfirst) != Alpha_shape_2::REGULAR || as.classify(vsecond) != Alpha_shape_2::REGULAR)
//            continue;

//        for(int m = 0; m < chains.size(); m++){
////        for(auto m : chains){

//            if(vfirst == chains[m].front()){
//                chains[m].push_front(vsecond);
//                edge_added = 1;
//                break;
//            }
//            else if(vfirst == chains[m].back()){
//                chains[m].push_back(vsecond);
//                edge_added = 1;
//                break;
//            }
//            else if(vsecond == chains[m].front()){
//                chains[m].push_front(vfirst);
//                edge_added = 1;
//                break;
//            }
//            else if(vsecond == chains[m].back()){
//                chains[m].push_back(vfirst);
//                edge_added = 1;
//                break;
//            }
//        }
//        if(!edge_added){
//            deque<Alpha_shape_2::Vertex_handle> pn;
//            pn.push_back(vfirst);
//            pn.push_back(vsecond);
//            chains.push_back(pn);
//        }

//    }


//    vector<Polygon_2> polys;

//    auto outer = chains.begin();
//    while(chains.size() > 0){


//        // TODO: removes chains to early, unclear why

//        // remove complete chains
//        if(outer->front() == outer->back()){
//            Polygon_2 poly;
//            for(auto v : *outer){
//                Alpha_shape_2::Vertex_handle vit = v;
//                poly.push_back(vit->point());
//            }
//            polys.push_back(poly);
//            outer = chains.erase(outer);
//        }


//        auto inner = chains.begin();
//        while(inner != chains.end()){

////        for(auto inner : chains){

//            if(*outer == *inner){
//                inner++;
//                continue;
//            }

//            // check back fits to front
//            Alpha_shape_2::Vertex_handle vstart = outer->front();
//            Alpha_shape_2::Vertex_handle vend = outer->back();

//            cout << inner->size() << endl;

//            auto f = inner->front();
//            if(vend == inner->front()){
//                auto it = inner->begin();
//                it++; // skip the first
//                while(it != inner->end())
//                    outer->push_back(*it++);
//                inner = chains.erase(inner);
//            }
//            else if(vend == inner->back()){
//                auto it = inner->end();
//                it--; // skip the first
//                while(it != inner->begin())
//                    outer->push_back(*it--);
//                inner = chains.erase(inner);
//            }
//            else if(vstart == inner->front()){
//                auto it = inner->begin();
//                it++; // skip the first
//                while(it != inner->end())
//                    outer->push_front(*it++);
//                inner = chains.erase(inner);
//            }
//            else if(vstart == inner->back()){
//                auto it = inner->end();
//                it--; // skip the first
//                while(it != inner->begin())
//                    outer->push_front(*it--);
//                inner = chains.erase(inner);
//            }
//            else{
//                inner++;
//            }
////            else if(vfirst == m.back()){
////                m.push_back(vsecond);
////            }
////            else if(vsecond == m.front()){
////                m.push_front(vfirst);
////            }
////            else if(vsecond == m.back()){
////                m.push_back(vfirst);
////            }
//        }
//        // circulate (i.e. start at the beginning again)
//        if(outer == chains.end())
//            outer = chains.begin();
//        outer++;


//    }

////    cout << "poly is simple: " << poly.is_simple() << endl;
////    assert(poly.is_simple());

////    ofstream out(dir.path+dir.write_file+dir.suffix+".off");
////    out << "OFF\n" << vertex_count << " 1 0\n";
////    out << pts_stream.str();
////    out << vertex_count << " " << poly_stream.str();

//    dir.suffix = "_segs" + dir.suffix;
//    exportEdges(dir,segments);

//    return polys[0];

//}






Polygon_2 alphaShapeToPolygon(dirHolder& dir, Alpha_shape_2& as, EPICK::Plane_3& plane){

    vector<Segment> segments;

//    stringstream pts_stream, poly_stream;

    vector<pair<Alpha_shape_2::Face_handle,int>> processed_edges;

    // find first border edge
    Alpha_shape_2::Alpha_shape_edges_iterator eit = as.alpha_shape_edges_begin();
    Alpha_shape_2::Vertex_handle vfirst;
    Alpha_shape_2::Vertex_handle vsecond;
    Alpha_shape_2::Face_handle face;
    int i;
    // collect all border edges (this includes hole borders
    for(eit = as.alpha_shape_edges_begin(); eit != as.alpha_shape_edges_end(); eit++){


        // see here how to access the two vertex handles of the edge: https://stackoverflow.com/questions/4837179/getting-a-vertex-handle-from-an-edge-iterator
        face = eit->first;
        i = eit->second;
        vfirst = face->vertex(face->cw(i));
        vsecond = face->vertex(face->ccw(i));

        if(as.classify(vfirst) == Alpha_shape_2::REGULAR && as.classify(vsecond) == Alpha_shape_2::REGULAR)
            processed_edges.push_back(*eit);
    }

    vector<pair<Polygon_2,double>> polys;
    // while not all border edges are used to assemble the polygon (and its holes) continue this loop
    while(processed_edges.size() > 1){

        // start at a random border edge
        Polygon_2 poly;

        face = processed_edges[0].first;
        i = processed_edges[0].second;
        vfirst = face->vertex(face->cw(i));
        vsecond = face->vertex(face->ccw(i));

        processed_edges.erase(processed_edges.begin());

        set<Alpha_shape_2::Vertex_handle> poly_handles;

//        int vertex_count = 0;
//        pts_stream << plane.to_3d(vsecond->point()) << endl;
//        poly_stream << vertex_count++ << " ";
        poly.push_back(vsecond->point());
        poly_handles.insert(vsecond);

        // remember the first starting point as a break criterion
        auto vvfirst = vfirst;
        while(vsecond != vvfirst){

            // circulate edges around vsecond
            Alpha_shape_2::Edge_circulator ecirc = as.incident_edges(vsecond,face);
            auto fecirc = ecirc;
            while(fecirc != ++ecirc){
                face = ecirc->first;
                int i = ecirc->second;
                vfirst = face->vertex(face->cw(i));
                vsecond = face->vertex(face->ccw(i));
                if(as.classify(vfirst) == Alpha_shape_2::REGULAR && as.classify(vsecond) == Alpha_shape_2::REGULAR){
                    auto current_edge = find(processed_edges.begin(),processed_edges.end(),make_pair(face,i));
                    processed_edges.erase(current_edge);
                    break;
                }
            }

            int size = poly_handles.size();
            poly_handles.insert(vsecond);

            // dirty hack to slightly move the point towards the previous point, if this vertex of the alpha shape already exists in the polygon
            // also the only reason why the poly_handles set is needed
            if(size == poly_handles.size()){
                auto mover = EPICK::Vector_2(vsecond->point(),poly[poly.size()-1]);
                auto vsecond_point = EPICK::Point_2(vsecond->point().x()+mover.x()*0.1,vsecond->point().y()+mover.y()*0.1);
                Segment seg(plane.to_3d(poly[poly.size()-1]), plane.to_3d(vsecond_point));
//                pts_stream << plane.to_3d(vsecond_point) << endl;
                poly.push_back(vsecond_point);
                segments.push_back(seg);
            }
            else{
                Segment seg(plane.to_3d(poly[poly.size()-1]), plane.to_3d(vsecond->point()));
//                pts_stream << plane.to_3d(vsecond->point()) << endl;
                poly.push_back(vsecond->point());
                segments.push_back(seg);
            }

//            poly_stream << vertex_count++ << " ";

            if(processed_edges.size() == 1)
                break;
        }

//        cout << "poly is simple: " << poly.is_simple() << endl;
        if(poly.is_simple()){
            // some holes are not simple, unclear why
//        assert(poly.is_simple());
            polys.push_back(make_pair(poly,abs(poly.area())));
        }
    }
    // poly with the biggest area should be the outer boundary
    // wholes are discarded here but could be add with something like
    // Polygon_with_holes_2(polys[0].first, polys.begin()++.first(), polys.end().first))
    sort(polys.begin(), polys.end(), sort_by_area);

    assert(polys[0].first.is_simple());

//    ofstream out(dir.path+dir.write_file+dir.suffix+".off");
//    out << "OFF\n" << vertex_count << " 1 0\n";
//    out << pts_stream.str();
//    out << vertex_count << " " << poly_stream.str();

//    exportEdges(dir,segments);

    return polys[0].first;

}




void makeAlphaShapes(dirHolder dir, dataHolder& data){

//    vector<shared_ptr<Alpha_shape_2>> alpha_shapes;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nAlpha shaping..."  << endl;
    // cannot actually do this with alpha shape, because _with_info base is not supported for alpha shape
    // thus there is absolutely no way to keep track of the 3D point a 2D point corresponds to.            

    RansacPlaneRange planes = data.efficient_ransac.planes();
    cout << "\nMake 2D convex hulls of inlier points per plane..." << endl;
    cout << "\t-export " << planes.size() << " polygons to " << dir.write_file+"_ransac_polygon_X.off" << endl;
    const int nc=8;
    array<CGAL::Color,nc> color_palette{CGAL::red(), CGAL::blue(), CGAL::green(), CGAL::orange(), CGAL::gray(), CGAL::purple(), CGAL::violet(), CGAL::yellow()};

    // alpha shape the projected inlier points, which should give me a 2-manifold;
    auto it = planes.begin();
    int count = 0;
    vector<Point> as_points;
    vector<EPICK::Point_2> projected_inlier_points;
    int pc = 0;
    while(it != planes.end()){
//        cout << "process plane: " << pc++ << endl;
        DT2 delaunay;
        vector<size_t>::const_iterator index_it = (*it)->indices_of_assigned_points().begin();
        EPICK::Plane_3 current_plane(*it->get());
        projected_inlier_points.clear();
        while (index_it != (*it)->indices_of_assigned_points().end()){
            // Retrieve point.
            const auto pt = *(data.point_set.begin() + (*index_it));
            // insert the point projected to plane and then plane projected to XY into the 2DT
            auto handle = delaunay.insert(current_plane.to_2d(data.point_set.point(pt)));
            // keep the 3D point for exporting
            handle->info() = current_plane.projection(data.point_set.point(pt));
            // insert the point into an alpha shape
            projected_inlier_points.push_back(current_plane.to_2d(data.point_set.point(pt)));
            // Proceed with the next point.
            index_it++;
        }
        // alpha shapping
        Alpha_shape_2 as(projected_inlier_points.begin(), projected_inlier_points.end());
        as.set_mode(Alpha_shape_2::REGULARIZED);
        Alpha_shape_2::Alpha_iterator opt = as.find_optimal_alpha(1);
        as.set_alpha(*opt);

        // stuff for exporting the alpha shape (not converting it to a simple polygon
//        as_points.clear();
//        auto asi = as.all_vertices_begin();
//        while(asi != as.all_vertices_end()){
//            if(as.classify(asi)==Alpha_shape_2::REGULAR || as.classify(asi)==Alpha_shape_2::SINGULAR)
//                as_points.push_back(current_plane.to_3d(asi->point()));
//            asi++;
//        }
//        dir.suffix = "_alpha_shape_" + to_string(count);
//        exportAlphaShape(dir,as_points);
//        data.ass_points.push_back(as_points);

        auto poly = alphaShapeToPolygon(dir,as,current_plane);
        data.alpha_polygons.push_back(poly);
        dir.suffix = "_polygon_" + to_string(count);
        exportPolygon(dir,poly,current_plane);

        // a vector of alpha_shapes in the form of a vector of shared pointers
//        shared_ptr<Alpha_shape_2> ptr2as(&as);
//        alpha_shapes.push_back(ptr2as);

//        dir.suffix = "_ransac_polygon_" + to_string(count);
//        exportConvexHull2d(dir,delaunay,color_palette[count%nc]);
        // Proceed with the next detected shape.
        count++, it++;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done after "<< duration.count() << "s" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// STRUCTURE /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PSS structure(dataHolder& data,
               meshProcessingOptions& options){

    cout << "\nStructure point set..." << endl;
    auto start = std::chrono::high_resolution_clock::now();

//    Efficient_ransac efficient_ransac = *data.ransac;

    auto planes = data.efficient_ransac.planes();
//    auto index_map = CGAL::Shape_detection::Point_to_shape_index_map<RansacTraits>(data.point_set.points(), planes);
    auto index_map = CGAL::Shape_detection::Point_to_shape_index_map<RansacTraits>(data.point_set, planes);

    data.point_set.collect_garbage();

    double structure_epsilon = options.structure_epsilon;
    if(!structure_epsilon > 0.0){
        CGAL::oriented_bounding_box(data.points, data.bounding_box);
        structure_epsilon = 0.01*sqrt(max(CGAL::squared_distance(data.bounding_box[0],data.bounding_box[2]),CGAL::squared_distance(data.bounding_box[0],data.bounding_box[4])));
    }
    cout << "\t-structure epsilon: " << structure_epsilon << endl;


    #define CGAL_PSP3_VERBOSE
    // init calls run() and clean()
    PSS pss(data.point_set,
         planes,
         structure_epsilon,
         CGAL::parameters::point_map(data.point_set.point_map()).
         normal_map(data.point_set.normal_map()).
         plane_map(CGAL::Shape_detection::Plane_map<RansacTraits>()).
         plane_index_map(index_map));


    for(const auto e : pss.m_edges)
        data.pss_segments.push_back(e.segment);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done after "<< duration.count() << "s" << endl;
    return pss;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// OPTIMIZATION //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcUnaries(dataHolder& data, PSS& pss){

    // TODO: what could in theory be done to speed this up even more is to snap the obvious once; like I am doing it in snapping
    // maybe with a very low threshold to only snap very obvious stuff
    // then do an edge collapse - not loosing any information, but getting rid of already obvious snaps
    // only then solve the rest with a graph-cut


    auto start = chrono::high_resolution_clock::now();
    cout << "\t-calculate unaries" << endl;

    data.point_set.collect_garbage();

    RansacPlaneRange planes = data.efficient_ransac.planes();
    auto index_map = CGAL::Shape_detection::Point_to_shape_index_map<RansacTraits>(data.point_set, planes);

    // mesh properties
    SurfaceMesh::Property_map<vertex_descriptor,int> vlabel = data.smesh.add_property_map<vertex_descriptor, int>("v:label", 0).first;
    SurfaceMesh::Property_map<vertex_descriptor,index_dist_map> vproxy = data.smesh.add_property_map<vertex_descriptor, index_dist_map>("v:proxy").first;
    SurfaceMesh::Property_map<vertex_descriptor,int> vpindex = data.smesh.add_property_map<vertex_descriptor, int>("v:pindex").first;
    // point set properties
    Point_set::Property_map<SurfaceMesh::Vertex_index> vpoint = data.point_set.property_map<SurfaceMesh::Vertex_index>("v:index").first;


    vector<pair<int,double>> plane,edge,corner;
    index_dist_map idmap;
    double dist;
    double min_dist;
    int fl = 0, pl = 0, el = 0, cl = 0;
    for(size_t idx = 0; idx < data.point_set.size(); idx++){

        plane.clear(), edge.clear(), corner.clear();

        SurfaceMesh::Vertex_index vidx = vpoint[idx];
        Point pt = data.smesh.point(vidx);
        int plane_index = get(index_map, idx);

        // plane dist
        if(plane_index == -1){
            int i = 0;
            for(const auto p : planes){
                EPICK::Plane_3 dp = *p;
                // if the free form point projects into the alpha shape, take the plane distance
                if(pointInAlphaShape(data.alpha_polygons[i],dp.to_2d(pt)))
                    dist = CGAL::squared_distance(dp, pt);
                // if not, take the distance to the border of the alpha shape as plane distance
                else{
                    min_dist = 1000000; // abuse min_dist to each proxy type as min_dist to all planes for free form points
                    auto vit = data.alpha_polygons[i].begin();
                    while(vit != data.alpha_polygons[i].end()){
                        dist = CGAL::squared_distance(dp.to_3d(*vit),pt);
                        if(dist < min_dist)
                            min_dist = dist;
                        vit++;
                    }
                    dist = min_dist;
                }
                plane.push_back(make_pair(i++,dist));
            }
            sort(plane.begin(), plane.end(), sort_by_dist);
            plane[0].second = sqrt(plane[0].second);
            idmap.plane = plane[0];
            plane_index = plane[0].first;
            min_dist = 0;
            vlabel[vidx] = 0; // init as free form
            fl++;
        }
        // not a free form point, so simply get the distance to the assigned plane
        // could be changed to get the distance to the alpha shape instead (if the point doesn't project into it)
        else{
            dist = CGAL::squared_distance(data.ransac_planes[plane_index],pt);
            idmap.plane = make_pair(plane_index,sqrt(dist));
            min_dist = dist;
            vlabel[vidx] = 1; // init as plane
            pl++;
        }
        vpindex[vidx] = plane_index;

        // edge dist
        for(int i = 0; i < pss.m_edges.size(); i++){
            if(!pss.m_edges[i].active)
                continue;
            // TODO: instead of conditioning on this, condition on whether the current point is in epsilon distance
            // to all alpha shapes (which planes) constitute this edge
            // this is to avoid points snapping to non existing edges (edges that have coplanar adjacent faces)
            auto eit = find(pss.m_edges[i].indices.begin(), pss.m_edges[i].indices.end(), idx);
            if(eit != pss.m_edges[i].indices.end()){
                dist = CGAL::squared_distance(pss.m_edges[i].segment, pt);
                edge.push_back(make_pair(i,dist));
            }
        }
        if(edge.empty()){ // only if no edge found, do exhaustive search
//            for(int i = 0; i < pss.m_edges.size(); i++){
//                dist = CGAL::squared_distance(pss.m_edges[i].support, pt);
//                edge.push_back(make_pair(i,dist));
//            }
            edge.push_back(make_pair(-1,1000));
        }
        sort(edge.begin(), edge.end(), sort_by_dist);
        dist = sqrt(edge[0].second) / 2;
        edge[0].second = dist;
        idmap.edge = edge[0];
        if(dist < min_dist){
            min_dist = dist;
            vlabel[vidx] = 2; // init as edge
            el++;
        }


        // corner dist
        for(int i = 0; i < pss.m_corners.size(); i++){
            if(!pss.m_corners[i].active)
                continue;
            auto cit = find(pss.m_corners[i].planes.begin(), pss.m_corners[i].planes.end(), plane_index);
            if(cit != pss.m_corners[i].planes.end()){
                dist = CGAL::squared_distance(pss.m_corners[i].support, pt);
                corner.push_back(make_pair(i,dist));
            }
        }
        if(corner.empty()){ // only if no corner found, do exhaustive search
            for(int i = 0; i < pss.m_corners.size(); i++){
                dist = CGAL::squared_distance(pss.m_corners[i].support, pt);
                corner.push_back(make_pair(i,dist));
            }
//            corner.push_back(make_pair(-1,1000));
        }
        sort(corner.begin(), corner.end(), sort_by_dist);
        dist = sqrt(corner[0].second) / 3;
        corner[0].second = dist;
        idmap.corner = corner[0];
        if(dist < min_dist){
            vlabel[vidx] = 3; // init as corner
            cl++;
        }

        // save everything in the proxy property which holds all necessary data for setting the unaries in the energy
        vproxy[vidx] = idmap;
//        cout << "plane: " << idmap.plane.second
//             << "; edge: " << idmap.edge.second
//             << "; corner: " << idmap.corner.second
//             << "-> label: " << vlabel[vidx] << endl;
//        cout << "done " << count++ << "/" << data.smesh.number_of_vertices() << endl;
    }
    cout << "\t-initial labels:" << endl;
    cout << "\t\t-free: " << fl << endl;
    cout << "\t\t-plane: " << pl << endl;
    cout << "\t\t-edge: " << el << endl;
    cout << "\t\t-corner: " << cl << endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done after " << duration.count() << "s" << endl;

}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "GCoptimization.h"
typedef double gtype;
void optimizeSnapping(dataHolder& data, PSS& pss, meshProcessingOptions& mo, runningOptions& ro){

    cout << "\nMesh snapping..." << endl;


    // calc unaries
    calcUnaries(data, pss);

    gtype free_cost = ro.free_weight;
    if(ro.free_weight == -1.0){
        free_cost = mo.ransac_epsilon;
    }
    double lambda = ro.reg_weight;
    SurfaceMesh::Property_map<vertex_descriptor,index_dist_map> vproxy = data.smesh.property_map<vertex_descriptor, index_dist_map>("v:proxy").first;
    SurfaceMesh::Property_map<vertex_descriptor,int> vindex = data.smesh.add_property_map<vertex_descriptor, int>("v:index").first;
    SurfaceMesh::Property_map<vertex_descriptor,int> vlabel = data.smesh.property_map<vertex_descriptor, int>("v:label").first;

    int num_vertices = data.smesh.number_of_vertices();
    int num_labels = 4;

    GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_vertices,num_labels);

    // first set up the array for data costs and set the initial label in the same loop
    gtype *data_term = new gtype[num_vertices*num_labels];

    SurfaceMesh::vertex_iterator vit;
    int idx = 0;
    for(vit = data.smesh.vertices().begin(); vit != data.smesh.vertices().end(); vit++){

        data_term[idx*4+0] = free_cost;
        data_term[idx*4+1] = vproxy[*vit].plane.second;
        data_term[idx*4+2] = vproxy[*vit].edge.second;
        data_term[idx*4+3] = vproxy[*vit].corner.second;

        // init vertices to proxy with shortest dist
        gc->setLabel(idx,vlabel[*vit]);

        // set an index for the neighborhood function
        vindex[*vit] = idx++;
    }


    if(!ro.optimization){
        delete gc;
        return;
    }
    auto start = chrono::high_resolution_clock::now();

    // next set up the array for smooth costs
    gtype *smooth_term = new gtype[num_labels*num_labels];
    //free
    smooth_term[0] = 0; //free,free
    smooth_term[1] = 1; //free,plane
    smooth_term[2] = 1; //free,edge
    smooth_term[3] = 1; //free,corner
    //plane
    smooth_term[4] = 1; //plane,free
    smooth_term[5] = 0; //plane,plane
    smooth_term[6] = 0.5; //plane,edge
    smooth_term[7] = 1; //plane,corner
    //edge
    smooth_term[8]  = 1; //edge,free
    smooth_term[9]  = 0.5; //edge,plane
    smooth_term[10] = 0; //edge,edge
    smooth_term[11] = 1; //edge,corner
    //corner
    smooth_term[12] = 1; //corner,free
    smooth_term[13] = 1; //corner,plane
    smooth_term[14] = 1; //corner,edge
    smooth_term[15] = 0; //corner,corner


//    smooth_term[0] = 1; //free,free
//    smooth_term[1] = 2; //free,plane
//    smooth_term[2] = 2; //free,edge
//    smooth_term[3] = 2; //free,corner
//    //plane
//    smooth_term[4] = 2; //plane,free
//    smooth_term[5] = 1; //plane,plane
//    smooth_term[6] = 1.5; //plane,edge
//    smooth_term[7] = 2; //plane,corner
//    //edge
//    smooth_term[8]  = 2; //edge,free
//    smooth_term[9]  = 1.5; //edge,plane
//    smooth_term[10] = 0.5; //edge,edge
//    smooth_term[11] = 1.5; //edge,corner
//    //corner
//    smooth_term[12] = 2; //corner,free
//    smooth_term[13] = 2; //corner,plane
//    smooth_term[14] = 1.5; //corner,edge
//    smooth_term[15] = 0.5; //corner,corner

//    smooth_term[0] = 1; //free,free
//    smooth_term[1] = 2; //free,plane
//    smooth_term[2] = 2.5; //free,edge
//    smooth_term[3] = 3; //free,corner
//    //plane
//    smooth_term[4] = 2; //plane,free
//    smooth_term[5] = 1; //plane,plane
//    smooth_term[6] = 2; //plane,edge
//    smooth_term[7] = 2; //plane,corner
//    //edge
//    smooth_term[8]  = 2.5; //edge,free
//    smooth_term[9]  = 2; //edge,plane
//    smooth_term[10] = 0.5; //edge,edge
//    smooth_term[11] = 2; //edge,corner
//    //corner
//    smooth_term[12] = 3; //corner,free
//    smooth_term[13] = 2; //corner,plane
//    smooth_term[14] = 2; //corner,edge
//    smooth_term[15] = 1; //corner,corner

    // TODO: need to find a way to make edge snapping cheap,
    // but not cheap enough so that plane points snap to no valid edges

    gc->setDataCost(data_term);
    gc->setSmoothCost(smooth_term);

    // set neighbors
    for(vit = data.smesh.vertices().begin(); vit != data.smesh.vertices().end(); vit++){

        auto first_hedge = data.smesh.halfedge(*vit);
        auto next_hedge = first_hedge;
        auto target = data.smesh.target(first_hedge);
        do{
            auto source = data.smesh.source(next_hedge);

//            double weight = -1;
            double weight = 1;
            if(vproxy[source].plane.first != vproxy[target].plane.first)
                weight = 1;
            if(vindex[source] > vindex[target])
                gc->setNeighbors(vindex[source], vindex[target], lambda*weight);
            next_hedge = data.smesh.next_around_target(next_hedge);

        }while(first_hedge != next_hedge);
    }



    // Optimization
    cout << "\t-optimize with free-form weight, lambda: " << free_cost << "," << ro.reg_weight << endl;
    cout << "\t-before optimization data + smoothness energy: " << endl;
    cout << "\t\t" << gc->giveDataEnergy() << " + " << gc->giveSmoothEnergy() << " = "  << gc->compute_energy() << endl;
    // use swap because it is good when you have two labels
//        gc->swap(-1);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
    gc->expansion(5);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
    // iterate over all cells and get their optimized label
    for(vit = data.smesh.vertices().begin(); vit != data.smesh.vertices().end(); vit++)
        vlabel[*vit] = gc->whatLabel(vindex[*vit]);
    cout << "\t-after optimization energy is " << gc->compute_energy() << endl;


    delete gc;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;

    // TODO: only thing left to do is iterate over the vertices again and snap them to their corresponding label
    // i.e. if label = 1, snap them to vproxy[vidx].plane.first (maybe make a EPICK::Plane_3 vector to make things easier, because there is the iterator index problem with planes again)
    // when label = 2, snap to pss.m_edges[vproxy[vidx].edge.first]
    // and so on...



}

void snapAndCollapse(dataHolder& data, PSS& pss){


    SurfaceMesh::Property_map<vertex_descriptor,index_dist_map> vproxy = data.smesh.property_map<vertex_descriptor, index_dist_map>("v:proxy").first;
    SurfaceMesh::Property_map<vertex_descriptor,int> vlabel = data.smesh.property_map<vertex_descriptor, int>("v:label").first;
    SurfaceMesh::Property_map<vertex_descriptor,CGAL::Color> vcolor = data.smesh.add_property_map<vertex_descriptor, CGAL::Color>("v:color").first;
    SurfaceMesh::Property_map<vertex_descriptor,Vertex_status> vstatus = data.smesh.add_property_map<vertex_descriptor, Vertex_status>("v:status", FREE).first;

    SurfaceMesh::vertex_iterator vit;
    for(vit = data.smesh.vertices().begin(); vit != data.smesh.vertices().end(); vit++){

        // free form
        if(vlabel[*vit] == 0){
            vcolor[*vit] = CGAL::white();
            vstatus[*vit] = FREE;
            continue;
        }

        // plane
        if(vlabel[*vit] == 1){
//            data.smesh.point(*vit) = data.ransac_planes[vproxy[*vit].plane.first].projection(data.smesh.point(*vit));
            vcolor[*vit] = CGAL::blue();
            vstatus[*vit] = PLANE;
            continue;
        }

        // edge
        if(vlabel[*vit] == 2){
//            data.smesh.point(*vit) = pss.m_edges[vproxy[*vit].edge.first].support.projection(data.smesh.point(*vit));
            vcolor[*vit] = CGAL::red();
            vstatus[*vit] = EDGE;
            continue;
        }

        // corner
        if(vlabel[*vit] == 3){
//            data.smesh.point(*vit) = pss.m_corners[vproxy[*vit].corner.first].support;
            vcolor[*vit] = CGAL::yellow();
            vstatus[*vit] = CORNER;
            continue;
        }
    }


}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// SNAPPING ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void snapWithLabel(dataHolder& data, PSS& pss){

    SurfaceMesh::Property_map<vertex_descriptor,index_dist_map> vproxy = data.smesh.property_map<vertex_descriptor, index_dist_map>("v:proxy").first;
    SurfaceMesh::Property_map<vertex_descriptor,int> vlabel = data.smesh.property_map<vertex_descriptor, int>("v:label").first;
    SurfaceMesh::Property_map<vertex_descriptor,CGAL::Color> vcolor = data.smesh.add_property_map<vertex_descriptor, CGAL::Color>("v:color").first;
    SurfaceMesh::Property_map<vertex_descriptor,Vertex_status> vstatus = data.smesh.add_property_map<vertex_descriptor, Vertex_status>("v:status", FREE).first;

    SurfaceMesh::vertex_iterator vit;
    for(vit = data.smesh.vertices().begin(); vit != data.smesh.vertices().end(); vit++){

        // free form
        if(vlabel[*vit] == 0){
            vcolor[*vit] = CGAL::white();
            vstatus[*vit] = FREE;
            continue;
        }

        // plane
        if(vlabel[*vit] == 1){
            data.smesh.point(*vit) = data.ransac_planes[vproxy[*vit].plane.first].projection(data.smesh.point(*vit));
            vcolor[*vit] = CGAL::blue();
            vstatus[*vit] = PLANE;
            continue;
        }

        // edge
        if(vlabel[*vit] == 2){
            data.smesh.point(*vit) = pss.m_edges[vproxy[*vit].edge.first].support.projection(data.smesh.point(*vit));
            vcolor[*vit] = CGAL::red();
            vstatus[*vit] = EDGE;
            continue;
        }

        // corner
        if(vlabel[*vit] == 3){
            data.smesh.point(*vit) = pss.m_corners[vproxy[*vit].corner.first].support;
            vcolor[*vit] = CGAL::yellow();
            vstatus[*vit] = CORNER;
            continue;
        }
    }
}
void snapMesh(dataHolder& data,
                PSS& pss,
                meshProcessingOptions& options){

    auto start = std::chrono::high_resolution_clock::now();

    cout << "\nSnap mesh..." << endl;

    data.point_set.collect_garbage();

    RansacPlaneRange planes = data.efficient_ransac.planes();
    auto index_map = CGAL::Shape_detection::Point_to_shape_index_map<RansacTraits>(data.point_set, planes);

    double snap_epsilon = options.snap_epsilon;
    if(!snap_epsilon > 0.0){
        CGAL::oriented_bounding_box(data.smesh, data.bounding_box);
        snap_epsilon = 0.01*sqrt(max(CGAL::squared_distance(data.bounding_box[0],data.bounding_box[2]),CGAL::squared_distance(data.bounding_box[0],data.bounding_box[4])));
    }
    cout << "\t-snap epsilon: " << snap_epsilon << endl;


    // get corners and edges calculated by pss
    auto edges = pss.m_edges;
    auto corners = pss.m_corners;

    double dist;
    bool snapped;

    SurfaceMesh::Property_map<vertex_descriptor,CGAL::Color> vcolors = data.smesh.add_property_map<vertex_descriptor, CGAL::Color>("v:color", CGAL::white()).first;
    SurfaceMesh::Property_map<vertex_descriptor,Vertex_status> vstatus = data.smesh.add_property_map<vertex_descriptor, Vertex_status>("v:status", FREE).first;
    SurfaceMesh::Property_map<vertex_descriptor,int> vpindex = data.smesh.add_property_map<vertex_descriptor, int>("v:pindex").first;
    Point_set::Property_map<SurfaceMesh::Vertex_index> vpoint  = data.point_set.property_map<SurfaceMesh::Vertex_index>("v:index").first;


    for(size_t idx = 0; idx < data.point_set.size(); idx++){
        snapped = false;

        SurfaceMesh::Vertex_index vidx = vpoint[idx];
        Point pt = data.smesh.point(vidx);
        int plane_index = get (index_map, idx);
        vpindex[vidx] = plane_index;

        // first check if it is an unassigned point, and if so push it back and continue
        if(plane_index == -1){
            vcolors[vidx] = CGAL::white();
            vstatus[vidx] = FREE;
            continue;
        }


        // check corner snap
        for(const auto c : corners){
            if(!c.active)
                continue;
            dist = sqrt(CGAL::squared_distance(c.support, pt));
            if(dist <= sqrt(3)*snap_epsilon){
                vcolors[vidx] = CGAL::yellow();
                data.smesh.point(vidx) = c.support;
                vstatus[vidx] = CORNER;
                snapped = true;
            }
            if(snapped)
                break;
        }
        if(snapped) // snapped to corner
            continue;

        // check edge snap
        for(const auto e : edges){
            if(!e.active)
                continue;
            dist = sqrt(CGAL::squared_distance(e.support, pt));
            if(dist <= sqrt(2)*snap_epsilon){
                // here I am just checking distance to a line (infinite),
                // what I would need to check is distance to segment, but there is no such thing in PSS
                // however, the indices vector of edges has the correct assignments, so can simply check there
                auto eit = find(e.indices.begin(), e.indices.end(), idx);
                if(eit != e.indices.end()){
                    vcolors[vidx] = CGAL::red();
                    data.smesh.point(vidx) = e.support.projection(pt);
                    vstatus[vidx] = EDGE;
                    snapped = true;
                }
            }
            if(snapped)
                break;
        }
        if(snapped) // snapped to line
            continue;

        // check plane snap to correct plane
        int i = 0;
        for(const auto pl : planes){
            if(i++ == plane_index){
                EPICK::Plane_3 dp = *pl;
                // here I am just checking distance to a plane (infinite),
                // what I would need to check is distance to polygon, but there is no such thing in PSS
                vcolors[vidx] = CGAL::blue();
                data.smesh.point(vidx) = dp.projection(pt);
                vstatus[vidx] = PLANE;
                snapped = true;
            }
            if(snapped)
                break;
        }
        if(snapped) // snapped to correct plane
            continue;

        // check plane snap to any plane within dist
        for(const auto pl : planes){
            EPICK::Plane_3 dp = *pl;
            dist = sqrt(CGAL::squared_distance(dp, pt));
            if(dist <= snap_epsilon){
                // here I am just checking distance to a plane (infinite),
                // what I would need to check is distance to polygon, but there is no such thing in PSS
                vcolors[vidx] = CGAL::blue();
                data.smesh.point(vidx) = dp.projection(pt);
                vstatus[vidx] = PLANE;
                snapped = true;
            }
            if(snapped)
                break;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done after "<< duration.count() << "s" << endl;

}
void snapPoint_set(dirHolder& dir, dataHolder& data,
                   PSS& pss,
          meshProcessingOptions& options){

    cout << "\nSnap point set..." << endl;

    data.point_set.collect_garbage();

    RansacPlaneRange planes = data.efficient_ransac.planes();
    auto index_map = CGAL::Shape_detection::Point_to_shape_index_map<RansacTraits>(data.point_set, planes);

    double snap_epsilon = options.snap_epsilon;
    if(!snap_epsilon > 0.0){
        CGAL::oriented_bounding_box(data.points, data.bounding_box);
        snap_epsilon = 0.01*sqrt(max(CGAL::squared_distance(data.bounding_box[0],data.bounding_box[2]),CGAL::squared_distance(data.bounding_box[0],data.bounding_box[4])));
    }
    cout << "\t-snap epsilon: " << snap_epsilon << endl;


    // get corners and edges calculated by pss
    auto edges = pss.m_edges;
    auto corners = pss.m_corners;

//    // sanity check for points in original point set being in epsilon to at least one of the two planes forming an edge
//    vector<Point> edge_points;
//    for(const auto e : edges){
//        if(!e.active)
//            continue;
//        for(const auto idx : e.indices){
//            if(idx < data.point_set.size()){
//                edge_points.push_back(data.point_set.point(idx));
//            }
//        };
//    }
//    dir.suffix = "_edge_points";
//    exportPLY(dir,edge_points);

    double dist;
    bool snapped;
    vector<Point> pcorner, pedge, pplane, pfree, pall;
    vector<vertex_info> icorner, iedge, iplane, ifree, iall;
    Point_set::Property_map<Vector> sensor_vecs = data.point_set.property_map<Vector>("sensor_vec").first;
    Point_set::Property_map<vector<Point>> sensor_positions  = data.point_set.property_map<vector<Point>>("sensor_positions").first;

    for(size_t idx = 0; idx < data.point_set.size(); idx++){
        const auto pt = data.point_set.point(idx);
        vertex_info inf;
        snapped = false;

        // first check if it is an unassigned point, and if so push it back and continue
        int plane_index = get (index_map, idx);
        if (plane_index == -1){
            inf.color = CGAL::white();
            inf.sensor_positions.push_back(sensor_positions[idx][0]);
            inf.sensor_vec = sensor_vecs[idx];

            pall.push_back(pt);
            iall.push_back(inf);
            pfree.push_back(pt);
            ifree.push_back(inf);
            continue;
        }

        // check corner snap
        for(const auto c : corners){
            if(!c.active)
                continue;
            dist = sqrt(CGAL::squared_distance(c.support, pt));
            if(dist <= sqrt(3)*snap_epsilon){
                inf.color = CGAL::yellow();
                inf.sensor_positions.push_back(sensor_positions[idx][0]);
                inf.sensor_vec = sensor_vecs[idx];

                pall.push_back(c.support);
                iall.push_back(inf);
                pcorner.push_back(c.support);
                icorner.push_back(inf);
                snapped = true;

            }
            if(snapped)
                break;
        }
        if(snapped) // snapped to corner
            continue;

        // check edge snap
        for(const auto e : edges){
            if(!e.active)
                continue;
            dist = sqrt(CGAL::squared_distance(e.support, pt));
            if(dist <= sqrt(2)*snap_epsilon){
                // here I am just checking distance to a line (infinite),
                // what I would need to check is distance to segment, but there is no such thing in PSS
                // however, the indices vector of edges has the correct assignments, so can simply check there
                auto eit = find(e.indices.begin(), e.indices.end(), idx);
                if(eit != e.indices.end()){
                    inf.color = CGAL::red();
                    inf.sensor_positions.push_back(sensor_positions[idx][0]);
                    inf.sensor_vec = sensor_vecs[idx];

                    pall.push_back(e.support.projection(pt));
                    iall.push_back(inf);
                    pedge.push_back(e.support.projection(pt));
                    iedge.push_back(inf);
                    snapped = true;
                }
            }
            if(snapped)
                break;
        }
        if(snapped) // snapped to line
            continue;

        // check snap to correct plane
        int i = 0;
        for(const auto pl : planes){
            if(i++ == plane_index){
                EPICK::Plane_3 dp = *pl;
                // here I am just checking distance to a plane (infinite),
                // what I would need to check is distance to polygon, but there is no such thing in PSS
                inf.color = CGAL::blue();
                inf.sensor_positions.push_back(sensor_positions[idx][0]);
                inf.sensor_vec = sensor_vecs[idx];

                pall.push_back(dp.projection(pt));
                iall.push_back(inf);
                pplane.push_back(dp.projection(pt));
                iplane.push_back(inf);
                snapped = true;
            }
            if(snapped)
                break;
        }
        if(snapped) // snapped to correct plane
            continue;

        // check snap to any plane
        for(const auto pl : planes){
            EPICK::Plane_3 dp = *pl;
            dist = sqrt(CGAL::squared_distance(dp, pt));
            if(dist <= snap_epsilon){
                // here I am just checking distance to a plane (infinite),
                // what I would need to check is distance to polygon, but there is no such thing in PSS
                inf.color = CGAL::blue();
                inf.sensor_positions.push_back(sensor_positions[idx][0]);
                inf.sensor_vec = sensor_vecs[idx];

                pall.push_back(dp.projection(pt));
                iall.push_back(inf);
                pplane.push_back(dp.projection(pt));
                iplane.push_back(inf);
                snapped = true;
            }
            if(snapped)
                break;
        }
    }

    exportOptions eo;
    eo.color = true;
    eo.sensor_position = true;
    dir.suffix = "_snapped";
    exportPLY(dir,pall, iall, eo);

    dir.suffix = "_snapped_corners";
    exportPLY(dir,pcorner,icorner,eo);

    dir.suffix = "_snapped_edges";
    exportPLY(dir,pedge,iedge,eo);

    dir.suffix = "_snapped_planes";
    exportPLY(dir,pplane,iplane,eo);

    dir.suffix = "_snapped_free";
    exportPLY(dir,pfree,ifree,eo);

    data.points.clear();
    data.infos.clear();
    data.points = pall;
    data.infos = iall;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// SIMPLIFICATION /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int removeSelfIntersections(dataHolder& data){


    auto fccmap =  data.smesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
    size_t components = CGAL::Polygon_mesh_processing::connected_components(data.smesh, fccmap);

    cout << "\t-number of components: " << components << endl;


    cout << "\nRemove self intersections..." << endl;
    data.smesh.collect_garbage();

    int self_intersect = PMP::does_self_intersect(data.smesh);
    cout << "\t-is closed? " << CGAL::is_closed(data.smesh) << endl;
    cout << "\t-does self intersect? " << self_intersect << endl;
    if(self_intersect){
        cout << "\t-number of mesh vertices " << data.smesh.number_of_vertices() << endl;
        PMP::experimental::remove_self_intersections(data.smesh,
                                                     PMP::parameters::preserve_genus(false).
                                                     apply_per_connected_component(true));
        self_intersect = PMP::does_self_intersect(data.smesh);
        cout << "\t-self intersections removed? " << !self_intersect << endl;
        cout << "\t-number of mesh vertices " << data.smesh.number_of_vertices() << endl;
        cout << "\t-is closed? " << CGAL::is_closed(data.smesh) << endl;

    }
    data.smesh.collect_garbage();
    return self_intersect;

}


// BGL property map which indicates whether an edge is marked as non-removable
struct Constrained_edge_map
  : public boost::put_get_helper<bool, Constrained_edge_map>
{
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;
  Constrained_edge_map(const CGAL::Unique_hash_map<key_type,bool>& aConstraints)
    : mConstraints(aConstraints)
  {}
  reference operator[](const key_type& e) const { return  is_constrained(e); }
  bool is_constrained(const key_type& e) const { return mConstraints.is_defined(e); }
private:
  const CGAL::Unique_hash_map<key_type,bool>& mConstraints;
};
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
void constrainedEdgeCollapse(dirHolder& dir, dataHolder& data){


    cout << "\nConstrained edge collapse..." << endl;
    cout << "\t-edges: " << data.smesh.number_of_edges() << endl;

    // get the surface mesh properties
    auto vstatus = data.smesh.property_map<vertex_descriptor, Vertex_status>("v:status").first;
    auto vpindex = data.smesh.property_map<vertex_descriptor, int>("v:pindex").first;
    auto vcolor = data.smesh.property_map<vertex_descriptor, CGAL::Color>("v:color").first;
    auto vproxy = data.smesh.property_map<vertex_descriptor, index_dist_map>("v:proxy").first;


    // collect constrained edges
    CGAL::Unique_hash_map<edge_descriptor, bool> constraint_hmap(false);
    Constrained_edge_map constraints_map(constraint_hmap);
    CGAL::Unique_hash_map<edge_descriptor, bool> crease_constraint_hmap(false);
    Constrained_edge_map crease_constraints_map(crease_constraint_hmap);


    // cost and placement policies
    // Garland&Heckbert simplification policies
    typedef typename SMS::GarlandHeckbert_policies<SurfaceMesh, EPICK>            GH_policies;
    typedef typename GH_policies::Get_cost                                        GH_cost;
    typedef typename GH_policies::Get_placement                                   GH_placement;
    GH_policies gh_policies(data.smesh);
//    const GH_cost& cost = gh_policies.get_cost();
//    const GH_placement& placement = gh_policies.get_placement();

//    // constrained placement; i.e. vertex for constrained edge stays in place, and thus contracts vertex of non-constrained edge
    const GH_cost& cost = gh_policies.get_cost();
    SMS::Constrained_placement<SMS::Midpoint_placement<SurfaceMesh>,
                               Constrained_edge_map > placement(constraints_map);
//     edle length cost and midpoint placement
//    auto cost = SMS::Edge_length_cost<SurfaceMesh>();
//    auto placement = SMS::Midpoint_placement<SurfaceMesh>();


    std::size_t nb_constrained_edges = 0;
    // gather non collapsable edges
    vector<edge_descriptor> removable_edges;
    for(edge_descriptor ed : edges(data.smesh)){

        halfedge_descriptor hedge = halfedge(ed, data.smesh);
        auto target = data.smesh.target(hedge);
        auto source = data.smesh.source(hedge);
        // detect collapsable edges
        if((vstatus[target] == PLANE &&
                vstatus[source] == PLANE &&
                vpindex[target] == vpindex[source])){
            // if source and target are both snapped to the same plane, edge can be collapsed
            vcolor[target] = vcolor[source] = CGAL::green();
            removable_edges.push_back(ed);
        }
        // maybe this part needs to go in a separate (subsequent) snapping
        else if(vstatus[target] == EDGE &&
                            vstatus[source] == EDGE &&
                            vproxy[source].edge.first == vproxy[target].edge.first){
            // if source and target are both snapped to the same crease, edge can be collapsed
            vcolor[target] = vcolor[source] = CGAL::orange();
            removable_edges.push_back(ed);
        }
        else if(vstatus[target] == CORNER &&
                            vstatus[source] == CORNER &&
                            vproxy[source].corner.first == vproxy[target].corner.first){
            // if source and target are both snapped to a corner, edge can be collapsed
            vcolor[target] = vcolor[source] = CGAL::purple();
            removable_edges.push_back(ed);
        }
        else{
            ++nb_constrained_edges;
            constraint_hmap[ed] = true;
        }
    }
    dir.suffix = "_collapsable";
    exportPLY(dir,data.smesh);
    cout << endl;

//    int i = 0;
//    for(auto ed : removable_edges){
//        cout << "try to collapse: " << ed << endl;
//        if(CGAL::Euler::does_satisfy_link_condition(ed,data.smesh)){
//            CGAL::Euler::collapse_edge(ed,data.smesh);
//            data.smesh.collect_garbage();
//        }
//    }


//    int border_edges = 0;
//    for(auto e : data.smesh.edges()){
//        if(data.smesh.is_border(e))
//            border_edges++;
//    }
//    cout << "\t-mesh border edges: " << border_edges << endl;

    // Contract the surface mesh as much as possible
    SMS::Count_stop_predicate<SurfaceMesh> stop(0);
    cout << "\t-collapsing as many non-feature edges of mesh as possible..." << endl;
    cout << "\t-plane removable: " << data.smesh.number_of_edges() - nb_constrained_edges << endl;
    int r = SMS::edge_collapse(data.smesh, stop,
                               CGAL::parameters::edge_is_constrained_map(constraints_map)
                                                .get_cost(cost)
                                                .get_placement(placement));
    cout << "\t-removed: " << r << endl;
    cout << "\t-edges: "<< data.smesh.number_of_edges() << endl;

    dir.suffix = "_constrained_collapsed";
    exportPLY(dir,data.smesh);

//    std::cout << "\t-check that no removable edge has been forgotten..." << std::endl;
//    r = SMS::edge_collapse(data.smesh,
//                           stop,
//                           CGAL::parameters::edge_is_constrained_map(constraints_map)
//                                            .get_cost(SMS::Edge_length_cost<SurfaceMesh>())
//                                            .get_placement(SMS::Midpoint_placement<SurfaceMesh>()));
//    assert(r == 0);
//    if(r == 0) {
//      std::cout << "\t-OK\n";
//    } else {
//      std::cout << "ERROR! " << r << " edges removed!\n";
//    }

}


#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>

void removeInnerVertices(dirHolder& dir, dataHolder& data){


    cout << "\nSimplify mesh by removing inner vertices..." << endl;

    auto vstatus = data.smesh.add_property_map<vertex_descriptor, Vertex_status>("v:status").first;
    auto vpindex = data.smesh.add_property_map<vertex_descriptor, int>("v:pindex").first;
    auto vcolor = data.smesh.property_map<vertex_descriptor, CGAL::Color>("v:color").first;

    cout << "point set size " << data.point_set.size() << endl;
    cout << "mesh vertices " << data.smesh.number_of_vertices() << endl;

    set<SurfaceMesh::Halfedge_index> removals;
    bool remove;
    auto vh = data.smesh.vertices().begin();
    while(vh != data.smesh.vertices().end()){

        if(vstatus[*vh] != PLANE){
            vh++;
            continue;
        }

        auto first_hedge = data.smesh.halfedge(*vh);
        const auto target = data.smesh.target(first_hedge);
        auto next_hedge = first_hedge;
        remove = true;
        while(remove){
            const auto source = data.smesh.source(next_hedge);
            if(vpindex[source] != vpindex[target])
                remove = false;
            next_hedge = data.smesh.next_around_target(next_hedge);
            if(next_hedge == first_hedge)
                break;
        }
        if(remove){
            data.smesh.remove_vertex(*vh);
            removals.insert(first_hedge);
            vcolor[*vh] = CGAL::green();
        }
        vh++;
    }

////    cout << "\t-collected " << data.smesh.number_of_removed_vertices() << " vertices to remove" << endl;
    int rc = removals.size();
    cout << "\t-collected " << rc << " vertices to remove" << endl;

    vector<SurfaceMesh::Halfedge_index> rev;
    for(auto r : removals)
        rev.push_back(r);

    dir.suffix = "_before_simplification";
    exportOFF(dir,data.smesh);

    data.smesh.remove_property_map(vcolor);
    data.smesh.remove_property_map(vpindex);
    data.smesh.remove_property_map(vstatus);
    data.smesh.collect_garbage();

    SurfaceMesh soup;
    vector<Point> points;
    for(const auto p : data.smesh.points())
        points.push_back(p);
    vector<vector<int>> faces;
    for(auto f : data.smesh.faces()){
        vector<int> face;
        auto v = data.smesh.vertices_around_face(data.smesh.halfedge(f)).begin();
        if(*v > points.size())
            continue;
        face.push_back(*v), v++;
        if(*v > points.size())
            continue;
        face.push_back(*v), v++;
        if(*v > points.size())
            continue;
        face.push_back(*v);
        faces.push_back(face);
    }

    bool oriented = PMP::orient_polygon_soup(points, faces);
    cout << "\t-is oriented: " << oriented << endl;

    cout << "\t-is manifold: " << PMP::is_polygon_soup_a_polygon_mesh(faces) << endl;

    PMP::polygon_soup_to_polygon_mesh(points,faces,soup);
//    cout << "soup to mesh successful? " << soup_worked << endl;
    dir.suffix = "_soup";
    exportOFF(dir,soup);

    SurfaceMesh new_mesh = data.smesh;

//    auto ncolor = new_mesh.add_property_map<vertex_descriptor, CGAL::Color>("v:color", CGAL::white()).first;
//    ncolor[data.smesh.target(rev[77])] = CGAL::green();
//    dir.suffix = "_tester";
//    exportOFF(dir,new_mesh);

    // from: https://stackoverflow.com/questions/60161225/how-do-i-remove-many-vertices-from-a-surface-mesh-in-cgal
    for(int i = 0; i < 163; i++){
        cout << "try edge " << i << "/" << rc;
        if(new_mesh.is_valid(rev[i])){
            cout << ": removed" << endl;
            CGAL::Euler::remove_center_vertex(rev[i],new_mesh);
        }
    }

    data.smesh.collect_garbage();
    dir.suffix = "_after_simplification";
    exportOFF(dir,new_mesh);

}













///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////// CONTROL FUNCTIONS /////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void simplify(dirHolder& dir, dataHolder& data, runningOptions ro, meshProcessingOptions options){

    // fix the random state to get deterministic RANSAC
    CGAL::get_default_random() = CGAL::Random(42);

    ////////////////////////////////////////////
    //////// NORMALS (needed for ransac) ///////
    ////////////////////////////////////////////
//    // Add normal property and estimate normal values
//    cout << "\nEstimate normals..." << endl;
//    data.point_set.add_normal_map();
//    int nb_neighbors = CGAL::estimate_global_k_neighbor_scale(data.point_set);
//    cout << "\t-estimated neighborhood size to " << nb_neighbors << endl;
//    CGAL::pca_estimate_normals<Concurrency_tag>(data.point_set, nb_neighbors, data.point_set.parameters());
//    // TODO: check which normals are actually used for mesh straightening, as
//    // I already put average face normals on the point set before this

    ///////////////////////
    //////// RANSAC ///////
    ///////////////////////
    ransac(dir, data, options);
    dir.suffix = "_ransac";
    exportPLY(dir,data.point_set);
//    exportRansacPlanes(dir,data);
    makeAlphaShapes(dir,data);
    dir.suffix = "_segmentation.seg";
    exportSegmentation(dir,data);

    ///////////////////////
    ////// STRUCTURE //////
    ///////////////////////
    auto pss = structure(data, options);
    exportStructured(dir, pss);
//    lineToSegment(data, pss);
    dir.suffix = "_segments";
    exportEdges(dir,data.pss_segments);

    ///////////////////////
    //////// SNAP  ////////
    ///////////////////////


    if(data.smesh.is_empty()){
        snapPoint_set(dir, data, pss, options);
    }
    else{
//        if(ro.optimization){
            optimizeSnapping(data,pss,options,ro);
            snapAndCollapse(data, pss); // so far this is just putting the label for the collapse
            dir.suffix = "_colored_by_label";
            exportPLY(dir,data.smesh);
//            snapWithLabel(data, pss);
//            if(ro.optimization)
//                dir.suffix = "_optimized";
//            else
//                dir.suffix = "_initial";
//            exportPLY(dir,data.smesh);
    }

    ////////////////////////////////////////////
    //////// REMOVE SELF-INTERSECTIONS /////////
    ////////////////////////////////////////////
//    if(removeSelfIntersections(data)){
//        cout << "WARNING, MESH HAS SELF-INTERSECTIONS" << endl;
//    }
//    else{
//        dir.suffix = "_fixed_self_intersections";
//        exportOFF(dir,data.smesh);
//    }

    ///////////////////////
    ////// SIMPLIFY ///////
    ///////////////////////
    if(options.simplify && !data.smesh.is_empty()){
        //    simplify(dir,data);
        constrainedEdgeCollapse(dir,data);
    }


    snapWithLabel(data, pss);
    if(ro.optimization)
        dir.suffix = "_optimized";
    else
        dir.suffix = "_initial";
    exportPLY(dir,data.smesh);

}


// TODO:
// 1) make a default constructor (without arguments) for PSS, so I can put it on the dataHolder
// 2) read building reconstruction literature a la, Kada, where they probably naively try to find polygon intersections
// hopefully they only fail because of incomplete lines, which is not really a problem in my approach
// 3) implement a polygon-polygon intersection (https://stackoverflow.com/questions/6195413/intersection-between-3d-flat-polygons)
// should be doable. Then I have proper segments, which are under- instead of overcomplete, compared to the
// ones I have right now.
// 4) benchmark my existing method; edges are not important; simplification variable is no of ransac planes (controlled by threshold)
// Then run constrained edge collapse and see if it is better than normal edge collapse without plane detection
// 5) maybe I can remove center vertex of lines; I should still have access to the points that snapped to one line
// 6) can we learn a network to edge-collapse a noisy mesh into a cleaned (noise reduced) simplified mesh by supervising with
// distance to ground truth mesh and in latter epochs minize(number_of_vertices)
// MeshCNN probably does not train their model with such a task; but big parts
// of their architecture could be reused.
// 7) Do a PCA of the vertex set to determine local distance to plane; then use this as a local free form weight


// NOTES:
// illustration of an edge collapse that needs edge flip:
// slide 44 here: http://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/08_Simplification.pdf


int main(int argc, char const *argv[]){

    cliParser ip("simplify");
    if(ip.parse(argc, argv))
        return 1;
    ip.getInput();
    ip.getSimplify();
    ip.getOutput();


    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----MESH SIMPLIFICATION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    meshProcessingOptions options;
    options.max_number_of_proxies = ip.ro.max_number_of_proxies;
    options.number_of_components_to_keep = ip.ro.number_of_components_to_keep;
    options.try_to_close = ip.ro.try_to_close;
    options.ransac_epsilon = ip.ro.ransac_epsilon;
    options.ransac_normal = ip.ro.ransac_normal;
    options.ransac_regularization = ip.ro.ransac_regularization;
    options.structure_epsilon = ip.ro.structure_epsilon;
    options.snap_epsilon = ip.ro.snap_epsilon;
    options.simplify = ip.ro.simplify;

    // read file depending on file type
    if(ip.dh.read_file_type == "off")
        importOff(ip.dh, data.smesh);
    else if(ip.dh.read_file_type == "ply"){
        importLidarMesh(ip.dh, data);
        if(data.facets.size() > 0){
            createSurfaceMesh(data, options);
            createVertexSet(data);
//            data.points.clear();
        }
        else
            createPointSet(data);
    }
    else{
        cerr << "Could not read file ending or file type not supported. Supported are .off and .ply files" << endl;
        return 1;
    }

    ip.dh.suffix = "_input";
    exportOFF(ip.dh,data.smesh);

    simplify(ip.dh, data, ip.ro, options);


    /////////////////////////////
    //////// EVALUATION  ////////
    /////////////////////////////
    if(ip.ro.evaluate_mesh){

        if(!ip.dh.gt_poly_file.empty()){
            if(ip.dh.gt_poly_file.substr(ip.dh.gt_poly_file.length() - 3 ) == "off"){
                if(importOff(ip.dh.path+ip.dh.gt_poly_file, data.gt_poly))
                    return 1;
            }
            else{
                cerr << "\not a valid ending for a ground truth file. choose either .mpu or .off" << endl;
                return 1;
            }
            ip.ro.ground_truth = 1;
        }
        if(ip.ro.ground_truth){
            double iou;
            calcIOU(data, iou);
            printMeshEvaluation(ip.dh,data,iou);
        }
        else
            printMeshEvaluation(ip.dh,data,-1.0);

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----MESH SIMPLIFICATION FINISHED in "<< duration.count() << "s -----\n" << endl;


    return 0;
}



