#include <base/cgal_typedefs.h>

#include <CGAL/IO/read_ply_points.h> // CGAL file that is not present in CGAL version < 4.11, so copied it to my IO folder

#include <IO/fileIO.h>
//#include <util/geometricOperations.h>
//#include <util/helper.h>
#include <CGAL/Point_set_3/IO.h>

#include <boost/filesystem.hpp>

#include <rPLY/ply.h>

#ifdef RECONBENCH
#include "modeling/shape_loader.h"
#endif



using namespace std;




//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////
#ifdef RECONBENCH
int importImplicit(dirHolder dir, dataHolder& data){

    ImplicitFunction* shape = ShapeLoader::load_shape(dir.path+dir.gt_poly_file);
    data.implicit = shape;

    return 0;


}
#endif


#include <CGAL/Orthogonal_incremental_neighbor_search.h>

typedef boost::tuple<Point, vertex_info>                            point_info;

typedef CGAL::Search_traits_3<EPICK> SearchTraits;
typedef CGAL::Search_traits_adapter<point_info,
  CGAL::Nth_of_tuple_property_map<0, point_info>,
  SearchTraits>                                               SearchTraitsAdapted;

typedef CGAL::Orthogonal_incremental_neighbor_search<SearchTraitsAdapted> Incremental_neighbor_search;
typedef Incremental_neighbor_search::Tree Incremental_Tree;
void concatenateData(dataHolder& data1, dataHolder& data2, int copyInfo)
{

    cout << "\nConcatenate data..." << endl;

    assert(data1.points.size()==data1.infos.size());
    assert(data2.points.size()==data2.infos.size());


    // two versions of this function, one with copy info, and one where I give a specific color per dataset
    if(copyInfo){
        Incremental_Tree tree(boost::make_zip_iterator(boost::make_tuple( data1.points.begin(),data1.infos.begin() )),
                  boost::make_zip_iterator(boost::make_tuple( data1.points.end(), data1.infos.end() ) ));

        for(int i = 0; i < data2.points.size(); i++)
        {
            Incremental_neighbor_search search(tree, data2.points[i]);
            Incremental_neighbor_search::iterator it = search.begin();

            // index of nearest neighbor in AP
            int idx = boost::get<1>(it->first).global_idx;
            if(data1.has_color)
                data2.infos[i].color = data1.infos[idx].color;
            if(data1.has_normal)
                data2.infos[i].normal = data1.infos[idx].normal;

            data1.points.push_back(data2.points[i]);
            data1.infos.push_back(data2.infos[i]);
        }

    }
    else{
//        std::array<unsigned char, 3> red = {240,128,128};
//        std::array<unsigned char, 3> blue = {30,144,255};

        for(int i = 0; i < data1.infos.size(); i++)
        {
            data1.infos[i].color = CGAL::red();
        }
        for(int i = 0; i < data2.infos.size(); i++)
        {
            data2.infos[i].color = CGAL::blue();
            data1.points.push_back(data2.points[i]);
            data1.infos.push_back(data2.infos[i]);
        }
    }
}

// removed this from helper.cpp, because it was the only function there which needed the IO library
int toXTensor(dataHolder& data){

    data.xpoints.clear();
    data.xnormals.clear();
    data.xgtnormals.clear();
    data.xsensor_positions.clear();
    for(int i = 0; i < data.points.size(); i++){
        data.xpoints.push_back(data.points[i].x());
        data.xpoints.push_back(data.points[i].y());
        data.xpoints.push_back(data.points[i].z());
        if(data.has_normal){
            data.xnormals.push_back(data.infos[i].normal.x());
            data.xnormals.push_back(data.infos[i].normal.y());
            data.xnormals.push_back(data.infos[i].normal.z());
        }
        if(data.has_gt_normal){
            data.xgtnormals.push_back(data.infos[i].gt_normal.x());
            data.xgtnormals.push_back(data.infos[i].gt_normal.y());
            data.xgtnormals.push_back(data.infos[i].gt_normal.z());
        }
        if(data.has_sensor){
            data.xsensor_positions.push_back(data.infos[i].sensor_positions[0].x());
            data.xsensor_positions.push_back(data.infos[i].sensor_positions[0].y());
            data.xsensor_positions.push_back(data.infos[i].sensor_positions[0].z());
        }
    }

//    for(int i = 0; i < data.gt_points.size(); i++){

//        data.xgt_points.push_back(data.gt_points[i].x());
//        data.xgt_points.push_back(data.gt_points[i].y());
//        data.xgt_points.push_back(data.gt_points[i].z());
//    }

    return 0;
}



//int importTransformationMatrix(dirHolder dir, dataHolder& data){

//    cout << "\nImport allignment matrix from " << dir.path+dir.transformation_file+".txt" << endl;

//    Eigen::MatrixXd  transformation_matrix(4, 4);


//    ifstream instream(dir.path+dir.transformation_file+".txt");
//    if(!instream){
//        cout << "\nERROR: File does not exist or is empty!" << endl;
//        return 1;
//    }

//    string line;
//    vector<string> lineAsVector;

//    for(int i = 0; i < 4; i++){
//        getline(instream, line);
//        splitString(line, lineAsVector, ' ');
//        transformation_matrix(i,0) = atof(lineAsVector[0].c_str());
//        transformation_matrix(i,1) = atof(lineAsVector[1].c_str());
//        transformation_matrix(i,2) = atof(lineAsVector[2].c_str());
//        transformation_matrix(i,3) = atof(lineAsVector[3].c_str());
//    }
//    data.transformation_matrix = transformation_matrix;
//    return 0;
//}



int importPLYPoints(const dirHolder dir, dataHolder& data)
{
    // this reads a file with sensor data, saved as sx, sy, sz
    // optional, see readColmapRecon() function for reading a file where sensor is saved as nx, ny, nz

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nRead PLY points..." << endl;
    cout << "\t-from " << dir.read_file << endl;

    std::string ifn = dir.path + dir.read_file + ".ply";

    // read Binary PLY with sensor
    Mesh_ply aMesh;
    Import_PLY(ifn.c_str(), &aMesh);

    assert(aMesh.mVertices.size() > 0);

    data.has_sensor = !(aMesh.mvCapture.size() == 0);
    data.has_normal = !(aMesh.mNormals.size() == 0);
    data.has_color = !(aMesh.mvColors.size() == 0);

    if(!data.has_sensor && !data.has_normal){
        cout << "\nWARNING: NEITHER NORMALS NOR SENSOR INFORMATION FOUND IN THE INPUT FILE!\n" << endl;
    }
    for(int i = 0; i < aMesh.mVertices.size(); i++){
        // parse points
        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
        data.points.push_back(pt);
        vertex_info vinf;

        // parse normals
        if(data.has_normal){
            Vector normal(aMesh.mNormals[i].x, aMesh.mNormals[i].y, aMesh.mNormals[i].z);
            vinf.normal = normal;
        }
        // sensor
        if(data.has_sensor){ // if there is a sensor location saved under scalar_x0
            Point sensor_position(aMesh.mvCapture[i].x, aMesh.mvCapture[i].y, aMesh.mvCapture[i].z);
            vinf.sensor_positions.push_back(sensor_position);
            Vector vec(vinf.sensor_positions[0].x() - pt.x(), vinf.sensor_positions[0].y() - pt.y(), vinf.sensor_positions[0].z() - pt.z());
            vinf.sensor_vec = vec;
//            double dist = 10000;
//            int new_sensor = 1;
//            while(dist != 0){
//                Point a_sensor = data.sensor_map.find(i).second();
//                dist=CGAL::squared_distance(a_sensor,sensor_position);
//            }

//            data.sensor_map

        }
        else if(data.has_normal){ // take the sensor location / normal that is saved under nx
            vinf.sensor_positions.push_back(Point(vinf.normal.x()+pt.x(), vinf.normal.y()+pt.y(), vinf.normal.z()+pt.z()));
            vinf.sensor_vec = vinf.normal;
        }

        // color
        unsigned char r,g,b;
        if(data.has_color){
            r = aMesh.mvColors[i].x * 255;
            g = aMesh.mvColors[i].y * 255;
            b = aMesh.mvColors[i].z * 255;
        }
        CGAL::Color col(r,g,b);
        vinf.color = col;
        // index, e.g. used in the concatenate data function
        vinf.global_idx = i;
        // save vertex info
        data.infos.push_back(vinf);
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-" << data.points.size() << " points read" << endl;
    if(data.has_color)
        cout << "\t-with colors" << endl;
    if(data.has_sensor)
        cout << "\t-with sensor" << endl;
    if(data.has_normal)
        cout << "\t-with normal" << endl;

    cout << "\t-done in " << duration.count() << "s" << endl;

    return 0;
}


//int importPLYMeshWithSensorStructure(dirHolder dir, dataHolder& data){

//    // TODO: use tinyply for importing ply meshes, and get rid of the custom version of rPLY which is currently used

//    string ifn=dir.path+dir.read_file+".ply";

//    std::cout << "\nRead mesh PLY file with (sensor) topology " << dir.read_file << std::endl;

//    auto start = chrono::high_resolution_clock::now();

//    // read Binary PLY with sensor
//    Mesh_ply aMesh;
//    Import_PLY(ifn.c_str(), &aMesh);

//    assert(aMesh.mVertices.size() > 0);

//    data.has_sensor = !(aMesh.mvCapture.size() == 0);
//    data.has_normal = !(aMesh.mNormals.size() == 0);
//    data.has_color = !(aMesh.mvColors.size() == 0);

//    if(!data.has_sensor && !data.has_normal){
//        cout << "\nWARNING: NEITHER NORMALS NOR SENSOR INFORMATION FOUND IN THE INPUT FILE!\n" << endl;
//    }
//    for(int i = 0; i < aMesh.mVertices.size(); i++){
//        // save points
//        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
//        data.points.push_back(pt);
//        // save infos
//        vertex_info vec_inf;
//        // sensor
//        if(aMesh.mvCapture.size()>0){
//            Vector vec(aMesh.mvCapture[i].x - pt.x(), aMesh.mvCapture[i].y - pt.y(), aMesh.mvCapture[i].z - pt.z());
//            vec_inf.sensor_vec = vec;
//            vec_inf.sensor_positions.push_back(Point(aMesh.mvCapture[i].x, aMesh.mvCapture[i].y, aMesh.mvCapture[i].z));
//        }
//        // color
//        unsigned char r,g,b;
//        if(aMesh.mvColors.size() > 0){
//            r = aMesh.mvColors[i].x * 255;
//            g = aMesh.mvColors[i].y * 255;
//            b = aMesh.mvColors[i].z * 255;
//            CGAL::Color col(r,g,b);
//            vec_inf.color = col;
//        }
//        // save vertex info
//        data.infos.push_back(vec_inf);
//    }

//    // save incident sensor triangles in each 3DT vertex
//    for(int i = 0; i < aMesh.mIndices.size()/3; i++){
//        array<size_t,3> poly;
//        size_t id0 = aMesh.mIndices[(i*3)+0];
//        size_t id1 = aMesh.mIndices[(i*3)+1];
//        size_t id2 = aMesh.mIndices[(i*3)+2];
//        poly[0]=id0;
//        poly[1]=id1;
//        poly[2]=id2;
//        data.facets.push_back(poly);

//        // save incident sensor triangles in each 3DT vertex
//        data.infos[id0].sensor_tet.push_back(i);
//        data.infos[id1].sensor_tet.push_back(i);
//        data.infos[id2].sensor_tet.push_back(i);
//        // basically just integrating the loop below in here.

//    }

//    auto stop = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
//    cout << "\t-points: " << data.points.size() << endl;
//    cout << "\t-facets: " << data.facets.size() << endl;
//    cout << "\t-in " << duration.count() << "s" << endl;
//}

int importPLYMesh(const dirHolder& dir, SurfaceMesh& import_mesh){
    cout << "\nRead file " << dir.read_file+".ply" << endl;
//    cout << "\nRead file " << dir.path+dir.read_file+".off" << endl;
    ifstream in(dir.path+dir.read_file+".ply");
    read_ply(in,import_mesh);
    if(!(import_mesh.number_of_vertices() > 0)){
        cerr << "File is empty!\n" << endl;
        return 1;
    }
    return 0;
}

int importMesh(const dirHolder& dir, dataHolder& data){

    // read file depending on file type
    if(dir.read_file_type == "off")
        return(importOFFMesh(dir, data.smesh));
    else if(dir.read_file_type == "ply")
        return(importPLYMesh(dir, data.smesh));
    else{
        cerr << "File type not specified or not support. Supported are .off and .ply files" << endl;
        return 1;
    }

}


int importOFFMesh(const dirHolder& dir, SurfaceMesh& import_mesh){

    cout << "Read file " << dir.read_file+".off" << endl;
    ifstream in(dir.path+dir.read_file+".off");
    in >> import_mesh;
    if(!(import_mesh.number_of_vertices() > 0)){
        cout << "File is empty!\n" << endl;
        return 1;
    }

    return 0;

}
int importOFFMesh(const dirHolder& dir, Polyhedron& import_poly){


    cout << "\nRead file " << dir.read_file+".off" << endl;
//    cout << "\nRead file " << dir.path+dir.read_file+".off" << endl;
    ifstream in(dir.path+dir.read_file+".off");
    in >> import_poly;
    if(!(import_poly.size_of_vertices() > 0)){
        cout << "File is empty!\n" << endl;
        return 1;
    }

    return 0;

}
int importOFFMesh(const string path, Polyhedron& import_poly){
    cout << "\nRead file " << path << endl;
//    cout << "\nRead file " << dir.path+dir.read_file+".off" << endl;
    ifstream in(path);
    in >> import_poly;
    if(!(import_poly.size_of_vertices() > 0)){
        cout << "File is empty!\n" << endl;
        return 1;
    }

    return 0;
};



//#include "cnpy.h"
#include <xtensor-io/xnpz.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>
int importNPZPoints(const dirHolder& dir, dataHolder& data){

    cout << "\nRead NPZ points..." << endl;
    cout << "\t-from " << dir.read_file << endl;

    /////// with xtensor
    auto a = xt::load_npz(dir.path+dir.read_file+".npz");

    if(a.find("points") == a.end()){
        cout << "\nERROR: No points array found in NPZ file!" << endl;
        return 1;
    }

    xt::xarray<double> normals;
    if(a.find("normals") != a.end()){
        normals = a["normals"].cast<double>();
        data.has_normal = true;
    }
    xt::xarray<double> sensor_pos;
    if(a.find("sensor_position") != a.end()){
        sensor_pos = a["sensor_position"].cast<double>();
        data.has_sensor = true;
    }
    auto points = a["points"].cast<double>();
    if(data.has_sensor)
        assert(points.size()==sensor_pos.size());
    if(data.has_normal)
        assert(points.size()==normals.size());

    data.points.clear();
    data.infos.clear();
    Point p,s;
    for(int i = 0; i < points.shape()[0]; i++){
        p = Point(points(i,0),points(i,1),points(i,2));
        data.points.push_back(p);
        s = Point(sensor_pos(i,0),sensor_pos(i,1),sensor_pos(i,2));
        vertex_info vinf;
        if(data.has_normal)
            vinf.normal = Vector(normals(i,0),normals(i,1),normals(i,2));
        if(data.has_sensor){
            vinf.sensor_positions.push_back(s);
            vinf.sensor_vec = s-p;
        }
        data.infos.push_back(vinf);
    }

    cout << "\t-" << data.points.size() << " points read" << endl;
    if(data.has_color)
        cout << "\t-with colors" << endl;
    if(data.has_sensor)
        cout << "\t-with sensor" << endl;
    if(data.has_normal)
        cout << "\t-with normal" << endl;

    return 0;



}



#ifdef OpenMVS
// fix with PI from here: https://github.com/cdcseacave/openMVS/issues/643
// otherwise nameclash with some CGAL headers
#pragma push_macro("PI")
#undef PI
#include "MVS.h"
#pragma pop_macro("PI")
using namespace MVS;

int importOMVSScene(const dirHolder dir, dataHolder& data){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nRead dense point cloud..." << endl;
    cout << "\t-from " << dir.path+"openMVS/densify_file.mvs" << endl;


    if ( !boost::filesystem::exists(dir.path+"openMVS/densify_file.mvs")){
        cout << "\nERROR: densify_file.mvs does not exist in the given path!" << std::endl;
        return 1;
    }

    Scene scene(0);
    scene.Load(dir.path+"openMVS/densify_file.mvs");


    // make a sensor map (for debbuging)
    FOREACH(k, scene.images){
        Image& imageData = scene.images[k];
        // skip invalid, uncalibrated or discarded images
        if (!imageData.IsValid())
            continue;
        // reset image resolution
//        imageData.ReloadImage(0, false);
        if(!imageData.ReloadImage(0, false)){
            #ifdef RECMESH_USE_OPENMP
            bAbort = true;
            #pragma omp flush (bAbort)
            continue;
//            #else
//            return EXIT_FAILURE;
            #endif
        }
        imageData.UpdateCamera(scene.platforms);
        // select neighbor views
        if (imageData.neighbors.IsEmpty()) {
            IndexArr points;
            scene.SelectNeighborViews(k, points);
        }
        Camera& camera = imageData.camera;
        data.sensor_map[k] = Point(camera.C.x,camera.C.y,camera.C.z);
    };


    // fetch points
    FOREACH(i, scene.pointcloud.points){
        // first check if the point has valid images attached
        const PointCloud::ViewArr& views = scene.pointcloud.pointViews[i];
        if(views.IsEmpty())
            continue;

        // read points
        const PointCloud::Point& P(scene.pointcloud.points[i]);
        Point current_point(P.x, P.y, P.z);
        data.points.push_back(current_point);

        vertex_info vinf;
        // read normal
        if(scene.pointcloud.normals.size()>0){
            const PointCloud::Normal& N(scene.pointcloud.normals[i]);
            Vector current_normal(N.x, N.y, N.z);
            vinf.normal = current_normal;
            data.has_normal = 1;
        }


        // read color
        const PointCloud::Color& C(scene.pointcloud.colors[i]);
        vinf.color = {C.r,C.g,C.b};

        // read views
        FOREACH(j, views) {
            const PointCloud::View viewID(views[j]);
            const Image& imageData = scene.images[viewID];
            if(!imageData.IsValid())
                continue;
            const Camera& camera = imageData.camera;
            vinf.sensor_positions.push_back(Point(camera.C.x, camera.C.y, camera.C.z));
        }
        // if no sensor is associated with point, delete point and go to next one (without adding the point info of course)
        if(vinf.sensor_positions.size() < 1){
            data.points.pop_back();
            continue;
        }
        vinf.sensor_vec = Vector(vinf.sensor_positions[0]-current_point);
        vinf.global_idx=i;
        data.infos.push_back(vinf);
    }

    // set has color and has normal
    data.has_color = 1;
    data.has_sensor = 1;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-" << data.points.size() << " points read" << endl;
    cout << "\t-in " << duration.count() << "s" << endl;

    return 0;
}
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// OUTPUT ////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int exportNPZ(const dirHolder& dir, dataHolder& data){

    if(data.xpoints.empty()){
        cerr << "\nERROR: xpoints is empty. Run toXTensor on dataHolder before exportNPZ." << endl;
        return 1;
    }
    else{
        cout << "\nExport points..." << endl;
        cout << "\t-to " << dir.write_file+dir.suffix+".npz" << endl;
        cout << "\t-" << data.xpoints.size()/3 << " points" << endl;
    }

    std::vector<std::size_t> shape = { data.points.size(), 3 };
    auto points = xt::adapt(data.xpoints, shape);
    xt::dump_npz(dir.path+dir.write_file+dir.suffix+".npz","points",points,true,false);

    if(data.has_normal){
        auto normals = xt::adapt(data.xnormals, shape);
        xt::dump_npz(dir.path+dir.write_file+dir.suffix+".npz","normals",normals,true,true);
    }
    if(data.has_gt_normal){
        auto gt_normals = xt::adapt(data.xgtnormals, shape);
        xt::dump_npz(dir.path+dir.write_file+dir.suffix+".npz","gt_normals",gt_normals,true,true);
    }
    if(data.has_sensor){
        auto sensor_position = xt::adapt(data.xsensor_positions, shape);
        xt::dump_npz(dir.path+dir.write_file+dir.suffix+".npz","sensor_position",sensor_position,true,true);
    }

    return 0;
}

//int export3DT(const dirHolder dir, dataHolder data){

//    boost::filesystem::path p(dir.write_file);
//    string outfile = p.stem().string();

//    cout << "\nExport 3DT..." << endl;
//    cout << "\t-to " << "dgnn/"+outfile+"_3dt.npz" << endl;

//    data.xverts.clear();
//    data.xfacets.clear();
//    data.xnfacets.clear();
//    data.xtets.clear();

//    for(Vertex_handle vh : data.Dt.finite_vertex_handles()){
//        data.xverts.push_back(vh->point().x());
//        data.xverts.push_back(vh->point().y());
//        data.xverts.push_back(vh->point().z());
//    }

//    Delaunay::Finite_facets_iterator fft;
//    int vidx;
//    Cell_handle c, mc;
//    for(fft = data.Dt.finite_facets_begin(); fft != data.Dt.finite_facets_end(); fft++){

//        c = fft->first; // cell
//        vidx = fft->second;

//        data.xnfacets.push_back(c->info().finite_idx);
//        mc = c->neighbor(vidx);
//        data.xnfacets.push_back(mc->info().finite_idx); // infinite cells will automatically have a finite_idx = -1

//        for(int j = vidx + 1; j <= vidx + 3; j++){
//                // so c->vertex() gives me the global vertex handle from the Dt
//                data.xfacets.push_back(c->vertex(j%4)->info().finite_idx);
//        }
//    }

//    for(Cell_handle ch : data.Dt.finite_cell_handles()){
//        data.xtets.push_back(ch->vertex(0)->info().finite_idx);
//        data.xtets.push_back(ch->vertex(1)->info().finite_idx);
//        data.xtets.push_back(ch->vertex(2)->info().finite_idx);
//        data.xtets.push_back(ch->vertex(3)->info().finite_idx);
//    }

//    cout << "\t-" << data.xverts.size()/3 << " vertices" << endl;
//    cout << "\t-" << data.xfacets.size()/3 << " facets" << endl;
//    cout << "\t-" << data.xtets.size()/4 << " tetrahedra" << endl;



//    std::vector<std::size_t> vshape = { data.Dt.number_of_vertices(), 3 };
//    auto verts = xt::adapt(data.xverts, vshape);
//    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","vertices",verts,true,false);

//    std::vector<std::size_t> fshape = { data.Dt.number_of_finite_facets(), 3 };
//    auto facets = xt::adapt(data.xfacets, fshape);
//    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","facets",facets,true,true);

//    std::vector<std::size_t> fnshape = { data.Dt.number_of_finite_facets(), 2 };
//    auto nfacets = xt::adapt(data.xnfacets, fnshape);
//    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","nfacets",nfacets,true,true);

//    std::vector<std::size_t> tshape = { data.Dt.number_of_finite_cells(), 4 };
//    auto tets = xt::adapt(data.xtets, tshape);
//    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","tetrahedra",tets,true,true);

//    return 0;

//}

int export3DT(const dirHolder dir, dataHolder& data){

    boost::filesystem::path p(dir.write_file);
    string outfile = p.stem().string();

    cout << "\nExport 3DT..." << endl;
    cout << "\t-to " << "dgnn/"+outfile+"_3dt.npz" << endl;

    data.xverts.clear();
    data.xvertsi.clear();
    data.xfacets.clear();
    data.xfacetsi.clear();
    data.xnfacets.clear();
    data.xtets.clear();
    data.xtetsi.clear();

    for(Vertex_handle vh : data.Dt.all_vertex_handles()){
//        cout << "vertex " << vh->info().global_idx << endl;
        data.xverts.push_back(vh->point().x());
        data.xverts.push_back(vh->point().y());
        data.xverts.push_back(vh->point().z());
        data.xvertsi.push_back(data.Dt.is_infinite(vh));
    }
//    auto shape = data.xvertsi.shape();
//    cout << "xvertsi shape " << xt::adapt(shape) << endl;

    Delaunay::All_facets_iterator fft;
    int vidx;
    Cell_handle c, mc;
    for(fft = data.Dt.all_facets_begin(); fft != data.Dt.all_facets_end(); fft++){

        c = fft->first; // cell
        vidx = fft->second;

        data.xnfacets.push_back(c->info().global_idx);
        mc = c->neighbor(vidx);
        data.xnfacets.push_back(mc->info().global_idx);

        for(int j = vidx + 1; j <= vidx + 3; j++){
                // so c->vertex() gives me the global vertex handle from the Dt
                data.xfacets.push_back(c->vertex(j%4)->info().global_idx);
        }
        data.xfacetsi.push_back(data.Dt.is_infinite(c));
    }

    for(Cell_handle ch : data.Dt.all_cell_handles()){
//        cout << "tets " << ch->info().global_idx << endl;
        data.xtets.push_back(ch->vertex(0)->info().global_idx);
        data.xtets.push_back(ch->vertex(1)->info().global_idx);
        data.xtets.push_back(ch->vertex(2)->info().global_idx);
        data.xtets.push_back(ch->vertex(3)->info().global_idx);
        data.xtetsi.push_back(data.Dt.is_infinite(ch));
    }

    cout << "\t-" << data.xverts.size()/3 << " vertices" << endl;
    cout << "\t-" << data.xfacets.size()/3 << " facets" << endl;
    cout << "\t-" << data.xtets.size()/4 << " tetrahedra" << endl;


    vector<size_t> vshape = { data.Dt.number_of_vertices()+1, 3 }; // number of vertices gives only finite number
    auto verts = xt::adapt(data.xverts, vshape);
    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","vertices",verts,true,false);
    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","inf_vertices",xt::adapt(data.xvertsi),true,true);


    vector<size_t> fshape = { data.Dt.number_of_facets(), 3 };
    auto facets = xt::adapt(data.xfacets, fshape);
    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","facets",facets,true,true);
    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","inf_facets",xt::adapt(data.xfacetsi),true,true);

    vector<size_t> fnshape = { data.Dt.number_of_facets(), 2 };
    auto nfacets = xt::adapt(data.xnfacets, fnshape);
    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","nfacets",nfacets,true,true);

    vector<size_t> tshape = { data.Dt.number_of_cells(), 4 };
    auto tets = xt::adapt(data.xtets, tshape);
    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","tetrahedra",tets,true,true);
    xt::dump_npz(dir.path+"dgnn/"+outfile+"_3dt.npz","inf_tetrahedra",xt::adapt(data.xtetsi),true,true);

    return 0;

}


///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// OFF ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void exportOFF(const dirHolder& dir, Polyhedron& poly)
{
    cout << "\nExport mesh..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".off" << endl;
    ofstream out(dir.path+dir.write_file+dir.suffix+".off");
    out << std::setprecision(15);
    out << poly;
    out.close();
}
void exportOFF(const dirHolder& dir, SurfaceMesh& out_mesh){
    cout << "\nExport mesh..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".off" << endl;
    ofstream out(dir.path+dir.write_file+dir.suffix+".off");
//    CGAL::set_binary_mode(out);
    out << std::setprecision(15);
    out << out_mesh;
    out.close();
}
void exportOFF(SurfaceMeshExact& out_mesh, std::string path){
    path = path + ".off";
    ofstream out(path);
    out << std::setprecision(15);
    out << out_mesh;
    out.close();
}
void exportOFF(Polyhedron& out_mesh, std::string path)
{
    path = path + ".off";
    ofstream out(path);
    out << std::setprecision(15);
    out << out_mesh;
    out.close();
    std::cout << "Exported to " << path << std::endl;
}
void exportOFF(Polyhedron_Exact& out_mesh, std::string path)
{
    path = path + ".off";
    ofstream out(path);
    out << std::setprecision(15);
    out << out_mesh;
    out.close();
    std::cout << "Exported to " << path << std::endl;
}
void exportOFF(Tetrahedron& in_tet, std::string path)
{
    Polyhedron out_poly;
    out_poly.make_tetrahedron(in_tet.vertex(0), in_tet.vertex(1), in_tet.vertex(2), in_tet.vertex(3));
    path = path + ".off";
    ofstream out(path);
    out << std::setprecision(15);
    out << out_poly;
    out.close();
}
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// PLY ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void printPLYHeader(fstream& fo, exportOptions& eo,
                    int nv, int nf=0,
                    int precision=8){


    if(eo.normals && eo.sensor_vec)
        std::cout << "\nprintPLYHeader: NORMALS AND SENSOR VECTORS BOTH WRITTEN TO PLY FILE WITH NX,NY,NZ\n" << std::endl;
    if(nf == 0 && eo.facetColor)
        std::cout << "\nprintPLYHeader: PROVIDE NUMBER OF FACETS IF YOU WANT TO EXPORT FACETS\n" << std::endl;

    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    if(eo.normals){
        fo << "property float nx" << std::endl;
        fo << "property float ny" << std::endl;
        fo << "property float nz" << std::endl;
    }
    if(eo.color){
        fo << "property uchar red" << std::endl;
        fo << "property uchar green" << std::endl;
        fo << "property uchar blue" << std::endl;
    }
    if(eo.sensor_vec){
        fo << "property float nx" << std::endl;
        fo << "property float ny" << std::endl;
        fo << "property float nz" << std::endl;
    }
//    else if(eo.sensor_vec && eo.normals){
//        fo << "property float sx" << std::endl;
//        fo << "property float sy" << std::endl;
//        fo << "property float sz" << std::endl;
//    }
    if(eo.sensor_position){
        fo << "property float sx" << std::endl;
        fo << "property float sy" << std::endl;
        fo << "property float sz" << std::endl;
    }
    if(eo.cam_index)
        fo << "property int camera_index" << std::endl;
    if(eo.score)
        fo << "property float score" << std::endl;
    if(nf > 0){
        fo << "element face " << nf << std::endl;
        fo << "property list uchar int vertex_indices" << std::endl;
        if(eo.facetColor){
            fo << "property uchar red" << std::endl;
            fo << "property uchar green" << std::endl;
            fo << "property uchar blue" << std::endl;
        }
    }
    fo << "end_header" << std::endl;
    fo << setprecision(precision);
};
int exportPLY(const dirHolder& dir, Point_set& points){

    cout << "\nExport points..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".ply" << endl;
    cout << "\t-" << points.size() << " points"  << endl;

    ofstream out(dir.path+dir.write_file+dir.suffix+".ply");

    return write_ply_point_set(out,points);

}
int exportPLY(const dirHolder& dir, vector<Point>& points){

    // make Point_set_3 from vector of points
    Point_set ps;
    ps.reserve(points.size());
    for(const auto p : points)
        ps.insert(p);

    cout << "\nExport point set..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".ply" << endl;
    cout << "\t-" << points.size() << " points"  << endl;

    ofstream out(dir.path+dir.write_file+dir.suffix+".ply");

    return write_ply_point_set(out,ps);

}
void exportPLY(const dirHolder& dir, vector<Point>& points, vector<vertex_info>& infos, exportOptions& eo){


    cout << "\nExport points..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".ply" << endl;

    if(eo.normals || eo.color || eo.sensor_vec || eo.sensor_position)
        assert(points.size()==infos.size());

    // get number of vertices and triangles of the triangulation
    int nv = points.size();

    fstream fo;
    fo.open(dir.path+dir.write_file+dir.suffix+".ply", std::fstream::out);

    printPLYHeader(fo, eo, nv);

    for(int i = 0; i < points.size(); i++){
        // print data to file
        assert(points.size() > 0);
        fo << points[i] << " ";                           // coordinates
        if(eo.normals){
//            assert(infos[i].normal.x());
            fo << infos[i].normal << " ";
        }
        if(eo.color){
//            assert(infos[i].color[0]);
//            fo << "255 0 0 ";
            fo << int(infos[i].color.r()) << " " << int(infos[i].color.g()) << " " << int(infos[i].color.b()) << " ";
        }
        // output the DIRECTION to the sensor (not the SENSOR LOCATION!!) if you want to visualize it in e.g. Meshlab
        if(eo.sensor_vec){
//            assert(infos[i].sensor_vec.x());
            fo << infos[i].sensor_vec << " ";
        }
        if(eo.sensor_position){
//            assert(infos[i].sensor_positions[0].x());
            fo << infos[i].sensor_positions[0] << " ";
        }
        fo << endl;
    }

    fo.close();

    cout << "\t-" << nv << " points"  << endl;
}
int exportPLY(const dirHolder& dir, SurfaceMesh& out_mesh){

    cout << "\nExport mesh..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".ply" << endl;
    ofstream out(dir.path+dir.write_file+dir.suffix+".ply");

    return write_ply(out,out_mesh);

}
void exportPLY_bin(const dirHolder& dir, SurfaceMesh& out_mesh){

    // this code is from here: https://piero.dev/2017/04/how-to-write-a-cgal-polyhedron-to-ply-or-any-other-format/

    Polyhedron poly;
    CGAL::copy_face_graph(out_mesh, poly);

    cout << "\nExport mesh..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".bin.ply" << endl;

    ofstream os(dir.path+dir.write_file+dir.suffix+".bin.ply", ios::out | ios::binary);

    string tests = "ply";

    os.write(tests.c_str(), tests.size());
    os.close();
}
typedef typename Polyhedron::Vertex_const_iterator VCI;
typedef typename Polyhedron::Facet_const_iterator FCI;
typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;
void exportPLY_manual(const dirHolder& dir, SurfaceMesh& out_mesh){

    // this code is from here: https://piero.dev/2017/04/how-to-write-a-cgal-polyhedron-to-ply-or-any-other-format/

    Polyhedron poly;
    CGAL::copy_face_graph(out_mesh, poly);

    cout << "\nExport mesh..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".ply" << endl;

//    std::filebuf fb;
//    fb.open(dir.path+dir.write_file+dir.suffix+".ply", std::ios::out | std::ios::binary);
//    std::ostream os(&fb);

    ofstream os(dir.path+dir.write_file+dir.suffix+".ply", ios::out | ios::binary);

    os << "ply\n"
       << "format ascii 1.0\n"
       << "element vertex " << poly.size_of_vertices() << "\n"
       << "property float x\n"
       << "property float y\n"
       << "property float z\n"
       << "element face " << poly.size_of_facets() << "\n"
       << "property list uchar int vertex_indices\n"
       << "end_header\n";

    for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++){
        os << it->point().x() << " " << it->point().y() << " " << it->point().z() << std::endl;
    }

    typedef CGAL::Inverse_index<VCI> Index;
    Index index(poly.vertices_begin(), poly.vertices_end());

    for( FCI fi = poly.facets_begin(); fi != poly.facets_end(); ++fi) {
        HFCC hc = fi->facet_begin();
        HFCC hc_end = hc;

        os << circulator_size(hc) << " ";
        do {
            os << index[VCI(hc->vertex())] << " ";
            ++hc;
        } while( hc != hc_end);

        os << "\n";
    }

//    fb.close();
    os.close();
}
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// CUSTOM ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////





void exportIsoValues(const dirHolder& dir, Delaunay& Dt, vector<double>& values){

    std::ofstream fo(dir.path+dir.write_file+"_isovalue.ply");
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << values.size() << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "end_header" << std::endl;
    int j = 0;
    for(auto vit = Dt.finite_vertices_begin(); vit != Dt.finite_vertices_end(); vit++){
//        fo << vit->point() << " 0 " << (int)(values[j++]*255) << " 0" << endl;
        if(values[j++]>0.5)
            fo << vit->point() << " 255 0 0" << endl;
        else
            fo << vit->point() << " 0 0 255" << endl;
    }
}




void exportCameraCenter(const dirHolder& dir, dataHolder& data){

    cout << "\nExport cameras locations to cameras.ply" << endl;
    // get number of cameras
    int nc = data.sensor_map.size();

    fstream fo;
    fo.open(dir.path+"cameras.ply", std::fstream::out);

    exportOptions eo;
    eo.cam_index = true;
    eo.color = true;
    printPLYHeader(fo, eo, nc);

    assert(data.sensor_map.size() > 0);
    for(int i = 0; i < nc; i++){
        // print data to file
        auto current = data.sensor_map.find(i);
        fo << current->second << " ";
        fo << "0 255 0 ";
        fo << current->first << " ";
        fo << endl;
    }
    fo.close();
    cout << "\t-" << nc << " cameras exported"  << endl;
}
void exportCellCenter(const dirHolder& dir, const Delaunay& Dt){

    Delaunay::size_type nc = Dt.number_of_finite_cells();

    fstream fo;
    fo.open(dir.path+dir.write_file+"_cellCenter.ply", fstream::out);
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << nc << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
//    fo << "property float inside_score" << std::endl;
//    fo << "property float outside_score" << std::endl;
//    fo << "property float inside" << std::endl;
//    fo << "property float outside" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(8);

    // calc min/max of scores for scaling
    Delaunay::Finite_cells_iterator cit;
//    std::vector<double> inside_scores;
//    std::vector<double> outside_scores;
//    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
//        inside_scores.push_back(cit->info().inside_score);
//        outside_scores.push_back(cit->info().outside_score);
//    }
//    double inside_min = *std::min_element(inside_scores.begin(), inside_scores.end());
//    double outside_min = *std::min_element(outside_scores.begin(), outside_scores.end());
//    double inside_max = *std::max_element(inside_scores.begin(), inside_scores.end());
//    double outside_max = *std::max_element(outside_scores.begin(), outside_scores.end());

    // print scores to file
    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        Point p1 = cit->vertex(0)->point();
        Point p2 = cit->vertex(1)->point();
        Point p3 = cit->vertex(2)->point();
        Point p4 = cit->vertex(3)->point();

        Point centroid = CGAL::centroid(p1,p2,p3,p4);


//        double inside_score = cit->info().inside_score;
//        double outside_score = cit->info().outside_score;
//        int green = 0;
//        int red  = int(255*(inside_score-inside_min)/(inside_max-inside_min));
//        int blue = int(255*(outside_score-outside_min)/(outside_max-outside_min));
//        if(inside_max == 0 || outside_max == 0)
//            red = 0;
//            blue = 0;
//        if(inside_score == outside_score){
//            green = 128;
//        }
//        fo << centroid << " " << red << " " << green << " " << blue << " "
//           << inside_score << " " << outside_score << " ";
////           << cit->info().inside_count << " " << cit->info().outside_count << " ";
//        if(cit->info().gc_label == 0)
//            fo << "255 0";
//        else
//            fo << "0 255";

        float os = cit->info().outside_score;
        float is = cit->info().inside_score;
        if(os > is)
            fo << centroid << " 0 0 " << (int)(os*255) << endl;
        else if(os < is)
            fo << centroid << " " << (int)(is*255) <<  " 0 0" << endl;
        else
            fo << centroid << " 0 255 0" << endl;



    }
    fo.close();
    std::cout << "\nExported colored cell center to " << dir.write_file+dir.suffix << std::endl;

}
void exportCellScore(const dirHolder& dir, const Delaunay& Dt){
    ofstream fo(dir.path+dir.write_file+"_cellScore.ply");
    Point p1, p2, p3, p4, centroid;
    Vector v;
    stringstream pstr, tstr;
    int vcount = 0;
    int fcount = 0;
    for(auto cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){

        if(Dt.is_infinite(cit))
            continue;

        stringstream cstr;
        float os = cit->info().outside_score;
        float is = cit->info().inside_score;
        if(os > is)
            cstr << "0 0 " << (int)(os*255) << endl;
        else if(os < is)
            cstr << (int)(is*255) <<  " 0 0" << endl;
        else
            cstr << "0 255 0" << endl;

        p1 = cit->vertex(0)->point();
        p2 = cit->vertex(1)->point();
        p3 = cit->vertex(2)->point();
        p4 = cit->vertex(3)->point();
        centroid = CGAL::centroid(p1,p2,p3,p4);

        for(int i = 0; i < 4; i++){

            Facet current_facet(cit, i);
            fcount++;

            tstr << "3 ";
            Triangle tri = Dt.triangle(current_facet);

            for(int j = 0; j < 3; j++){
                p1 = tri.vertex(j);
                v = centroid - p1;
                p1 = p1 + v*0.05;
                pstr << p1 << endl;
                tstr << vcount++ << " ";
            }
            tstr << cstr.str();
        }
    }
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << vcount << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "element face " << fcount << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "end_header" << std::endl;
    fo << pstr.str() << tstr.str();
    cout << "\nExported colored cell to " << dir.write_file+"_cellScore.ply" << endl;

}

void exportConvexHull(const dirHolder& dir, const Delaunay& Dt)
{


    // TODO: add cellScore output as points to the convexHull so one can visualize both in parallel, should work


    // this used to be exportColoredFacets, but exportConvexHull was a bit obsolete, because this also exports
    // the convex hull, just with colored facets.

    auto start = std::chrono::high_resolution_clock::now();

    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();

    fstream fo;
    fo.open(dir.path+dir.write_file+"_convexHull.ply", fstream::out);
    exportOptions eo;
    eo.normals = true; eo.color = true; eo.facetColor = true;
    printPLYHeader(fo, eo,
                   nv, nf);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates, color and normal of the vertex to the PLY file
    int index = 0;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        // reset the vertex index here, because I need to know the order of exactly this loop here
        // for the indexing of the facets in the PLY file
        vft->info().global_idx = index;
        // print data to file
        // coordinates
        fo  << vft->point() << " "
        // normal
            << vft->info().sensor_vec << " "
        // color
            << int(vft->info().color[0]) <<  " " << int(vft->info().color[1]) <<  " " << int(vft->info().color[2]) <<  std::endl;
        index++;
    }

    // Save the facets to the PLY file
    int vidx;
    Cell_handle c;
    Delaunay::Finite_facets_iterator fft;
    for(fft = Dt.finite_facets_begin() ; fft != Dt.finite_facets_end() ; fft++){

        // get vertex and cell index that describes the facet
        // facet fft is represented by std::pair(cell c, int vidx). vidx is the vertex opposite to the cell.
        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite.
        // see: https://doc.cgal.org/latest/Triangulation_3/index.html
        c = fft->first;         // cell
        vidx = fft->second;     // vertex index
        fo << 3 << ' ';
        // fix the orientation
        for(int j = vidx + 1; j <= vidx + 3; j++)
            fo << c->vertex(j%4)->info().global_idx << ' ';

//        float os = c->info().outside_score;
//        float is = c->info().inside_score;
//        fo << " " << (int)(is*255) << " 0 " << (int)(os*255);

//        float facet_score = c->info().facet_weights[vidx];
//        fo << " 0 0 " << (int)(facet_score*255);

        int i = c->info().gc_label;
        Cell_handle d = c->neighbor(i);
        int j = d->info().gc_label;

        if(i == 1 && j == 1)
            fo << "212 239 223";
        else
            fo << "82 190 128";
//        else(i == 0 && j == 0)
//            fo << "82 190 128";
//        else
//            fo << "255 0 0";

        fo << std::endl;
    }
    fo.close();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "\nExported convex hull with colored facets to " << dir.write_file+dir.suffix << " in " << duration.count() << "s" <<  std::endl;
}

void exportInterface(const dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& eo){

    auto start = std::chrono::high_resolution_clock::now();

    cout << "\nExport interface..." << endl;
    cout << "\t-to " <<  dir.write_file+".ply" << endl;

    data.remaining_points.clear();
    data.remaining_facets.clear();

    // get the remaining points and facets
    std::set<Vertex_handle> remaining_set;
    std::vector<std::vector<Vertex_handle>> remaining_facets;
    Delaunay::Finite_facets_iterator fft;
    int deletedFaceCount = 0;
    int clabel, mlabel, vidx;
    Cell_handle c,m;
    Facet current_facet;
    for(fft = data.Dt.finite_facets_begin(); fft != data.Dt.finite_facets_end(); fft++){


        // get vertex and cell index that describes the facet
        // facet fft is represented by std::pair(cell c, int vidx). vidx is the vertex opposite to the cell.
        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite.
        // see: https://doc.cgal.org/latest/Triangulation_3/index.html

        c = fft->first; // cell
        // this is needed to make sure the orientation check succeds
        // because it will fail for infinite cells
        if(data.Dt.is_infinite(c)){
            current_facet = data.Dt.mirror_facet(*fft);
            c = current_facet.first;
        }
        else
            current_facet = *fft;

        ////////// THIS HAS TO BE WRITTEN AFTER THE IF-ELSE ABOVE!!!!!!!!! ////////////////////
        vidx = current_facet.second;     // cell based vertex index
        ///////////////////////////////////////////////////////////////////////////////////////

        clabel = c->info().gc_label;
        m = data.Dt.mirror_facet(current_facet).first;
        mlabel = m->info().gc_label;

        // simply export cell center points here to see what is going on.

        // check if two neighbouring cells have the same label, and if so continue to next face
        if(clabel == mlabel){
            deletedFaceCount++;
            continue;
        }

        // if label of neighbouring cells is not the same...
        // if opposite vertex vidx is 2, we start at j = vidx + 1 = 3, 3%4 = 3
        // next iteration: j = 4, 4%4 = 0, next iteration: j = 5, 5%4 = 1;
        // so we exactely skip 2 - the opposite vertex.


        // this gets the triangle vertex handles. however the orientation of the triangle is arbitrary
        vector<Vertex_handle> tri;
//        for(int j = vidx + 1; j < vidx + 4; j++){        // should be changed to this
        // explaination here: https://stackoverflow.com/questions/7938311/cgal-help-getting-triangles-coordinates-from-delaunay-triangulation
        for(int j = vidx + 1; j <= vidx + 3; j++){
            // so c->vertex() gives me the global vertex handle from the Dt
            tri.push_back(c->vertex(j%4));
            remaining_set.insert(c->vertex(j%4));
        }

        // fix the orientation
        // check if the fourth point of the outside/inside tetrahedron is on
        // the positive/negative side, to get consistent orientation of the triangle
        EPICK::Plane_3 p(tri[0]->point(), tri[1]->point(), tri[2]->point());
        if(p.has_on_positive_side(c->vertex(vidx)->point()) && clabel == 0){
            Vertex_handle temp = tri[1];
            tri[1] = tri[2];
            tri[2] = temp;
        }
        else if(p.has_on_negative_side(c->vertex(vidx)->point()) && clabel == 1){
            Vertex_handle temp = tri[1];
            tri[1] = tri[2];
            tri[2] = temp;
        }
        remaining_facets.push_back(tri);

    } // end of finite facet iterator

    int bnv = data.Dt.number_of_vertices();
    int nv = remaining_set.size();

    int bnf = data.Dt.number_of_finite_facets();
    int nf = bnf - deletedFaceCount;

    ///////// file output ////////
    // create PLY output file for outputting the triangulation, with point coordinates, color, normals and triangle facets

    fstream fo;
    fo.open(dir.path+dir.write_file+dir.suffix+".ply", fstream::out);
    printPLYHeader(fo, eo, nv, nf);

    std::set<Vertex_handle>::iterator vit;
    vidx = 0;
    for(vit = remaining_set.begin(); vit != remaining_set.end(); vit++){
        // reset the vertex index here, because I need to know the order of exactly this loop here
        // for the indexing of the facets in the PLY file
        (*vit)->info().global_idx = vidx++;
        // also create a remaining_points vector
        data.remaining_points.push_back((*vit)->point());
        // print data to file
        // coordinates
        fo  << (*vit)->point() << " ";
        // normal
        if(eo.normals)
            fo << (*vit)->info().normal << " ";
        if(eo.color)
            fo << int((*vit)->info().color[0]) <<  " " << int((*vit)->info().color[1]) <<  " " << int((*vit)->info().color[2]) << " ";
        if(eo.sensor_vec)
            fo << (*vit)->info().sensor_vec << " ";
        if(eo.sensor_position)
            fo << (*vit)->info().sensor_positions[0] << " ";
        fo << std::endl;
    }

    std::vector<std::vector<Vertex_handle>>::iterator fit;
    for(fit = remaining_facets.begin(); fit != remaining_facets.end(); fit++){
        // start printed facet line with a 3
        fo << 3 << ' ';
        std::array<size_t,3> poly;
        for(int j = 0; j < 3; j++){
            // fill the "remaining polygon"-vector
            poly[j] = (*fit)[j]->info().global_idx;
            // print the indicies of each cell to the file
            fo << (*fit)[j]->info().global_idx << ' ';
        }
        data.remaining_facets.push_back(poly);
        // put a color on the face, so that in Meshlab I can activate the color by facet mode, to compare with the "colored facet file"
        if(eo.facetColor)
            fo << "0 0 200";
        fo << std::endl;
    }
    fo.close();

    if(options.optimization){
        std::cout << "\t-Before optimization (points, facets): (" << bnv << ", " << bnf << ")" << std::endl;
        std::cout << "\t-After optimization (points, facets): (" << nv << ", " << nf << ")" << std::endl;
    }
    else{
        std::cout << "\t-Before pruning (points, facets): (" << bnv << ", " << bnf << ")" << std::endl;
        std::cout << "\t-After pruning (points, facets): (" << nv << ", " << nf << ")" << std::endl;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "\t-done in " << duration.count() << "s" << std::endl;
}




