#ifdef OpenMVS

#include "IO/mvIO.h"
#include <base/cgal_typedefs.h>
// fix with PI from here: https://github.com/cdcseacave/openMVS/issues/643
// otherwise nameclash with some CGAL headers
#pragma push_macro("PI")
#undef PI
#include "MVS.h"
#pragma pop_macro("PI")
#include <boost/filesystem.hpp>

using namespace MVS;
using namespace std;

int loadOMVSScene(dirHolder dir, dataHolder& data){

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
