#include <exe/sure.h>

#include <exe/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

#include <util/helper.h>
#include <util/vectorArithmetic.h>

#include <processing/meshProcessing.h>
#include <processing/edgeManifoldness.h>
#include <processing/optimizationTet.h>
#include <processing/optimizationFacet.h>

#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>
#include <processing/evaluation.h>
#include <processing/rayTracingFacet.h>


#include <learning/learning.h>
#include <learning/learningMath.h>
#include <learning/learningRayTracing.h>
#include <learning/learningRayTracingGroundTruth.h>
#include <learning/learningIO.h>

#ifdef Open3D
#include "open3d/Open3D.h"
#include "open3d/geometry/TetraMesh.h"
#endif

#include <boost/filesystem.hpp>
using namespace boost::filesystem;


int main(int argc, char const *argv[]){


    cliParser ip("sample");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getSample())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----SAMPLE SURFACE RECONSTRUCTION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;



    dataHolder data;
    meshProcessingOptions options;
    options.number_of_components_to_keep = ip.ro.number_of_components_to_keep;
    options.try_to_close = ip.ro.try_to_close;
    options.try_to_make_manifold = ip.ro.try_to_make_manifold;

    // read file depending on file type
    if(ip.dh.read_file_type == "off")
        importOff(ip.dh, data.smesh);
    else if(ip.dh.read_file_type == "ply"){
        importLidarMesh(ip.dh, data);
        if(data.facets.size() > 0){
            createSurfaceMesh(data, options);
        }
    }
    else{
        cerr << "Could not read file ending or file type not supported. Supported are .off and .ply files" << endl;
        return 1;
    }



    ////////// sample surface mesh //////////
    sampleMesh(data, ip.eo);
    ip.dh.suffix = "_sampled";
    // turn of color and sensor_vec export, because mesh sampling obviously does not have that
    exportOptions exportO;
    exportO.sensor_vec = false;
    exportO.sensor_position = false;
    exportO.color = false;
    exportO.normals = false;
    exportPLY(ip.dh, data.points, data.infos, exportO);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----SAMPLE SURFACE RECONSTRUCTION FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



