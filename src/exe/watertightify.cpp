#include <IO/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

#ifdef OpenMVS
#include <IO/mvIO.h>
#endif

#include <util/helper.h>
#include <util/geometricOperations.h>

#include <processing/meshProcessing.h>
#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>
#include <processing/evaluation.h>

#include <boost/filesystem.hpp>
using namespace boost::filesystem;



int main(int argc, char const *argv[]){

    cliParser ip("scan");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getScan())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----WATERTIGHTIFY MESH-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;

    // 1. load the watertight mesh
    if(importOFFMesh(ip.dh.path+ip.dh.read_file+".off", data.gt_poly))
        return 1;

    if(ip.ro.gt_type == "closed"){
        // 2. scan the mesh
        if(scanObjectClosed(data, ip.ro))
            return 1;
    }
    else{
        if(scanObjectOpen(data, ip.ro))
            return 1;
    }

    // 3. estimate normals for the scan
    data.has_normal = false;
    if(!ip.ro.normal_method.empty()){
        if(estimateNormals(data, ip.ro))
            return 1;
    }

    // 4. export the scan with sensor vectors and normals

    if(ip.dh.write_file.empty())
        ip.dh.write_file = ip.dh.read_file;

    exportOptions eo;
    eo.sensor_vec = false;
    eo.sensor_position = false;
    eo.normals = true;
    if(ip.eo.toply){
//        eo.sensor_vec = true;
//        eo.normals = false;
        exportPLY(ip.dh, data.points, data.infos, eo);
    }

    eo.normals = true;
    if(ip.eo.tonpz){
        toXTensor(data);
        exportNPZ(ip.dh,data);
    }



    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----WATERTIGHTIFY MESH FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



