#include <exe/inputParser.h>
#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <IO/ethIO.h>
#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>

#ifdef OpenMVS
#include <IO/mvIO.h>
#endif

int main(int argc, char const *argv[]){


    cliParser ip("eth");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getETH())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----DELAUNAY-GRAPH-CUT-BASED SURFACE RECONSTRUCTION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    if(ip.dh.write_file.empty())
        ip.dh.write_file = ip.dh.read_file;

//    // load bounding box
//    cout << "\nRead bounding box from " << ip.dh.path+ip.dh.read_file+"_bb.ply" << endl;
//    std::ifstream in(ip.dh.path+ip.dh.read_file+"_bb.ply");
//    CGAL::read_ply(in,data.sm_obb);
//    if(data.sm_obb.is_empty()){
//        cout << "\nERROR: empty bounding box" << endl;
//        return 1;
//    }

    // make and export bounding box
    importLidarPoints(ip.dh,data);
    makeOrientedBoundingBox(data,ip.ro);
    std::ofstream out(ip.dh.path+ip.dh.write_file+"_obb.ply");
    CGAL::write_ply(out,data.sm_obb);

    data.points.clear();
    data.infos.clear();
    exportOptions eo;


    // load and clip ground truth
    if(importETH3D(ip.dh,data)){
        cout << "\nERROR: in importETH3D()" << endl;
        return 1;
    }
//    ip.dh.suffix = "_gt";
//    eo.normals = false;
//    eo.sensor_position = false;
//    eo.sensor_vec = true;
//    exportPLY(ip.dh,data.points,data.infos,eo);

//    clipWithBoundingBox(data.points,data.infos,data.sm_obb);
    ip.dh.suffix = "_gt_clipped";
    eo.normals = false;
    eo.sensor_position = false;
    eo.sensor_vec = true;
    exportPLY(ip.dh,data.points,data.infos,eo);



    runningOptions ro;
    ro.normal_method = ip.ro.normal_method;
    ro.normal_neighborhood = ip.ro.normal_neighborhood;
    ro.normal_orientation = ip.ro.normal_orientation;
    estimateNormals(data,ro);

    ip.dh.suffix = "_gt_normals";
    eo.normals = true;
    eo.sensor_position = false;
    eo.sensor_vec = false;
    exportPLY(ip.dh,data.points,data.infos,eo);


//    // load and clip dense cloud
//    if(loadOMVSScene(ip.dh, data))
//        return 1;
//    clipWithBoundingBox(data.points,data.infos,data.sm_obb);
//    ip.dh.suffix = "_dense_clipped";
//    exportPLY(ip.dh,data.points,data.infos,eo);

    // make closed volume by clipping bounding box and mesh


    // save the closed volume to data.gt_poly




    // clip dense cloud with bounding box


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----DELAUNAY-GRAPH-CUT-BASED SURFACE RECONSTRUCTION FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



