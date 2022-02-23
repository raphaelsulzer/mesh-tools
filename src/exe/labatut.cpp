#include <exe/labatut.h>
#include <exe/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <IO/ttIO.h>
#include <IO/ethIO.h>
#ifdef COLMAP
#include <IO/colmapIO.h>
#endif
#include <util/helper.h>
#include <util/geometricOperations.h>

#include <processing/meshProcessing.h>
#include <processing/graphCut.h>
#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>
#include <processing/evaluation.h>
#include <processing/rayTracingFacet.h>


#ifdef Open3D
#include "open3d/Open3D.h"
#include "open3d/geometry/TetraMesh.h"
#endif

#include <boost/filesystem.hpp>
using namespace boost::filesystem;


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
///////////////////////// CONTROL FUNCTIONS /////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
int runLabatut(dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& exportO){

    // TOOD: add an assert for checking that there are no duplicate points in the input for the 3DT, because the ray tracing crashes
    // if there are. Maybe it would however be best to change my data.points+data.infos vector to a CGAL::Point_set_3
    // with property maps for each info.



    ///////////////////////////////
    ///////// FILE NAMING /////////
    ///////////////////////////////
    options.scoring="_"+options.scoring;
    dir.rw_string = double2string(options.area_reg_weight+options.angle_reg_weight+options.cc_reg_weight+options.sv_reg_weight);
    if(options.percentage_of_outliers > 0.0)
        dir.ol_string = double2string(options.percentage_of_outliers*100);
    if(dir.read_file.empty())
        dir.read_file = dir.write_file;
    if(dir.write_file.empty())
        dir.write_file = dir.read_file+dir.ol_string+options.scoring+dir.rw_string;
    else
        dir.write_file+=options.scoring+dir.rw_string;

    ///////////////////////////////
    ///////// IMPORT DATA /////////
    ///////////////////////////////
    // reconstruction input
    if(options.data_source == "ply"){
        importPLYPoints(dir, data);
    }
    else if(options.data_source == "npz"){
        importNPZ(dir, data);
    }
    #ifdef COLMAP
    else if(options.data_source == "colmap"){
        readColmapFiles(dir, data);
    }
    #endif
    #ifdef OpenMVS
    else if(options.data_source == "omvs"){
        if(importOMVSScene(dir, data))
            return 1;
    }
    #endif
    else{
        cout << "ERROR: not a valid reconstruction input" << endl;
        return 1;
    }


    ///////////////////////////////////
    ///////// PREPROCESS DATA /////////
    ///////////////////////////////////

//    if(options.scale > 0.0)
//        scalePointSet(dir, data, options.scale);

    if(options.scale > 0.0)
        standardizePointSet(dir, data, options.scale);

//    if(options.scoring == "_rt")
//        orderSensors(data);
//    else
//        cout << "\nSensors not ordered for ray-tracing" << endl;

    if(exportO.cameras){
        dir.suffix = "_cameras";
        exportCameraCenter(dir, data);
    }
    // export the scanned points
    if(exportO.scan){
        dir.suffix = "_scan";
        exportPLY(dir, data.points, data.infos, exportO);
    }


    //////////////////////////////////////
    /////// DELAUNAY TRIANGULATION ///////
    //////////////////////////////////////
    if(options.Dt_epsilon > 0.0)
        makeAdaptiveDelaunayWithInfo(data, options.Dt_epsilon);
    else
        makeDelaunayWithInfo(data);

    if(options.labatut_sigma == -1.0 && (options.scoring == "_rt" || options.scoring == "_clrt"))
        options.labatut_sigma = pcaKNN(data.Dt);
    else{
        cout << "\nd parameter of Labatu set to " << options.labatut_sigma << endl;
    }

    //////////////////////////////////////
    ///////////// RAY TRACING ////////////
    //////////////////////////////////////
    auto rc = processing::RayCaster(data, options);
    rc.run(1);



    ////////////////////////////
    /////// OPTIMIZATION ///////
    ////////////////////////////
    if(options.optimization == 1){
        options.gco_iterations = -1;
        if(options.score_type == "mine" || options.scoring == "_cl" || options.scoring == "_cs" || options.scoring == "_clcs" || options.scoring == "_lrtcs")
            graphCutTet(data, options);
        else
            graphCutFacet(data,options);
    }
    else{
        cout << "\nLabel cells without optimization by taking max score" << endl;
        cout << "\t-infinite cells will be labelled outside" << endl;
        for(auto cit = data.Dt.all_cells_begin(); cit != data.Dt.all_cells_end(); cit++){
            if(data.Dt.is_infinite(cit)){
                cit->info().gc_label = 1;
                continue;
            }
            // this means untraversed cells and 50/50 cells will be labelled as inside
            cit->info().gc_label = cit->info().outside_score > cit->info().inside_score ? 1 : 0;
        }
    }

    ////////////////////////////
    ////////// EXPORT //////////
    ////////////////////////////

    //// export interface
    if(exportO.interface)
        exportInterface(dir, data, options, exportO);

    #ifdef OpenMVS
    if(options.clean_mesh)
        omvsCleanMesh(dir,data);
    #endif

//    if(exportO.coloredFacets)
//        exportColoredFacets(dir, data.Dt, options.optimization);
    if(exportO.cellScore){
        exportCellCenter(dir, data.Dt);
        exportCellScore(dir, data.Dt);
    }
    if(exportO.convexHull){
        exportConvexHull(dir, data.Dt);
    }
    ////////// create surface mesh ////////////
    // set export options
    meshProcessingOptions mpOptions;
    mpOptions.try_to_close = options.try_to_close;
    mpOptions.try_to_make_manifold = options.try_to_make_manifold;
    mpOptions.number_of_components_to_keep = options.number_of_components_to_keep;
    mpOptions.factor_for_removing_large_faces = options.factor_for_removing_large_faces;
    if(exportO.mesh){
        dir.suffix = "_mesh";
        createSurfaceMesh(data, mpOptions);
        exportPLY(dir, data.smesh);
    }


    ////////// sample surface mesh //////////
    if(exportO.sampling){
        data.points.clear();
        data.infos.clear();
        if(data.smesh.is_empty()){
            createSurfaceMesh(data, mpOptions);
        }
        sampleMesh(data, exportO);
        dir.suffix = "_sampled";
        // turn off color and sensor_vec export, because mesh sampling obviously does not have that
        exportO.sensor_vec = false;
        exportO.sensor_position = false;
        exportO.color = false;
        exportO.normals = false;
        exportPLY(dir, data.points, data.infos, exportO);
    }

    if(options.evaluate_mesh){
        if(data.smesh.is_empty()){
            createSurfaceMesh(data, mpOptions);
        }
        if(options.ground_truth){
            double iou;
            if(calcIOU(data, exportO.sampling_method_option, iou))
                return 1;
            printMeshEvaluation(dir,data,iou);
        }
        else
            printMeshEvaluation(dir,data,-1.0);

    }

    return 0;
}



int main(int argc, char const *argv[]){


    cliParser ip("labatut");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getLabatut())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----LABATUT SURFACE RECONSTRUCTION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    if(runLabatut(ip.dh, data, ip.ro, ip.eo))
        return 1;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----LABATUT SURFACE RECONSTRUCTION FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



