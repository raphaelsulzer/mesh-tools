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
#include <processing/edgeManifoldness.h>
#include <processing/graphCut.h>

#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>
#include <processing/evaluation.h>
#include <processing/rayTracingFacet.h>

#include <learning/learning.h>
#include <learning/learningMath.h>
#include <learning/learningRayTracing.h>
#include <learning/learningRayTracingGroundTruth.h>
#include <learning/learningIO_bin.h>


#ifdef Open3D
#include "open3d/Open3D.h"
#include "open3d/geometry/TetraMesh.h"
#endif

#include <CGAL/optimal_bounding_box.h>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;


/////////////////////////////////////////////////////////////////////
///////////////////////// CONTROL FUNCTIONS /////////////////////////
/////////////////////////////////////////////////////////////////////
int extractFeatures(dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& exportO){


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
        dir.write_file = dir.read_file;

    ///////////////////////////////
    ///////// IMPORT DATA /////////
    ///////////////////////////////
    // sampling input
    if(options.data_source == "scan"){
        if(!options.ground_truth){
            cout << "\nERROR: you need to load a ground truth polygon (with -g) for scanning it" << endl;
            return 1;
        };
        if(scanObjectClosed(data, options))
            return 1;
    }
    else if(options.data_source == "ply"){
        importPLYPoints(dir, data);
    }
    else if(options.data_source == "npz"){
        importNPZ(dir, data);
    }
//    #ifdef COLMAP
//    else if(options.data_source == "colmap"){
//        readColmapFiles(dir, data);
//    }
//    #endif
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


    // ground truth
    if(!dir.gt_poly_file.empty()){
        if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "off"){
            if(importOFFMesh(dir.path+dir.gt_poly_file, data.gt_poly))
                return 1;
//            CGAL::Polygon_mesh_processing::keep_largest_connected_components(data.gt_poly, 1);
        }
        #ifdef RECONBENCH
        else if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "mpu"){
            if(importImplicit(dir, data))
                return 1;
        }
        #endif
        else{
            cerr << "\nnot a valid ending for a ground truth file. choose either .mpu or .off" << endl;
            return 1;
        }
        assert(data.gt_poly.size_of_vertices() > 0);
        options.ground_truth = 1;
    }

    ///////////////////////////////////
    ///////// PREPROCESS DATA /////////
    ///////////////////////////////////
    if(!dir.transformation_file.empty()){
        // import translation matrix
//        if(importTransformationMatrix(dir,data))
//            return 1;
        if(applyTransformationMatrix(data))
            return 1;
    }
    // calc gt obb
    if(!options.gt_isclosed){
        cout << "\nMake bounding box of open GT for cropping input..." << endl;
        assert(options.ground_truth);
        // get centroid of ground truth
        getOrientedBoundingBox(data.gt_poly,data.bb_array,data.bb_surface_mesh);
        data.gt_centroid = getBoundingBoxCentroid(data.bb_array,1000.0);
        cout << "\t-ground truth centroid for ray tracing: " << data.gt_centroid << "-1000" << endl;
        clipWithBoundingBox(data.points,data.infos,data.bb_surface_mesh);

    }

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
    if(exportO.toply){
        exportO.color = data.has_color;
        exportO.normals = data.has_normal;
        exportPLY(dir, data.points, data.infos, exportO);
    }

    if(exportO.tonpz){
        toXTensor(data);
        exportNPZ(dir,data);
    }

    //////////////////////////////////////
    /////// DELAUNAY TRIANGULATION ///////
    //////////////////////////////////////
    if(options.Dt_epsilon > 0.0)
        makeAdaptiveDelaunayWithInfo(data, options.Dt_epsilon);
    else
        makeDelaunayWithInfo(data);

//    cout << "Mean edge length after scaling: " << calcMeanEdgeLength(data) << endl;


    ///////////////////////////////////
    /////// TETRAHEDRON SCORING ///////
    ///////////////////////////////////
    // make an index, necessary for graphExport e.g.
    options.make_global_cell_idx=1;
    indexDelaunay(data, options);

    if(!options.ground_truth){
        // this constructs the features which are exported to the cell ray graph,
        // important to note the difference between score, and features, where score is not used in any way in the learning.
        if(options.scale == 0.0)
            cout << "\nConsider turning on scaling with --sn, if your learning data is not yet scaled to a unit cube!" << endl;
        assert(data.has_sensor);
        learning::rayTracing(data.Dt, options);
        // this is only for making a reconstruction of a lrt, but has no influence on the features or the learning
//        learning::aggregateScoreAndLabel(data.Dt);
    }
    else{
        assert(data.has_sensor);
        learning::rayTracing(data.Dt, options);
        if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "off"){
            if(options.gt_isclosed){
                if(labelObjectWithClosedGroundTruth(dir,data, options.number_of_points_per_cell))
                    return 1;
            }
            else{
                if(labelObjectWithOpenGroundTruth(dir,data, options.number_of_points_per_cell))
                    return 1;
            }
        }
        #ifdef RECONBENCH
        if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "mpu"){
            labelObjectWithImplicit(data, options.number_of_points_per_cell);
        }
        #endif

    }





    ////////////////////////////
    ////////// EXPORT //////////
    ////////////////////////////


    //// export features and labels
    // check if labels directory exists, if not create it
    path lpath(dir.path);
    // this shitty problem here is not present on the laptop
    lpath /= string("gt");
    if(!is_directory(lpath))
        create_directory(lpath);
    // export the features and labels
    auto ge = graphExporter(dir, data.Dt, options);
    ge.run(true);


    ///// export the 3DT as npz
    options.make_global_cell_idx=0;
    options.make_finite_cell_idx=1;
    options.make_finite_vertex_idx=1;
    indexDelaunay(data, options);
    export3DT(dir,data);

    if(exportO.interface){
        for(auto cit = data.Dt.all_cells_begin(); cit != data.Dt.all_cells_end(); cit++){
            if(data.Dt.is_infinite(cit)){
                cit->info().gc_label = 1;
                continue;
            }
            // this means untraversed cells and 50/50 cells will be labelled as inside
            cit->info().gc_label = cit->info().outside_score > cit->info().inside_score ? 1 : 0;
        }
        exportInterface(dir, data, options, exportO);
    }

    return 0;
}



int main(int argc, char const *argv[]){


    cliParser ip("feat");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getFeat())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----FEATURE EXTRACTION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    if(extractFeatures(ip.dh, data, ip.ro, ip.eo))
        return 1;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----FEATURE EXTRACTION FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



