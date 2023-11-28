#include <IO/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <IO/ttIO.h>
#include <IO/ethIO.h>
#ifdef COLMAP
#include <IO/colmapfileIO.h>
#endif
#include <util/helper.h>
#include <util/geometricOperations.h>

#include <processing/meshProcessing.h>
#include <processing/edgeManifoldness.h>
#include <processing/graphCut.h>

#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>
#include <processing/evaluation.h>

#include <learning/learning.h>
#include <learning/learningMath.h>
#include <learning/rayTracing.h>
#include <learning/learningIO_bin.h>


#include <CGAL/optimal_bounding_box.h>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;


/////////////////////////////////////////////////////////////////////
///////////////////////// CONTROL FUNCTIONS /////////////////////////
/////////////////////////////////////////////////////////////////////
int extractFeatures(dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& exportO){

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////// LOAD DATA /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // SAMPLING INPUT
    if(options.data_source == "ply"){
        importPLYPoints(dir, data);
    }
    else if(options.data_source == "npz"){
        importNPZPoints(dir, data);
    }
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

    // GROUND TRUTH INPUT
    if(!dir.gt_poly_file.empty()){
        if(dir.gt_poly_file.substr(dir.gt_poly_file.length() - 3 ) == "off"){
            if(importOFFMesh(dir.path+dir.gt_poly_file, data.gt_poly))
                return 1;
        }
        else{
            cerr << "\nnot a valid ending for a ground truth file. choose either .mpu or .off" << endl;
            return 1;
        }
        assert(data.gt_poly.size_of_vertices() > 0);
        options.ground_truth = 1;
    }

    // EVAL POINTS INPUT (list of uniformly distributed points in bbox; ConvONet style)
    if(!dir.occ_file.empty()){
        if(importOccPoints(dir,data))
            return 1;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////// PROCESS DATA ////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // DELAUNAY TRIANGULATION
    if(options.Dt_epsilon > 0.0)
        makeAdaptiveDelaunayWithInfo(data, options.Dt_epsilon);
    else
        makeDelaunayWithInfo(data);

    // INDEX DELAUNAY, necessary for graphExport e.g.
    options.make_global_cell_idx=1;
    options.make_global_vertex_idx=1;
    options.make_finite_cell_idx=0;
    options.make_finite_vertex_idx=0;
    indexDelaunay(data, options);

    // RAY TRACING
    assert(data.has_sensor);
    learning::rayTracing(data.Dt, options);

    // TETRAHEDRON SCORING
    if(options.ground_truth){
        if(options.gt_type == "closed"){
            if(labelObject(dir,data, options.number_of_points_per_cell))
                return 1;
        }
        else if(options.gt_type == "real"){
            if(labelReal(dir,data, options.number_of_points_per_cell))
                return 1;
        }
        else if(options.gt_type == "srd"){
            // label the synthetic room dataset
            if(labelSRD(dir,data, options.number_of_points_per_cell))
                return 1;
        }
        else{
            cout << "\nERROR: not a valid type for ground truth mesh." << endl;
            return 1;
        }
    }

    // GET TET INDEX PER EVAL POINT
    if(!dir.occ_file.empty())
        point2TetraIndex(data);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////// EXPORT //////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // CHECK IF DGNN DIRECTORY EXISTS, if not create it
    path lpath(dir.path);
    lpath /= string("dgnn");
    if(!is_directory(lpath))
        create_directory(lpath);

    // EXPORT _eval.npz file, with eval points, occupancy and point2tetIndex
    if(!dir.occ_file.empty())
        exportOccPoints(dir,data);

    // EXPORT features (and labels)
    auto ge = graphExporter(dir, data.Dt, options);
    ge.run(true);

    // EXPORT the 3DT as npz
    export3DT(dir,data);



    // EXPORT interface from ground truth labels, cameras and/or points
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
    // export cell scores
    if(exportO.cellScore){
        exportCellCenter(dir, data.Dt);
        exportCellScore(dir, data.Dt);
    }
    if(exportO.delaunay){
        export3DT(dir, data.Dt);
    }
    // cameras
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
    if(!ip.dh.path.empty())
        cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    if(extractFeatures(ip.dh, data, ip.ro, ip.eo))
        return 1;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----FEATURE EXTRACTION FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



