#include <IO/inputParser.h>
#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <learning/learningIO_bin.h>
#include <processing/meshProcessing.h>
#include <processing/graphCut.h>

int main(int argc, char const *argv[]){


    cliParser ip("occ2mesh");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getOcc2Mesh())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----OCC2MESH-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////// IMPORT SCAN ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // sampling input
    if(ip.ro.data_source == "ply"){
        importPLYPoints(ip.dh, data);
    }
    else if(ip.ro.data_source == "npz"){
        importNPZPoints(ip.dh, data);
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////// DELAUNAY TRIANGULATION /////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(ip.ro.Dt_epsilon > 0.0)
        makeAdaptiveDelaunayWithInfo(data, ip.ro.Dt_epsilon);
    else
        makeDelaunayWithInfo(data);



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////// PREDICTION ///////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(importPrediction(ip.dh,data,ip.ro))
        return 1;




    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////// OPTIMIZATION ///////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(ip.ro.optimization == 1){
        ip.ro.gco_iterations = -1;
        graphCutTet(data, ip.ro);
        cout << "\t-Infinite cells forced outside? " << ip.ro.closed_prior << endl;
    }
    else{
        cout << "\nLabel cells without optimization by taking max score..." << endl;
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


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////// EXPORT //////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(ip.eo.interface)
        exportInterface(ip.dh, data, ip.ro, ip.eo);

    #ifdef OpenMVS
    if(options.clean_mesh)
        omvsCleanMesh(dir,data);
    #endif

//    if(exportO.coloredFacets)
//        exportColoredFacets(dir, data.Dt, options.optimization);
    if(ip.eo.cellScore){
        exportCellCenter(ip.dh, data.Dt);
        exportCellScore(ip.dh, data.Dt);
    }
    if(ip.eo.convexHull){
        exportConvexHull(ip.dh, data.Dt);
    }
    ////////// create surface mesh ////////////
    // set export options
    meshProcessingOptions mpOptions;
    mpOptions.try_to_close = ip.ro.try_to_close;
    mpOptions.try_to_make_manifold = ip.ro.try_to_make_manifold;
    mpOptions.number_of_components_to_keep = ip.ro.number_of_components_to_keep;
    mpOptions.factor_for_removing_large_faces = ip.ro.factor_for_removing_large_faces;
    if(ip.eo.mesh){
        ip.dh.suffix = "_mesh";
        createSurfaceMesh(data, mpOptions);
        exportPLY(ip.dh, data.smesh);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----OCC2MESH FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



