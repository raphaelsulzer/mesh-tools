#include <IO/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <processing/meshProcessing.h>
#include <processing/evaluation.h>



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////// CONTROL FUNCTIONS /////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




int main(int argc, char const *argv[]){


    cliParser ip("eval");
    if(ip.parse(argc, argv))
        return 1;
    ip.getInput();
    ip.getEval();
    ip.getOutput();


    auto start = std::chrono::high_resolution_clock::now();

    cout << "\n-----EVALUATE MESH-----\n" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    meshProcessingOptions options;
    options.number_of_components_to_keep = ip.ro.number_of_components_to_keep;
    options.try_to_close = ip.ro.try_to_close;
    options.try_to_make_manifold = ip.ro.try_to_make_manifold;

    // points to sample
    int sample_points = ip.eo.sampling_method_option;


    // read file depending on file type
    if(ip.dh.read_file_type == "off")
        importOFFMesh(ip.dh, data.smesh);
    else if(ip.dh.read_file_type == "ply"){
        importPLYMesh(ip.dh, data.smesh);
    }
    else{
        cout << "File type not support. Supported are .off and .ply files\n" << endl;
        return 1;
    }

    if(!ip.dh.gt_poly_file.empty()){
        if(ip.dh.gt_poly_file.substr(ip.dh.gt_poly_file.length() - 3 ) == "off"){
            if(importOFFMesh(ip.dh.path+ip.dh.gt_poly_file, data.gt_poly))
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
        if(calcIOU(data.gt_poly,data.smesh, sample_points, iou))
            return 1;
        printMeshEvaluation(ip.dh,data,iou);
    }
    else
        printMeshEvaluation(ip.dh,data,-1.0);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----EVALUATE MESH FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



