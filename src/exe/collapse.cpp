#include <IO/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

#include <util/helper.h>
#include <util/geometricOperations.h>
#include <processing/meshProcessing.h>
#include <processing/evaluation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////// CONTROL FUNCTIONS /////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




int main(int argc, char const *argv[]){


    cliParser ip("collapse");
    if(ip.parse(argc, argv))
        return 1;
    ip.getInput();
    ip.getCollapse();
    ip.getOutput();


    auto start = std::chrono::high_resolution_clock::now();

    cout << "\n-----EDGE COLLAPSE-----\n" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    meshProcessingOptions options;
    options.keep_edges = ip.ro.keep_edges;
    options.number_of_components_to_keep = ip.ro.number_of_components_to_keep;
    options.try_to_close = ip.ro.try_to_close;

    cout << "keep edges " << options.keep_edges << endl;

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

    edgeCollapse(data, options);

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
        int sample_points = ip.eo.sampling_method_option;
        calcIOU(data.gt_poly,data.smesh, sample_points, iou);
        printMeshEvaluation(ip.dh,data,iou);
    }
    else
        printMeshEvaluation(ip.dh,data,-1.0);

    if(ip.dh.write_file.empty()){
        ip.dh.write_file = ip.dh.read_file;
        ip.dh.suffix = double2string(ip.ro.keep_edges);
    }
    exportOFF(ip.dh, data.smesh);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----EDGE COLLAPSE FINISHED in "<< duration.count() << "s -----\n" << endl;

}



