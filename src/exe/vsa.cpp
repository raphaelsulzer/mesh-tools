#include <exe/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

#include <util/helper.h>
#include <util/vectorArithmetic.h>
#include <processing/meshProcessing.h>
#include <processing/evaluation.h>

#include "boost/tuple/tuple.hpp"


#include <random>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
namespace VSA = CGAL::Surface_mesh_approximation;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor        face_descriptor;



//void circulateAroundFacet()



void printMDL(dataHolder& data){


    typedef boost::graph_traits<SurfaceMesh>::face_descriptor        face_descriptor;
    typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor      vertex_descriptor;

    SurfaceMesh::Property_map<face_descriptor,CGAL::Color> colors;
    SurfaceMesh::Property_map<face_descriptor,Vector> fnormals;
    SurfaceMesh::Property_map<face_descriptor,bool> processed;


    colors = data.smesh.add_property_map<face_descriptor, CGAL::Color>("f:color", CGAL::white()).first;
    fnormals = data.smesh.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;
    processed = data.smesh.add_property_map<face_descriptor, bool>("f:processed", false).first;
    // compute normals
    PMP::compute_face_normals(data.smesh, fnormals);

    array<CGAL::Color,7> color_palette{CGAL::red(), CGAL::blue(), CGAL::green(), CGAL::orange(), CGAL::gray(), CGAL::purple(), CGAL::violet()};
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> distrib(0, 6);



    int coplanar = 0;
    int skipped = 0;

    for(auto f : data.smesh.faces()){

        if(processed[f])
            continue;

        // get first opposite halfedge
        SurfaceMesh::Halfedge_index first_opposite_halfedge = data.smesh.opposite(data.smesh.halfedge(f));
        SurfaceMesh::Halfedge_index next_opposite_halfedge = first_opposite_halfedge;

        do{
            SurfaceMesh::Face_index current_neighbor = data.smesh.face(next_opposite_halfedge);
            if (current_neighbor == SurfaceMesh::null_face()){
                skipped++;
                next_opposite_halfedge = data.smesh.next(next_opposite_halfedge);
            }
            else{
//                cout << CGAL::scalar_product(fnormals[current_neighbor],fnormals[f]) << endl;
                if(0.99999999999 < CGAL::scalar_product(fnormals[current_neighbor],fnormals[f]) < 1.00000000001){
                    colors[f] = color_palette[distrib(gen)];
                    colors[current_neighbor] = color_palette[distrib(gen)];
                    coplanar++;
                    processed[f] = true;
                    processed[current_neighbor] = true;
                    // TODO:
                    // maybe the result is correct, because maybe every single face has a neighbor which is coplanar, would make sense as an output of the algorithm
                    // so what I actually need to do is do a region growing or a recursive neighborhood search and color connected planar components in a random CGAL::Color();
                    // for random color, simply define a couple of cgal colors, put them in an array and take a random element of the array for coloring
                }
                next_opposite_halfedge = data.smesh.next(next_opposite_halfedge);
            }
        }
        while(first_opposite_halfedge != next_opposite_halfedge);

    }
    cout << "skipped faces: " << skipped << endl;
    cout << "coplanar faces: " << coplanar << endl;

}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////// CONTROL FUNCTIONS /////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char const *argv[]){


    cliParser ip("vsa");
    if(ip.parse(argc, argv))
        return 1;
    ip.getInput();
    ip.getVsa();
    ip.getOutput();


    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----VARIATIONAL SHAPE APPROXIMATION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
    meshProcessingOptions options;
    options.number_of_components_to_keep = 1;
    options.try_to_close = 1;

    // read file depending on file type
    if(ip.dh.read_file_type == "off")
        importOff(ip.dh, data.smesh);
    else if(ip.dh.read_file_type == "ply"){
        importLidarMesh(ip.dh, data);
        createSurfaceMesh(data, options);
        data.points.clear();
    }
    else{
        cerr << "File type not specified or not support. Supported are .off and .ply files" << endl;
        return 1;
    }


//    printMDL(data);
//    exportPLY_CGAL(ip.dh,data.smesh);

    options.max_number_of_proxies = ip.ro.max_number_of_proxies;
    options.number_of_components_to_keep = ip.ro.number_of_components_to_keep;
    options.try_to_close = ip.ro.try_to_close;
    variationalShapeApproximation(data, options);
//    printMDL(data);

    if(!ip.dh.gt_poly_file.empty()){
        if(ip.dh.gt_poly_file.substr(ip.dh.gt_poly_file.length() - 3 ) == "off"){
            if(importOff(ip.dh.path+ip.dh.gt_poly_file, data.gt_poly))
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
        if(calcIOU(data, iou))
            return 1;
        printMeshEvaluation(ip.dh,data,iou);
    }
    else
        printMeshEvaluation(ip.dh,data,-1.0);

    if(ip.dh.write_file.empty()){
        ip.dh.write_file = ip.dh.read_file + "_vsa";
        ip.dh.suffix = double2string(ip.ro.max_number_of_proxies);
    }
//    exportOFF(ip.dh, data.smesh);
    exportPLY(ip.dh,data.smesh);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----VARIATIONAL SHAPE APPROXIMATION FINISHED in "<< duration.count() << "s -----\n" << endl;

}



