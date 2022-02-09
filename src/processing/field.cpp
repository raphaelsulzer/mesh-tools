#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <util/helper.h>
#include <util/geometricOperations.h>


#include <processing/field.h>


#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <boost/foreach.hpp>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <util/helper.h>

#include <random>


using namespace std;


int makeFieldOnDelaunayVertices(dataHolder& data, string mode){

    auto start = std::chrono::high_resolution_clock::now();

    cout << "\nMaking field..." << endl;
    cout << "\t-using mode: " << mode << endl;

    if(data.gt_poly.size_of_vertices() == 0){
        cout << "\nERROR: you need to load a ground truth polygon (with -g) for scanning it" << endl;
        return 1;
    }
    if(!data.gt_poly.is_closed()){
        cout << "\nERROR: ground truth is not closed" << endl;
        return 1;
    }

    // Initialize the point-in-polyhedron tester
    // Construct AABB tree with a KdTree
    SurfaceMesh smesh;
    CGAL::copy_face_graph(data.gt_poly, smesh);
    Tree tree(faces(smesh).first, faces(smesh).second, smesh);
    tree.accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    Point p;
    double val;
    Delaunay::Finite_vertices_iterator aci;

    const Point_inside inside_tester(tree);

    for(aci = data.Dt.finite_vertices_begin(); aci != data.Dt.finite_vertices_end(); aci++){

        p = aci->point();

        if(inside_tester(p) == CGAL::ON_BOUNDED_SIDE){
            val=-1.0;
        }
        else if(inside_tester(p) == CGAL::ON_UNBOUNDED_SIDE){
            val=1.0;
        }
        else{
            val=0.0;
        }

        if(mode == "occ"){
            aci->info().field_value = val;
        }
        else if(mode == "sdf"){
            aci->info().field_value = sqrt(tree.squared_distance(p))*val;
        }
        else{
            cout << "\nERROR: not a valid field mode!" << endl;
            return 1;
        }
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-in "<<duration.count() << "s" << endl;

    return 0;
}
