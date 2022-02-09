#include "poisson/poisson.h"
#include "base/cgal_typedefs.h"
#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>
#include <IO/fileIO.h>
#include <util/helper.h>

#include <CGAL/poisson_surface_reconstruction.h>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// POISSON RECONSTRUCTION ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
PoissonReconstructor::PoissonReconstructor(dirHolder a, exportOptions b)
{
    dir = a;
    eo = b;
}


void PoissonReconstructor::makePoissonMesh(dataHolder& data){
    auto start = chrono::high_resolution_clock::now();

    cout << "\nPoisson surface reconstruction of " << data.points.size() << " points" << endl;

    double average_spacing = CGAL::compute_average_spacing<Concurrency_tag>
      (data.pVP, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()));
    cout << "\t-computed average spacing of " << average_spacing << endl;
    // TODO: split this into seperate functions, so I can pass the Manifold_with_boundary_tag
    // https://doc.cgal.org/latest/Surface_mesher/group__PkgSurfaceMesher3FunctionsMakeMesh.html#Cross_link_anchor_1509

    CGAL::poisson_surface_reconstruction_delaunay
        (data.pVP.begin(), data.pVP.end(),
         CGAL::First_of_pair_property_map<PointVectorPair>(),
         CGAL::Second_of_pair_property_map<PointVectorPair>(),
         data.poly, average_spacing);

    // TODO: https://doc.cgal.org/4.14.3/Surface_mesher/Surface_mesher_2mesh_an_implicit_function_8cpp-example.html
    // add this to the function above, if it complains that point is already inserted
    // CGAL::Non_manifold_tag()

    auto stop = chrono::high_resolution_clock::now();
    auto full_duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done in " << full_duration.count() << "s" << endl;
}


void PoissonReconstructor::run(runningOptions options){


    ////////////// READ DATA ///////////////
    dataHolder data;
    importLidarPoints(dir, data);

    if(options.standardize)
        standardizePointSet(dir, data);

    if(options.subsample_n_points > 0)
        sampleFun(data.points, options.subsample_n_points);

    if(options.subsample_grid_spacing > 0)
        gridSimplify(data.points, options.subsample_grid_spacing);

    // estimate normals including estimating neighborhood size and normal orientation
    estimateNormals(data,0,2);


    if(eo.scan){
        dir.write_file = dir. read_file;
        dir.suffix = "_normals";
        eo.normals = true;
        exportPoints(dir, data, eo);
    }

    makePoissonMesh(data);

    int erased_components = data.poly.keep_largest_connected_components(1);
    exportOFF(data.poly, dir.path+dir.read_file+"_poisson_1");


}
