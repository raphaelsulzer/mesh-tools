#include <exe/sure.h>
#include <exe/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <IO/ttIO.h>
#include <IO/ethIO.h>
#ifdef OpenMVS
#include <IO/mvIO.h>
#endif
#ifdef COLMAP
#include <IO/colmapIO.h>
#endif
#include <util/helper.h>
#include <util/vectorArithmetic.h>

#include <processing/meshProcessing.h>
#include <processing/edgeManifoldness.h>
#include <processing/graphCut.h>

#include <processing/pointSetProcessing.h>
#include <processing/normalAndSensorProcessing.h>
#include <processing/evaluation.h>
#include <processing/rayTracingFacet.h>
#include <processing/field.h>

#include <learning/learning.h>
#include <learning/learningMath.h>
#include <learning/learningRayTracing.h>
#include <learning/learningRayTracingGroundTruth.h>
#include <learning/learningIO.h>

#include <CGAL/refine_mesh_3.h>

#ifdef Open3D
#include "open3d/Open3D.h"
#include "open3d/geometry/TetraMesh.h"
#endif

#include <CGAL/optimal_bounding_box.h>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;

#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <processing/field.h> // for the tree defitition


typedef Delaunay::Geom_traits GT;
typedef GT::FT FT;

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
//typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::FT FT;
//typedef FT (*Function)(Point);
typedef std::function<FT (Point)> Function;
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

FT signed_distance(Tree& tree, int& call_counter, Point p){

    const Point_inside inside_tester(tree);

    call_counter++;

    if(inside_tester(p) == CGAL::ON_BOUNDED_SIDE){
        return -1.0*sqrt(tree.squared_distance(p));
    }
    else if(inside_tester(p) == CGAL::ON_UNBOUNDED_SIDE){
        return 1.0*sqrt(tree.squared_distance(p));
    }
    else{
        return 0.0;
    }

}

FT occupancy(Tree& tree, int& call_counter, Point p){

    const Point_inside inside_tester(tree);

    call_counter++;

    if(inside_tester(p) == CGAL::ON_BOUNDED_SIDE){
        return -1.0;
    }
    else if(inside_tester(p) == CGAL::ON_UNBOUNDED_SIDE){
        return 1.0;
    }
    else{
        return 0.0;
    }

}



int exportBoissonnatSurface(dirHolder dir, dataHolder& data, runningOptions options){

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "\nExtract Isosurface with Boissonnat..." << endl;
    cout << "\t-use field " << options.field << endl;
    cout << "\t-lower bound facet angle: " << options.options[0] << endl;
    cout << "\t-upper bound radius surface Delaunay balls: " << options.options[1] << endl;
    cout << "\t-upper bound for facet center-center distances: " << options.options[2] << endl;

    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation


    SurfaceMesh smesh;
    CGAL::copy_face_graph(data.gt_poly, smesh);
    Tree tree(faces(smesh).first, faces(smesh).second, smesh);
    tree.accelerate_distance_queries();


    // idea to use std::bind to pass additional arguments to the oracle function is from here:
    //http://www.alecjacobson.com/weblog/?p=4034&cpage=1#comment-1028470
    // additionally, needed to add the std::ref() to the arguments, because bind automatically copies or moves its arguments
    // see here: https://stackoverflow.com/questions/26187192/how-to-bind-function-to-an-object-by-reference
    // and since AABBTree cannot be copied or moved it doesn't even compile without it
    // alternative would be to use a lambda expresion instead of std::bind():
    // [&tree](Point p)->FT{ return signed_distance(tree,p);},

    int counter=0;
    std::function<FT (Point)> func;

    if(options.field == "sdf"){
        func = std::bind(&signed_distance,std::ref(tree),std::ref(counter),std::placeholders::_1);
    }
    else if(options.field == "occ"){
        func = std::bind(&occupancy,std::ref(tree),std::ref(counter),std::placeholders::_1);
    }
    else{
        cout << "\nERROR: not a valid field mode!" << endl;
        return 1;
    }

    Surface_3 surface(
                    func,             // pointer to oracle function return signed distance or occupancy
                    Sphere_3(CGAL::ORIGIN, 75.0*75.0)); // bounding sphere with squared radius. using 40 here because reconbench objects are 75 mm cube, so max 37.5 radius
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(options.options[0],  // angular bound
                                                   options.options[1],  // radius bound
                                                   options.options[2]); // distance bound

    // criteria see here: https://doc.cgal.org/latest/Surface_mesher/classCGAL_1_1Surface__mesh__default__criteria__3.html
//    a lower bound on the minimum angle in degrees of the surface mesh facets.
//    an upper bound on the radius of surface Delaunay balls. A surface Delaunay ball is a ball circumscribing a facet, centered on the surface and empty of vertices. Such a ball exists for each facet of the current surface mesh. Indeed the current surface mesh is the Delaunay triangulation of the current sampling restricted to the surface which is just the set of facets in the three dimensional Delaunay triangulation of the sampling that have a Delaunay surface ball.
//    an upper bound on the center-center distances of the surface mesh facets. The center-center distance of a surface mesh facet is the distance between the facet circumcenter and the center of its surface Delaunay ball.
    //  topology	is the set of topological constraints which have to be verified by each surface facet. See section Delaunay Refinement for further details. Note that if one parameter is set to 0, then its corresponding criteria is ignored.

    cout << "\t-defined field" << endl;


    // meshing surface
    int init = 20;
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(),init);

    cout << "\t-made " << counter << " calls to oracle" << endl;
    std::ofstream out(dir.path+dir.write_file+"_"+options.field+"_"+to_string(counter)+".off");
    CGAL::output_surface_facets_to_off(out,c2t3);

    cout << "\t-make surface mesh" << endl;

//    SurfaceMesh sm;
//    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
//    out << sm << std::endl;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-done after " << duration.count() << "s" << endl;
}


#ifdef Open3D
#include "open3d/Open3D.h"
#include "open3d/geometry/TetraMesh.h"
void exportMITSurface(dirHolder& dir, Delaunay& Dt, bool mean, double iso_value){

    //// to Open3d tetrahmesh
    cout << "\nExport iso-mesh..." << endl;
    cout << "\t-at iso value " << iso_value << endl;
    cout << "\t-to " << dir.path+dir.write_file+dir.suffix << endl;
    vector<Eigen::Vector3d> points;
    vector<double> values;
    vector<Eigen::Vector4i, open3d::utility::Vector4i_allocator> tets;
    int i = 0;
    for(auto vit = Dt.finite_vertices_begin(); vit != Dt.finite_vertices_end(); vit++){

        if(mean){
            std::vector<Cell_handle> inc_cells;
            Dt.incident_cells(vit, std::back_inserter(inc_cells));
            auto cit = inc_cells.begin();
            double value = 0;
            double weight;
            int cell_count = 0;
            while(cit != inc_cells.end()){
                Cell_handle c = *cit;
    //            Point p1 = c->vertex(0)->point();
    //            Point p2 = c->vertex(1)->point();
    //            Point p3 = c->vertex(2)->point();
    //            Point p4 = c->vertex(3)->point();
    //            Point centroid = CGAL::centroid(p1,p2,p3,p4);
    //            weight = CGAL::squared_distance(centroid,vit->point());
    //            weight = exp(-weight*0.01);
    //            if(data.Dt.is_infinite(c))
    //                value+=10;
    //            else
    //                value+=c->info().outside_score*weight;

    //            double vol = data.Dt.tetrahedron(c).volume()/4.0;
    //            value+=c->info().outside_score*vol;
    //            value+=c->info().gc_label;
                value+=c->info().outside_score;
                cell_count++;
                cit++;
            }
            value/=cell_count;
            values.push_back(value);
        }
        else
            values.push_back(vit->info().field_value);

        points.push_back(cgal2Eigen(vit->point()));
        vit->info().global_idx = i++;
    }
    for(auto cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        Eigen::Vector4i tet;
        for(int i = 0; i < 4; i++){
            tet[i] = cit->vertex(i)->info().global_idx;
        }
        tets.push_back(tet);
    }
    auto tetmesh = open3d::geometry::TetraMesh(points,tets);

    auto start = std::chrono::high_resolution_clock::now();
    auto trimesh = tetmesh.ExtractTriangleMesh(values,iso_value);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "\t-TIME: MIT with 1M points in seconds: " << duration.count() << "ms" << endl;


    open3d::io::WriteTriangleMeshToPLY(dir.path+dir.write_file+dir.suffix, *trimesh, true, false, false, false, false, false);

//    exportIsoValues(dir, Dt, values);
}
#endif

int exportMINESurface(dirHolder dir, Delaunay& Dt){

    return 0;





}




int main(int argc, char const *argv[]){


    cliParser ip("iso");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getIso())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----ISOSURFACE EXTRACTION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;


    dataHolder data;

    // load a ground truth
    if(importOff(ip.dh.path+ip.dh.gt_poly_file, data.gt_poly))
        return 1;

    if(ip.ro.method == "boissonnat")
        exportBoissonnatSurface(ip.dh,data,ip.ro);
#ifdef Open3D
    else if(ip.ro.method == "mit"){
        // import a sampling of the ground truth for triangulation
        if(importLidarPoints(ip.dh, data))
            return 1;

        // TODO: add additional points in the bounding sphere of the ground truth
        if(ip.ro.add_points > 0){
            Bbox bb = CGAL::Polygon_mesh_processing::bbox(data.gt_poly);
            Point bb_center = Point((bb.xmin()+bb.xmax())/2,
                                    (bb.ymin()+bb.ymax())/2,
                                    (bb.zmin()+bb.zmax())/2);

            double radius = std::sqrt(pow(bb.xmin()-bb_center.x(),2)+
                             pow(bb.ymin()-bb_center.y(),2)+
                             pow(bb.zmin()-bb_center.z(),2));

            CGAL::Random random(42);
            CGAL::Random_points_in_sphere_3<Point> points_on_sphere(radius*1.1, random); // for aims
            CGAL::cpp11::copy_n(points_on_sphere, ip.ro.add_points, std::back_inserter(data.points));
            for(int i = 0; i < ip.ro.add_points; i++){
                vertex_info vinf;
                data.infos.push_back(vinf);
            }
        }

        makeDelaunayWithInfo(data);
        makeFieldOnDelaunayVertices(data,ip.ro.field);

        string n_calls = to_string(data.Dt.number_of_vertices());
        ip.dh.suffix = "_"+ip.ro.field+"_"+n_calls+".ply";
        exportMITSurface(ip.dh, data.Dt, 0, ip.ro.value);
    }
#endif
    else if(ip.ro.method == "mine")
        exportMINESurface(ip.dh, data.Dt);
    else{
        cout << "\nERROR: not a valid method for isosurface extraction." << endl;
        return 1;
    }

    if(ip.eo.scan){
        ip.eo.normals = false;
        ip.eo.sensor_vec = false;
        ip.eo.sensor_position = false;
        ip.eo.color = false;
        ip.dh.suffix = "_scan";
        exportPLY(ip.dh, data.points, data.infos, ip.eo);
    }
    if(ip.eo.convexHull){
        exportConvexHull(ip.dh, data.Dt);
    }


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----ISOSURFACE EXTRACTION in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



