#include <IO/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>

#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <boost/filesystem.hpp>
using namespace boost::filesystem;



int main(int argc, char const *argv[]){


    cliParser ip("clf");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----MAKE AN ORIENTED BOUNDING BOX AND SAVE AS PLY-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;

    importLidarPoints(ip.dh,data);
    std::array<Point, 8> obb_points;
    CGAL::oriented_bounding_box(data.points,obb_points);
    SurfaceMesh obb_sm;
    if(ip.ro.scale > 0.0){
        double mx = 0, my = 0, mz = 0;
        for(const auto p : obb_points){
            mx+=p.x();
            my+=p.y();
            mz+=p.z();
        }
        Point centroid(mx/8,my/8,mz/8);
        std::array<Point, 8> obb_scaled;
        int i = 0;
        for(const auto p : obb_points){
            Vector mover(ip.ro.scale*(p.x() - centroid.x()),
                         ip.ro.scale*(p.y() - centroid.y()),
                         ip.ro.scale*(p.z() - centroid.z()));
            obb_scaled[i++] = p + mover;
        }
        CGAL::make_hexahedron(obb_scaled[0], obb_scaled[1], obb_scaled[2], obb_scaled[3],
                              obb_scaled[4], obb_scaled[5], obb_scaled[6], obb_scaled[7], obb_sm);

    }
    else{
        CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                              obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb_sm);
    }
    CGAL::Polygon_mesh_processing::triangulate_faces(obb_sm);


    if(ip.dh.write_file.empty())
        ip.dh.write_file = ip.dh.read_file;
    std::ofstream out(ip.dh.path+ip.dh.write_file+"_obb.ply");
    CGAL::write_ply(out,obb_sm);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----MAKE AN ORIENTED BOUNDING BOX AND SAVE AS PLY FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



