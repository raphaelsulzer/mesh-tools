#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <learning/learningRayTracingGroundTruth.h>
#include <util/vectorArithmetic.h>
#include <processing/rayTracingTet.h>

using namespace std;

namespace learning{


////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
int nextGtCells(Delaunay& Dt,
                  Cell_handle& current_cell, Ray ray, vertex_info current_info, int oppositeVertex,
                  int current_image, bool inside)
{
    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // source               The "Delaunay" point of the current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point

    if(!Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int cellBasedVertexIndex = (oppositeVertex+i)%4;
            Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
            Point intersectionPoint;
            Point rayO = ray.source();
            Vector rayV = ray.to_vector();
            bool result = processing::rayTriangleIntersection(rayO, rayV, tri, intersectionPoint);
            // check if there is an intersection between the current ray and current triangle
            if(result){

                // add features to cell for learning
                if(inside)
                    current_cell->info().gt_inside += 1;
                else
                    current_cell->info().gt_outside += 1;

                // the if stops the traversal once I hit the cell with the sensor center
                if(!Dt.tetrahedron(current_cell).has_on_positive_side(current_info.sensor_positions[current_image])){
                    // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                    Facet mirror_fac = Dt.mirror_facet(std::make_pair(current_cell, cellBasedVertexIndex));
                    Cell_handle newCell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    nextGtCells(Dt,
                                  newCell, ray, current_info, newIdx,
                                  current_image, inside);
                }
                break;
            }
        }
    }
    return 0;
}

void firstGtCell(Delaunay& Dt, Point current_point, vertex_info current_info, int current_image,
               bool inside){

    // make the (inside) ray
    Ray ray(current_point, current_info.sensor_positions[current_image]);
    if(inside)
        ray = ray.opposite();

    // i already know that the current_cell intersects with inside and outside ray of the current vertex
    Cell_handle current_cell = Dt.locate(current_point);
    Facet current_facet;
    Triangle tri;
    Point intersectionPoint;
    Vector rayV;

    // iterate over the facets of this cell
    for(std::size_t i=0; i < 4; i++){

        current_facet = std::make_pair(current_cell, i);
        if(Dt.is_infinite(current_facet))
            continue;
        tri = Dt.triangle(current_facet);
        rayV = ray.to_vector();

        // check if there is an intersection between the current ray and current triangle
        bool intersect = processing::rayTriangleIntersection(current_point, rayV, tri, intersectionPoint);
        if(intersect){

            ///////// calculate the score for the current cell based on the distance /////////
            // add features to cell for learning
            if(inside)
                current_cell->info().gt_inside_first += 1;
            else
                current_cell->info().gt_outside_first += 1;

            // go to next cell
            if(!inside){
                Facet mirror_fac = Dt.mirror_facet(current_facet);
                Cell_handle newCell = mirror_fac.first;
                int newIdx = mirror_fac.second;
                // the if stops the traversal once I hit the cell with the sensor center
                if(!Dt.tetrahedron(current_cell).has_on_positive_side(current_info.sensor_positions[current_image]))
                    nextGtCells(Dt,newCell, ray, current_info, newIdx, current_image, inside);
                // TODO: maybe put an else here and increase the score of this cell even further, because it should DEFINITELY be outside
            }
            break;
        }// end of intersection result
    }// end of going around the facets of a cell
} // end of first cell

void countGtRays(dataHolder& data){

    std::cout << "\nVerify cells with ground truth scan..." << std::endl;
    assert(data.gt_points.size() > 0);
    assert(data.gt_points.size() == data.gt_infos.size());

    auto start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < data.gt_points.size(); i++){
        Point current_point = data.gt_points[i];
        vertex_info current_info = data.gt_infos[i];
        // collect outside votes
        firstGtCell(data.Dt, current_point, current_info, 0, 0);
        // collect inside votes
        firstGtCell(data.Dt, current_point, current_info, 0, 0);

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "\t-Ray tracing done in " << full_duration.count() << "s" << std::endl;
}



// end of namespace learning
}






