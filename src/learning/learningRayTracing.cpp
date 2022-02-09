#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <learning/learningRayTracing.h>
#include <util/geometricOperations.h>

using namespace std;

namespace learning{

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
int traverseCells(Delaunay& Dt,
                  Cell_handle& current_cell, Ray ray, double firstDist, double maxTraversalDist, int oppositeVertex,
                  int image, bool inside)
{
    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point


    if(!Dt.is_infinite(current_cell)){

        Point intersectionPoint;
        double dist;
        Point rayO = ray.source();
        Vector rayV = ray.to_vector();

        // push back the first intersection
        if(inside){
            current_cell->info().cb_facet_inside_first.push_back(firstDist);
            current_cell->info().fb_facet_inside[oppositeVertex].push_back(firstDist);
        }
        else{
            if(dist < maxTraversalDist){
                current_cell->info().cb_facet_outside_first.push_back(firstDist);
                current_cell->info().fb_facet_outside[oppositeVertex].push_back(firstDist);
            }
            else{
                current_cell->info().cb_facet_last_first.push_back(firstDist);
                current_cell->info().fb_facet_last[oppositeVertex].push_back(firstDist);
            }
        }

        // iterate over the faces of the current cell to find the second intersection
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int cellBasedVertexIndexz = (oppositeVertex+i)%4;
            Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndexz);
            bool intersect = rayTriangleIntersection(rayO, rayV, tri, intersectionPoint);
            // check if there is an intersection between the current ray and current triangle
            if(intersect){
                dist = sqrt(CGAL::squared_distance(intersectionPoint, rayO));
                // push back the second intersection
                if(inside){
                    current_cell->info().cb_facet_inside_second.push_back(dist);
                    current_cell->info().fb_facet_inside[cellBasedVertexIndexz].push_back(dist);
                    break; // stop the traversal in the second cell toward the inside
                }
                else{
                    if(dist < maxTraversalDist){
                        current_cell->info().cb_facet_outside_second.push_back(dist);
                        current_cell->info().fb_facet_outside[cellBasedVertexIndexz].push_back(dist);
                    }
                    else{
                        current_cell->info().cb_facet_last_second.push_back(dist);
                        current_cell->info().fb_facet_last[cellBasedVertexIndexz].push_back(dist);
                        break;
                    }
                }


                Facet mirror_fac = Dt.mirror_facet(std::make_pair(current_cell, cellBasedVertexIndexz));
                Cell_handle newCell = mirror_fac.first;
                int newIndex = mirror_fac.second;

                traverseCells(Dt,
                              newCell, ray, firstDist, maxTraversalDist, newIndex,
                              image, inside);
                break;
            }// end of intersect
        }// end of facet iteration
    }// end of finite check
    return 0;
}

void firstCell(Delaunay& Dt,
               Delaunay::Finite_vertices_iterator& vit, Ray ray, int image,
               bool inside){

    // make the inside ray
    if(inside)
        ray = ray.opposite();



    std::vector<Cell_handle> inc_cells;
    Dt.finite_incident_cells(vit, std::back_inserter(inc_cells));
    Point intersectionPoint;
    Point rayO = vit->point();
    Vector rayV = ray.to_vector();
    // set max traversal dist as sensor - point - distance
    double maxTraversalDist = sqrt(CGAL::squared_distance(rayO, vit->info().sensor_positions[image]));
    for(std::size_t i=0; i < inc_cells.size(); i++){

        Cell_handle current_cell = inc_cells[i];
        int cellBasedVertexIndex = current_cell->index(vit);
        Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
        // check if there is an intersection between the current ray and current triangle
        bool intersect = rayTriangleIntersection(rayO, rayV, tri, intersectionPoint);
        if(intersect){
            double dist = sqrt(CGAL::squared_distance(intersectionPoint, rayO));
            // the if stops the traversal once I hit the cell with the sensor center
            if(inside){
                current_cell->info().cb_vertex_inside.push_back(dist);
                current_cell->info().fb_vertex_inside[cellBasedVertexIndex].push_back(dist);
            }
            else{
                if(dist < maxTraversalDist){
                    current_cell->info().cb_vertex_outside.push_back(dist);
                    current_cell->info().fb_vertex_outside[cellBasedVertexIndex].push_back(dist);
                }
                else{
                    current_cell->info().cb_vertex_last.push_back(dist);
                    current_cell->info().fb_vertex_last[cellBasedVertexIndex].push_back(dist);
                    break;
                }
            }

            Cell_handle newCell = current_cell->neighbor(cellBasedVertexIndex);
            // the new cell can be accessed with neighbor function
            // the new cellBasedVertexIndex for starting point in traverseCells is the same as the current cellBasedVertexIndex
            // see: https://doc.cgal.org/latest/TDS_3/index.html#fig__TDS3figcomborient
            traverseCells(Dt, newCell, ray, dist, maxTraversalDist, cellBasedVertexIndex, image, inside);
            break;
        }// end of intersection result
    } // end of going around all the cells of a vertex
} // end of first cell



void rayTracing(Delaunay& Dt, runningOptions ro){



    int max_number_of_images = ro.number_of_rays;
    cout << "\nStart learning ray tracing..." << endl;
    cout << "\t-Trace " << max_number_of_images << " ray(s) to every point" << endl;

    auto start = std::chrono::high_resolution_clock::now();

    Delaunay::Finite_vertices_iterator vit;
    for(vit = Dt.finite_vertices_begin() ; vit != Dt.finite_vertices_end() ; vit++){
        int current_number_of_images = vit->info().sensor_positions.size();
        int number_of_images;
        if(max_number_of_images < 0 || max_number_of_images > current_number_of_images){
            number_of_images = current_number_of_images;
        }
        else{
            number_of_images = max_number_of_images;
        }
        for(int r = 0; r < number_of_images; r++){

            // construct ray from input point to sensor
            Ray ray(vit->point(), vit->info().sensor_positions[r]);

            // collect outside votes
            firstCell(Dt, vit, ray, r, 0);
            // collect inside votes
            firstCell(Dt, vit, ray, r, 1);
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "\t-Ray tracing done in " << full_duration.count() << "s" << std::endl;
}


// end of namespace learning
}






