#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <labatut/rayTracingFacet.h>
#include <util/geometricOperations.h>

using namespace std;


///////// Ray tracing implemented according to "Robust and efficient surface reconstruction from range data" by Labatut et al. 2009:
/// Line-of-sight (outside) and ray (inside) information is accumulated on the oriented facets of a 3D Delaunay triangulation.
/// Only parameter is sigma, which controls the weight per ray and the maximal inside traversal at the same time.


namespace processing{


void RayCaster::traverseOutside(Vertex_handle vit, Cell_handle& current_cell,  int oppositeVertex, int sensor_idx){
    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // source               The "Delaunay" point of the current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point


    if(!data_.Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int cellBasedVertexIndex = (oppositeVertex+i)%4;
            Triangle tri = data_.Dt.triangle(current_cell, cellBasedVertexIndex);
            Point intersectionPoint;
            Point source = vit->point();
            Vector rayV = Vector(vit->info().sensor_positions[sensor_idx] - source);
            bool intersection = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            if(intersection){
                // the if stops the traversal once I hit the cell with the sensor center
                if(!data_.Dt.tetrahedron(current_cell).has_on_positive_side(vit->info().sensor_positions[sensor_idx])){
                    double dist2 = CGAL::squared_distance(intersectionPoint, source);
                    // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                    Facet mirror_fac = data_.Dt.mirror_facet(std::make_pair(current_cell, cellBasedVertexIndex));
                    Cell_handle next_cell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
//                    cout << "new index " << newIdx << " previous index " << cellBasedVertexIndex << endl;
                    // set the facet weight
//                    current_cell->info().facet_weights[cellBasedVertexIndex] +=(alpha * ComputeDistanceProb(dist2));
                    current_cell->info().facet_weights[cellBasedVertexIndex] += ComputeDistanceProb(dist2);
                    traverseOutside(vit, next_cell, newIdx, sensor_idx);
                }
                else
                    current_cell->info().outside_score+=alpha_vis;
                break;
            }
        }
    }
    // put outside score of infinite cell very high
    else{
        // TODO: what I am doing here is not correct, the if statement should be activated, but it crashes if it is
//        if(data_.Dt.tetrahedron(current_cell).has_on_positive_side(vit->info().sensor_positions[sensor_idx])){
        // TODO: that must be the problem!!!!!!!!!
            current_cell->info().outside_score+=alpha_vis;
//        }
    }
    return;
}


void RayCaster::outside(Vertex_handle vit, int sensor_idx){

    std::vector<Cell_handle> inc_cells;
    data_.Dt.incident_cells(vit, std::back_inserter(inc_cells));

    for(std::size_t i=0; i < inc_cells.size(); i++){
        Cell_handle current_cell = inc_cells[i];

        if(!data_.Dt.is_infinite(current_cell))
        {
            int cellBasedVertexIndex = current_cell->index(vit);

            Triangle tri = data_.Dt.triangle(current_cell, cellBasedVertexIndex);
            Point intersectionPoint;
            Point source = vit->point();
            Vector rayV = Vector(vit->info().sensor_positions[sensor_idx] - source);

            // check if there is an intersection between the current ray and current triangle
            bool intersect = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            if(intersect){
                if(!data_.Dt.tetrahedron(current_cell).has_on_positive_side(vit->info().sensor_positions[sensor_idx])){
                    double dist2 = CGAL::squared_distance(intersectionPoint, source);
                    Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                    Facet mirror_fac = data_.Dt.mirror_facet(fac);
                    Cell_handle next_cell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    // set the facet weight
//                    current_cell->info().facet_weights[cellBasedVertexIndex] +=(alpha * ComputeDistanceProb(dist2));
                    current_cell->info().facet_weights[cellBasedVertexIndex] += ComputeDistanceProb(dist2);

                    traverseOutside(vit, next_cell, newIdx, sensor_idx);
                }
                else
                    current_cell->info().outside_score+=alpha_vis;
                break;
            }// end of intersection result
        }// end of finite cell check
        else{
            // TODO: what I am doing here is not correct, the if statement should be activated, but it crashes if it is
//            if(data_.Dt.tetrahedron(current_cell).has_on_positive_side(vit->info().sensor_positions[sensor_idx])){
                current_cell->info().outside_score+=alpha_vis;
//            }
        }
    }
}



void RayCaster::inside(Vertex_handle vit, int sensor_idx, int set_sink_weight){

    std::vector<Cell_handle> inc_cells;
    data_.Dt.incident_cells(vit, std::back_inserter(inc_cells));

    for(std::size_t i=0; i < inc_cells.size(); i++){
        Cell_handle current_cell = inc_cells[i];

        if(!data_.Dt.is_infinite(current_cell))
        {
            int cellBasedVertexIndex = current_cell->index(vit);

            Triangle tri = data_.Dt.triangle(current_cell, cellBasedVertexIndex);
            Point intersectionPoint;
            Point source = vit->point();
            Vector rayV = -Vector(vit->info().sensor_positions[sensor_idx] - source);

            // check if there is an intersection between the current ray and current triangle
            bool intersect = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            if(intersect){
                double dist2 = CGAL::squared_distance(intersectionPoint, source);
                current_cell->info().facet_weights[cellBasedVertexIndex] += ComputeDistanceProb(dist2);
                if(set_sink_weight)
                    current_cell->info().inside_score+=alpha_vis;
                break;
            }
        }// end of finite cell check
        else{
            // TODO: what I am doing here is not correct, the if statement should be activated, but it crashes if it is
//            if(data_.Dt.tetrahedron(current_cell).has_on_positive_side(vit->info().sensor_positions[sensor_idx])){
//                const int nvis=vit->info().sensor_positions.size();
//                const double alpha = ComputeVisibilityProb(nvis*nvis);
//                current_cell->info().outside_score+=alpha_vis;
//                double dist2 = CGAL::squared_distance(vit->info().sensor_positions[sensor_idx], vit->point());
//                current_cell->info().outside_score += ComputeDistanceProb(dist2);

//            }
        }
    }
} // end of going around all the cells of a vertex



void RayCaster::run(int set_sink_weight){

    string scoreType = options_.score_type;
    cout << "\nRay tracing..." << endl;
    cout << "\t-Trace " << options_.number_of_rays << " ray(s) to every point" << endl;
//    cout << "\t-Score type: " << options_.score_type << endl;

    auto start = chrono::high_resolution_clock::now();

    Delaunay::Finite_vertices_iterator vit;
    for(vit = data_.Dt.finite_vertices_begin() ; vit != data_.Dt.finite_vertices_end() ; vit++){
        int current_number_of_images = vit->info().sensor_positions.size();
        int number_of_images;
        if(options_.number_of_rays < 0 || options_.number_of_rays > current_number_of_images)
            number_of_images = current_number_of_images;
        else
            number_of_images = options_.number_of_rays;
        for(int sensor_idx = 0; sensor_idx < number_of_images; sensor_idx++){
            // construct ray from input point to sensor
            // collect outside votes
            outside(vit, sensor_idx);
            // collect inside votes
            inside(vit, sensor_idx, set_sink_weight);
        }
    }
    if(scoreType == "bodis"){
        double gamma = 2.0;
        Delaunay::All_cells_iterator aci;
        for(aci = data_.Dt.all_cells_begin(); aci != data_.Dt.all_cells_end(); aci++){
            aci->info().inside_score = 1 - (exp(-aci->info().inside_score/gamma));
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto full_duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done in " << full_duration.count() << "s" << endl;
}


// end of namespace rayTracing
}







