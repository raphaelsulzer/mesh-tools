#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <labatut/rayTracingFacet.h>
#include <util/geometricOperations.h>

using namespace std;


///////// Ray tracing implemented according to "Robust and efficient surface reconstruction from range data" by Labatut et al. 2009:
/// Line-of-sight (outside) and ray (inside) information is accumulated on the oriented facets of a 3D Delaunay triangulation.
/// Relevant parameters are
/// alpha, which controls the weight per line-of-sight / ray
/// sigma, which controls the softening of the weight around the input point. Trades accuracy vs complexity of the surface.
/// tau, which controls the maximal inside traversal.
/// In Labatut et al. tau = sigma, however, I find that this removes the possibility of controlling the complexity of the surface with
/// sigma, as with higher sigma, the inside traversal goes too far, and shrinks the surface.


namespace processing{


void RayCaster::traverseOutside(Vertex_handle vit, Cell_handle& current_cell,  int oppositeVertex, int sensor_idx, double tdist2 =0){

    if(!data_.Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int cellBasedVertexIndex = (oppositeVertex+i)%4;
            Triangle tri = data_.Dt.triangle(current_cell, cellBasedVertexIndex);
            Point intersectionPoint;
            Point source = vit->point();
            Vector rayV = Vector(vit->info().sensor_positions[sensor_idx] - source);
//            Tetrahedron tet;
            bool intersection = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            if(intersection){
                // the if stops the traversal once I hit the cell with the sensor center
//                const Tetrahedron tet = data_.Dt.tetrahedron(current_cell);
//                cout<<CGAL::volume(tet.vertex(0),tet.vertex(1),tet.vertex(2),tet.vertex(3))<<endl;
//                if(tet.is_degenerate())
//                    continue;
//                cout << current_cell->info().global_idx << endl;
                double dist2 = CGAL::squared_distance(intersectionPoint, source);
                if(dist2 <= tdist2)
                    break;
                tdist2=dist2;
                if(!data_.Dt.tetrahedron(current_cell).has_on_bounded_side(vit->info().sensor_positions[sensor_idx])){
                    // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                    Facet mirror_fac = data_.Dt.mirror_facet(std::make_pair(current_cell, cellBasedVertexIndex));
                    Cell_handle next_cell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    current_cell->info().facet_weights[cellBasedVertexIndex] += ComputeDistanceProb(dist2);
                    traverseOutside(vit, next_cell, newIdx, sensor_idx, tdist2);
                }
                else
                    //Labatut et al:
                    // "Note that the only tetrahedra to get a non zero-weighted link to
                    // the source are those (possibly infinite ones) containing the laser
                    // sources or sensors optical centers."
                    current_cell->info().outside_score+=alpha_vis;
                break;
            }
        }
    }
    else{
        // This means infinite cell on the outside get source weight
        // this one is set no matter if the shapes is closed or open.
        // the reason is that this cell is entered through a line-of-sight, thus it has strong outside evidence.
        // It departs from Labatut et al. (see citation above), but if it is not set, closed objects, for which the sensor
        // is not inside any cell do not get any source weight and the algorithm doesn't reconstruct anything.
        // One way to circumvent it would be to add the sensor positions to the triangulation.
        if(options_.closed_prior)
            current_cell->info().outside_score+=alpha_vis;
    }
    return;
}


void RayCaster::traverseInside(Vertex_handle vit, Cell_handle& current_cell,  int oppositeVertex, int sensor_idx){

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
                if(!data_.Dt.tetrahedron(current_cell).has_on_bounded_side(vit->info().sensor_positions[sensor_idx])){
                    double dist2 = CGAL::squared_distance(intersectionPoint, source);
                    // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                    Facet mirror_fac = data_.Dt.mirror_facet(std::make_pair(current_cell, cellBasedVertexIndex));
                    Cell_handle next_cell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    current_cell->info().facet_weights[cellBasedVertexIndex] += ComputeDistanceProb(dist2);
                    if(sqrt(dist2)>=options_.labatut_tau)
                        next_cell->info().inside_score+=alpha_vis;
                    else
                        traverseInside(vit, next_cell, newIdx, sensor_idx);
                }
                else
                    //Labatut et al:
                    // "Note that the only tetrahedra to get a non zero-weighted link to
                    // the source are those (possibly infinite ones) containing the laser
                    // sources or sensors optical centers."
                    cout << "cell with sensor in traverse" << endl;
                    current_cell->info().outside_score+=alpha_vis;
                break;
            }
        }
    }
    else{
        // this means infinite cell on the inside get sink weight
        // this one is set no matter if the shapes is closed or open.
        // the reason is that this cell is entered through a ray, thus it has strong inside evidence.
//        current_cell->info().inside_score+=alpha_vis;
    }
    return;
}


void RayCaster::outside(Vertex_handle vit, int sensor_idx){

    std::vector<Cell_handle> inc_cells;
    data_.Dt.incident_cells(vit, std::back_inserter(inc_cells));
    double tdist2=0;
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
                double dist2 = CGAL::squared_distance(intersectionPoint, source);
                if(dist2 <= tdist2)
                    break;
                tdist2=dist2;
                if(!data_.Dt.tetrahedron(current_cell).has_on_bounded_side(vit->info().sensor_positions[sensor_idx])){
                    double dist2 = CGAL::squared_distance(intersectionPoint, source);
                    Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                    Facet mirror_fac = data_.Dt.mirror_facet(fac);
                    Cell_handle next_cell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    // set the facet weight
                    current_cell->info().facet_weights[cellBasedVertexIndex] += ComputeDistanceProb(dist2);
                    traverseOutside(vit, next_cell, newIdx, sensor_idx);
                }
                else
                    //Labatut et al:
                    // "Note that the only tetrahedra to get a non zero-weighted link to
                    // the source are those (possibly infinite ones) containing the laser
                    // sources or sensors optical centers."
                    cout << "cell with sensor in outside" << endl;
                    current_cell->info().outside_score+=alpha_vis;
                break;
            }// end of intersection result
        }// end of finite cell check
    }
}



void RayCaster::inside(Vertex_handle vit, int sensor_idx){

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
                Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                Facet mirror_fac = data_.Dt.mirror_facet(fac);
                Cell_handle next_cell = mirror_fac.first;
                int newIdx = mirror_fac.second;
                // set the facet weight
                current_cell->info().facet_weights[cellBasedVertexIndex] += ComputeDistanceProb(dist2);
                if(sqrt(dist2)>=options_.labatut_tau)
                    next_cell->info().inside_score+=alpha_vis;
                else
                    traverseInside(vit, next_cell, newIdx, sensor_idx);
                break;
            }
        }// end of finite cell check
    }
} // end of going around all the cells of a vertex



void RayCaster::run(){

//    int i = 0;
//    int q;
//    for(auto vit = data_.Dt.finite_vertices_begin() ; vit != data_.Dt.finite_vertices_end() ; vit++){
//        cout << "Point " << i++ << " : " << vit->point() << endl;
//        cout << "Sensor " << i << " : " << vit->info().sensor_positions[0] << endl;
//        cout << "Vector " << i << " : " << vit->info().sensor_vec << endl;
//        q=5;
//    }

    cout << "\nRay tracing..." << endl;
    cout << "\t-Trace " << options_.number_of_rays << " ray(s) to every point" << endl;
    cout << "\t-Labatut's \u03B1_vis set to " << options_.labatut_alpha << endl;
    cout << "\t-Labatut's \u03C3 set to " << options_.labatut_sigma << endl;
    cout << "\t-Max. inside traversal (beyond first cell) \u03C4 set to " << options_.labatut_tau << endl;
    cout << "\t-Infinite cells forced outside? " << options_.closed_prior << endl;

    auto start = chrono::high_resolution_clock::now();

    Delaunay::Finite_vertices_iterator vit;
    for(vit = data_.Dt.finite_vertices_begin() ; vit != data_.Dt.finite_vertices_end() ; vit++){
        int current_number_of_images = vit->info().sensor_positions.size();
        if(current_number_of_images == 0)
            cout << "\nWARNING: Point with no sensor!" << endl;
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
            inside(vit, sensor_idx);
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto full_duration = chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done in " << full_duration.count() << "s" << endl;
}


// end of namespace rayTracing
}







