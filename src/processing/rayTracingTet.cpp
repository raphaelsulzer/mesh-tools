#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <processing/rayTracingTet.h>
#include <util/vectorArithmetic.h>

// TODO: check if I can replace the whole thing with the nearest_vertex(const Point& p, Cell_handle start) function,
// which can be found in the Delaunay_triangulation_3.h file in /usr/lib/CGAL

using namespace std;

namespace processing{

// TODO: finally save facet intersections instead of cells
// could also use the COLMAP code for this. make a new branch first
// and then show the stuff on Berger comparing the Labatu
// I need to show / see for myself what the LEARNING actually does in a CONTROLLED ENVIRONMENT
// even if OccNet is better here, they will not be better on large scale MVS!!

bool rayTriangleIntersection(const Point& rayOrigin,
                           const Vector& rayVector,
                           const Triangle& inTriangle,
                           Point& outIntersectionPoint){

    // implemented after: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    using namespace Eigen;
    const double EPSILON = 0.0000001;

    // init "Eigen" vectors
    Vector3d rayO(rayOrigin.x(), rayOrigin.y(), rayOrigin.z());
    Vector3d rayV(rayVector.x(), rayVector.y(), rayVector.z());

    Vector3d vertex0(inTriangle.vertex(0).x(), inTriangle.vertex(0).y(), inTriangle.vertex(0).z());
    Vector3d vertex1(inTriangle.vertex(1).x(), inTriangle.vertex(1).y(), inTriangle.vertex(1).z());
    Vector3d vertex2(inTriangle.vertex(2).x(), inTriangle.vertex(2).y(), inTriangle.vertex(2).z());

    Vector3d edge1, edge2, h, s, q;
    double a;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = rayV.cross(edge2);
    a = edge1.dot(h);
    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.
    double f = 1.0/a;
    s = rayO - vertex0;
    double u = f * s.dot(h);    // barycentric coordinate u
    if (u < 0.0 || u > 1.0)     // if u is not between 0 and 1, intersection point does not lie in triangle
        return false;
    q = s.cross(edge1);
    double v = f * rayV.dot(q); // barycentric coordinate v
    if (v < 0.0 || u + v > 1.0) // if v is not between 0 and 1, intersection point does not lie in triangle
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * edge2.dot(q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayOrigin + rayVector * t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}


////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
pair<double, double> wasureScore(double dist2, Vertex_handle vh, bool inside){

    double score_inside;
    double score_outside;
    // noise
    // TODO: if not noise estimation is done before, sigma_d should be set to 1
    double sigma_d = vh->info().sigma;
    // scene thickness
    // this parameter basically weighs the inside weight; thus for weakly supported surfaces it is good to keep this high
    double sigma_o = 1;
    // scale of the outside area?? // good for fontaine dataset is 1.0
    // outside weights
    double sigma_e = 1;
    // not to be confused with the following, which means if I am walking inside/outside
    if(inside){
        score_inside = (1 - 0.5*exp(-pow((sqrt(dist2)/sigma_d),2)))*exp(-pow((sqrt(dist2)/sigma_o),2));
        score_outside = 0.5*exp(-pow((sqrt(dist2)/sigma_d),2));
    }
    else{
        score_outside = (1 - 0.5*exp(-pow((sqrt(dist2)/sigma_d),2)))*exp(-pow((sqrt(dist2)/sigma_e),2));
        score_inside = 0.5*exp(-pow((sqrt(dist2)/sigma_d),2));
    }

    // first element is the outside score, second the inside score
    pair<double,double> scores(score_outside, score_inside);
    return scores;
}

pair<double, double> labatutScore(double dist2, Vertex_handle vh, bool inside){

    double score_inside;
    double score_outside;

    double sigma = vh->info().sigma;
    int alpha = vh->info().alpha;
//    int alpha = 1;
    if(inside){
        score_inside = alpha*(1-exp(-dist2/(2*sigma*sigma)));
        score_outside = 0;
    }
    else{
        score_outside = alpha*(1-exp(-dist2/(2*sigma*sigma)));;
        score_inside = 0;
    }

    // first element is the outside score, second the inside score
    pair<double,double> scores(score_outside, score_inside);
    return scores;
}


int traverseCells(Delaunay& Dt,
                  Cell_handle& current_cell, Ray ray, Vertex_handle vit, int oppositeVertex,
                  int image, bool inside,
                  string scoreType)
{
    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // source               The "Delaunay" point of the current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point

    // if(processed.find(current_cell) != processed.end()){
    //     return 0;
    // }

    if(!Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int cellBasedVertexIndex = (oppositeVertex+i)%4;
            Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
            Point intersectionPoint;
            Point rayO = ray.source();
            Vector rayV = ray.to_vector();
            bool result = rayTriangleIntersection(rayO, rayV, tri, intersectionPoint);
            // check if there is an intersection between the current ray and current triangle
            if(result){
                // get the distance between the source of the ray and the intersection point with the current cell
                // TODO:: maybe just pass this dist to the next traversel, and if the next dist is not bigger than break,
                // ...this could potentially be a much faster way than the processed-set
                double dist2 = CGAL::squared_distance(intersectionPoint, rayO);
                ////// weighting by how many cameras see a point
//                dist2*=vit->info().sensor_positions.size();


                ////// calculate the score for the current cell based on the distance //////
                std::pair<double, double> score;
                if(scoreType == "mine" || scoreType == "bodis")
                    score = labatutScore(dist2, vit, inside);
                else if(scoreType == "wasure")
                    score = wasureScore(dist2, vit, inside);
                else
                    cout << "not a valid score type" << endl;
                // add score to the cell for graph-cut
                current_cell->info().inside_score += score.second;
                current_cell->info().outside_score += score.first;
                // add features to cell for learning
                if(inside){
                    current_cell->info().inside_dist += sqrt(dist2);
                    current_cell->info().inside_count += 1;
                }
                else{
                    current_cell->info().outside_dist += sqrt(dist2);
                    current_cell->info().outside_count += 1;
                }
                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                Facet mirror_fac = Dt.mirror_facet(std::make_pair(current_cell, cellBasedVertexIndex));
                Cell_handle newCell = mirror_fac.first;
                int newIdx = mirror_fac.second;
                // the if stops the traversal once I hit the cell with the sensor center
                if(!Dt.tetrahedron(current_cell).has_on_positive_side(vit->info().sensor_positions[image])){
                    traverseCells(Dt,
                                  newCell, ray, vit, newIdx,
                                  image, inside,
                                  scoreType);
                }
                break;
            }
        }
    }
    // put outside score of infinite cell very high
    else{
        Point current_point = vit->point();
        Point current_sensor = vit->info().sensor_positions[image];

        Cell_handle locate_cell = Dt.locate(current_sensor);
        if(locate_cell != current_cell)
            return 0;

        double dist2 = sqrt(CGAL::squared_distance(current_sensor, current_point));
        std::pair<double, double> score;
        if(scoreType == "mine" || scoreType == "bodis")
            score = labatutScore(dist2, vit, inside);
        else if(scoreType == "wasure")
            score = wasureScore(dist2, vit, inside);
        else
            cout << "not a valid score type" << endl;
        // add score to the cell for graph-cut
        current_cell->info().inside_score += score.second;
        current_cell->info().outside_score += score.first;
        // add features to cell for learning
        if(inside){
            current_cell->info().inside_dist += sqrt(dist2);
            current_cell->info().inside_count += 1;
        }
        else{
            current_cell->info().outside_dist += sqrt(dist2);
            current_cell->info().outside_count += 1;
        }

//        // check if the current cell is the one the ray is coming from
//        Cell_handle locate_cell = Dt.locate(vit->info().sensor_positions[image]);
//        if(!inside && (locate_cell == current_cell)){
//            // TODO: put the actual score here, with dist = dist between sensor and point
//            // and actually I could also traverse the infinite cells until I am at the one with the sensor
//            current_cell->info().outside_score+=1;
//        }
////        if(inside)
////            current_cell->info().inside_score+=infiniteScore;
////        else
////            current_cell->info().outside_score+=infiniteScore;
    }
    return 0;
}

void firstCell(Delaunay& Dt,
               Delaunay::Finite_vertices_iterator& vit, Ray ray, int image,
               bool inside,
               string scoreType){

    // ray constructed from point origin to (end of) normal
    // introduces a ray r with source p and with a direction given by v.
//    Ray ray(vit->point(), vit->info().sensor_vec);

    // make the inside ray
    // in fact MicMac saves the camera normals pointing away from the camera,
    // so I take the opposite ray for outside traversal and the normal ray for inside
    if(inside)
        ray = ray.opposite();

    // vector of incident cells to the current vertex (from vertex iterator vit)
    // TODO: rewrite this to directly use the circulator as an iterator, should be easy!
    // and give speed up!
    std::vector<Cell_handle> inc_cells;
    Dt.incident_cells(vit, std::back_inserter(inc_cells));
    // 1. for every cell of incident cells, check if facet(cell, vertex) intersects with vertex normal
    // so this is checking in all directions of a vertex, but we will only have an intersection in (maximum) one direction
    // why only in one direction? because I'm only checking the OPPOSITE facade. It of course also intersects with the bordering facets
    // of all the neighbouring cells
    for(std::size_t i=0; i < inc_cells.size(); i++){
        Cell_handle current_cell = inc_cells[i];

        if(!Dt.is_infinite(current_cell))
        {
            int cellBasedVertexIndex = current_cell->index(vit);
//            Facet current_facet = std::make_pair(current_cell, cellBasedVertexIndex);
//            if(Dt.is_infinite(current_facet))
//                continue;
//            Triangle tri = Dt.triangle(current_facet);
            Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
            Point intersectionPoint;
            Point source = vit->point();
            Vector rayV = ray.to_vector();

            // check if there is an intersection between the current ray and current triangle
            bool intersect = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            if(intersect){
                double dist2;
//                if(inside){
//                    Sphere sph = Sphere(current_cell->vertex(0)->point(), current_cell->vertex(1)->point(),
//                                         current_cell->vertex(2)->point(), current_cell->vertex(3)->point());
//                    dist2 = CGAL::squared_distance(sph.center(), source);
//                }
//                else{
                    // get the distance between the source of the ray and the intersection point with the current cell
                    dist2 = CGAL::squared_distance(intersectionPoint, source);
//                }

                ////// weighting by how many cameras see a point
//                dist2*=vit->info().sensor_positions.size();
                ///////// calculate a weight of the current cell, based on the angles between sensor and triangles /////////
//                Point p1,p2,p3;
//                Vector fn;
//                double angle = 10;
//                for(int f = cellBasedVertexIndex+1; f < cellBasedVertexIndex+4; f++){
//                    // getting correctly oriented facet normals
//                    // implemented from here:
//                    // http://cgal-discuss.949826.n4.nabble.com/normal-vector-of-a-facet-td1580004.html#a1590886
//                    Facet fac(current_cell, f);
//                    // points on the facet
//                    p1 = fac.first->vertex((fac.second+1)&3)->point();
//                    p2 = fac.first->vertex((fac.second+2)&3)->point();
//                    p3 = fac.first->vertex((fac.second+3)&3)->point();

//                    fn = ( fac.second % 2 == 1) ?
//                                       CGAL::normal(p1, p2, p3):
//                                       CGAL::normal(p1, p3, p2);
//                    // still the angle can be between 0-90 or 180-270
//                    double temp_angle = std::acos((fn * rayV) / std::sqrt(fn.squared_length()*rayV.squared_length()));
//                    double new_angle;
//                    if(temp_angle > M_PI)
//                        new_angle = temp_angle - (M_PI+M_PI/2);
//                    else
//                        new_angle = -temp_angle + M_PI/2;
//                    angle = std::min(new_angle,angle);
//                }
//                dist2*=angle;

                ///////// calculate the score for the current cell based on the distance /////////
                std::pair<double, double> score;
                if(scoreType == "mine" || scoreType == "bodis")
                    score = labatutScore(dist2, vit, inside);
                else if(scoreType == "wasure")
                    score = wasureScore(dist2, vit, inside);
                else
                    cout << "not a valid score type" << endl;

                // for now only using volume weight on the secondary cells, since it has the tendency
                //  to oversmooth regions with little data support
//                double vol = Dt.tetrahedron(current_cell).volume();
                // add score to the cell for graph-cut
                current_cell->info().inside_score += score.second;
                current_cell->info().outside_score += score.first;
                // add features to cell for learning
                if(inside){
                    current_cell->info().inside_dist += sqrt(dist2);
                    current_cell->info().inside_count += 1;
                }
                else{
                    current_cell->info().outside_dist += sqrt(dist2);
                    current_cell->info().outside_count += 1;
                }

                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                // now from this new cell that I am in (get it from mirror_fac), iterate over all the triangles that are not the mirror triangle
                // and check if there is an intersection
                // this should be entered again at if(!Dt.is_infinite(current_cell)), since like this I can check if the cell is not already the infinite cell
                // so start from there to put this into a function
                if(!inside){
                    Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                    Facet mirror_fac = Dt.mirror_facet(fac);
                    Cell_handle newCell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    // go to next cell
                    // the if stops the traversal once I hit the cell with the sensor center
                    if(!Dt.tetrahedron(current_cell).has_on_positive_side(vit->info().sensor_positions[image])){
                        traverseCells(Dt,
                                      newCell, ray, vit, newIdx,
                                      image, inside,
                                      scoreType);
                    }
                    // TODO: maybe put an else here and increase the score, because this cell should DEFINITELY be outside
                }
                // if there was a match, break the loop around this vertex, so it can go to the next one
                // this is only done for speed reasons, it shouldn't have any influence on the label
                // because a ray can only hit more than one facet of a cell if it hits another point of the cell
                // in which case it goes THROUGH a facet of the cell, in which case the intersection is not a point
                // but an edge.
                // it does however make a difference if this is turn on or not. why??
                break;
            }// end of intersection result
        }// end of finite cell check
        else{         // put outside score of infinite cell very high
            // TODO: put the actual score here, with dist = dist between sensor and point
            // and actually I could also traverse the infinite cells until I am at the one with the sensor

            Point current_point = vit->point();
            Point current_sensor = vit->info().sensor_positions[image];
            Point current_point_extended = current_point + 0.0001*ray.to_vector();

            double dist2;

            if(inside){
                // what this does is saying if the ray crosses this current infinite cell
                // give it a score
                // implemented: if a point right behind/in front of the
                Cell_handle locate_cell = Dt.locate(current_point_extended);
                if(locate_cell != current_cell)
                    continue;
                // TODO: think about what distance to put here
                dist2 = 1;
            }
            else{
                Cell_handle locate_cell = Dt.locate(current_sensor);
                if(locate_cell != current_cell)
                    continue;
                dist2 = sqrt(CGAL::squared_distance(current_sensor, current_point));
            }

            std::pair<double, double> score;
            if(scoreType == "mine" || scoreType == "bodis")
                score = labatutScore(dist2, vit, inside);
            else if(scoreType == "wasure")
                score = wasureScore(dist2, vit, inside);
            else
                cout << "not a valid score type" << endl;
            // add score to the cell for graph-cut
            current_cell->info().inside_score += score.second;
            current_cell->info().outside_score += score.first;
            // add features to cell for learning
            if(inside){
                current_cell->info().inside_dist += sqrt(dist2);
                current_cell->info().inside_count += 1;
            }
            else{
                current_cell->info().outside_dist += sqrt(dist2);
                current_cell->info().outside_count += 1;
            }
            // TODO: what weight to give if I am traversing to inside?

            // this is basically useless, because it will just give same inside and outside score.
//        if(inside)
//                current_cell->info().inside_score+=infiniteScore;
//            else
//                current_cell->info().outside_score+=infiniteScore;
        }
    } // end of going around all the cells of a vertex
} // end of first cell


// if it intersect the cell, add a pair with vertexIndex+nCells and dist to intersection to the cell


void rayTracing(Delaunay& Dt, runningOptions ro){

    int max_number_of_images = ro.number_of_rays;
    string scoreType = ro.score_type;
    cout << "\nRay tracing..." << endl;
    cout << "\t-Trace " << max_number_of_images << " ray(s) to every point" << endl;
    cout << "\t-Score type: " << scoreType << endl;

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
            if(!ro.tetTracing)
                firstCell(Dt, vit, ray, r, 0, scoreType);
            // collect inside votes
            firstCell(Dt, vit, ray, r, 1, scoreType);
        }
    }

    if(scoreType == "bodis"){
        double gamma = 2.0;
        Delaunay::All_cells_iterator aci;

        for(aci = Dt.all_cells_begin(); aci != Dt.all_cells_end(); aci++){
            aci->info().inside_score = 1 - (exp(-aci->info().inside_score/gamma));
        }

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "\t-done in " << full_duration.count() << "s" << std::endl;
}


// end of namespace rayTracing
}






