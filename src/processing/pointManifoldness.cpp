//#include <base/cgal_typedefs.h>
//#include <IO/fileIO.h>
//#include <processing/rayTracingTet.h>
//#include <processing/tetIntersection.h>
//#include <processing/edgeManifoldness.h>
//#include <unordered_set>


//////////////////////////////////////
////////// Point manifoldness ////////
//////////////////////////////////////


//// TODO: could also implement this with kind of a recursive function that always looks in neighbourhood
//// so start with a seed cell, go to next cell and see if it has a different label and if it is connected to the previous cell
//// ... if I have more than two regions, point is not manifold

//// get non manifold points
//int isManifoldPoint(Delaunay& Dt, Vertex_handle p){

//    std::vector<Cell_handle> incident_cells;
//    Dt.incident_cells(p, std::back_inserter(incident_cells));
//    std::vector<Cell_handle>::iterator cc;
//    // sort cells in inside and outside cells
//    std::vector<Cell_handle> inside_cells;
//    std::vector<Cell_handle> outside_cells;
//    for(cc = incident_cells.begin(); cc != incident_cells.end(); cc++){
//        if(getCellLabel(*cc) == 0)
//            inside_cells.push_back(*cc);
//        else
//            outside_cells.push_back(*cc);
//    }
//    // if one of them is empty, point is manifold
//    int is = inside_cells.size();
//    int os = outside_cells.size();
//    std::vector<Cell_handle> region;
//    std::vector<Cell_handle> rest;
//    std::vector<Cell_handle> stack;
//    if(is < 2 || os < 2){
//        return 1;
//    }
//    else if(is <= os){
//        region.push_back(inside_cells[0]);
//        for(int i = 1; i < is; i++)
//            rest.push_back(inside_cells[i]);
//    }
//    else if(is > os){
//        region.push_back(outside_cells[0]);
//        for(int i = 1; i < os; i++)
//            rest.push_back(outside_cells[i]);
//    }
//    int new_stack = 1000;
//    int old_stack = 1001;
//    do{
//        Cell_handle first_rest_cell = rest[0];
//        bool first_rest_cell_added = 0;
//        for(auto r = region.begin(); r != region.end(); r++){
//            if(first_rest_cell->has_neighbor(*r)){
//                region.push_back(first_rest_cell);
//                rest.erase(rest.begin());
//                first_rest_cell_added = 1;
//                break;
//            }
//        }
//        if(!first_rest_cell_added){
//            stack.push_back(first_rest_cell);
//            rest.erase(rest.begin());
//        }
//        if(rest.size() == 0){
//            if(stack.size()==0)
//                return 1;
//            rest = stack;
//            old_stack=new_stack;
//            new_stack=stack.size();
//            stack.clear();
//        }
//    }while(old_stack > new_stack);
//    if(rest.size()==0 && stack.size()==0)
//        return 1;
//    else
//        return 0;
//}

//// get non manifold points
//int isManifoldPoint(Delaunay& Dt, Delaunay::Finite_vertices_iterator& p, std::vector<Cell_handle>& incident_cells){

//    Dt.incident_cells(p, std::back_inserter(incident_cells));
//    std::vector<Cell_handle>::iterator cc;
//    // sort cells in inside and outside cells
//    std::vector<Cell_handle> inside_cells;
//    std::vector<Cell_handle> outside_cells;
//    for(cc = incident_cells.begin(); cc != incident_cells.end(); cc++){
//        if(getCellLabel(*cc) == 0)
//            inside_cells.push_back(*cc);
//        else
//            outside_cells.push_back(*cc);
//    }
//    // if one of them is empty, point is manifold
//    int is = inside_cells.size();
//    int os = outside_cells.size();
//    std::vector<Cell_handle> region;
//    std::vector<Cell_handle> rest;
//    std::vector<Cell_handle> stack;
//    if(is < 2 || os < 2){
//        return 1;
//    }
//    else if(is <= os){
//        region.push_back(inside_cells[0]);
//        for(int i = 1; i < is; i++)
//            rest.push_back(inside_cells[i]);
//    }
//    else if(is > os){
//        region.push_back(outside_cells[0]);
//        for(int i = 1; i < os; i++)
//            rest.push_back(outside_cells[i]);
//    }
//    int new_stack = 1000;
//    int old_stack = 1001;
//    do{
//        Cell_handle first_rest_cell = rest[0];
//        bool first_rest_cell_added = 0;
//        for(auto r = region.begin(); r != region.end(); r++){
//            if(first_rest_cell->has_neighbor(*r)){
//                region.push_back(first_rest_cell);
//                rest.erase(rest.begin());
//                first_rest_cell_added = 1;
//                break;
//            }
//        }
//        if(!first_rest_cell_added){
//            stack.push_back(first_rest_cell);
//            rest.erase(rest.begin());
//        }
//        if(rest.size() == 0){
//            if(stack.size()==0)
//                return 1;
//            rest = stack;
//            old_stack=new_stack;
//            new_stack=stack.size();
//            stack.clear();
//        }
//    }while(old_stack > new_stack);
//    if(rest.size()==0 && stack.size()==0)
//        return 1;
//    else
//        return 0;
//}



//// get non manifold edges
//int isPointManifoldClique(Delaunay& Dt, Delaunay::Finite_vertices_iterator& p){


//    // first check if input point is manifold
//    if(!isManifoldPoint(Dt, p))
//        return 0;

//    // now assemble all other points in a set
//    std::vector<Cell_handle> incident_cells;
//    Dt.incident_cells(p, std::back_inserter(incident_cells));
//    std::vector<Cell_handle>::iterator cc;
//    std::unordered_set<Vertex_handle> vertices;
//    for(cc = incident_cells.begin(); cc != incident_cells.end(); cc++){

//        // add its 4th face neighbour
//        int idx = (*cc)->index(p);

//        // switch case
//        // if idx 0, other vertices to add are 1,2,3, and so on
//        std::vector<int> other_idx(3);
//        switch(idx){
//        case 0: other_idx[0] = 1; other_idx[1] = 2; other_idx[2] = 3;
//            break;
//        case 1: other_idx[0] = 0; other_idx[1] = 2; other_idx[2] = 3;
//            break;
//        case 2: other_idx[0] = 0; other_idx[1] = 1; other_idx[2] = 3;
//            break;
//        case 3: other_idx[0] = 0; other_idx[1] = 1; other_idx[2] = 2;
//        }
//        for(int i = 0; i < 3; i++)
//            vertices.insert((*cc)->vertex(other_idx[i]));
//    }
//    // now check all of them
//    std::unordered_set<Vertex_handle>::iterator vs;
//    for (vs = vertices.begin(); vs != vertices.end(); vs++){
//        if(!isManifoldPoint(Dt, *vs))
//            return 0;
//    }
//    return 1;
//}

//int isManifoldClique(Delaunay& Dt, Delaunay::Finite_edges_iterator& e){

//    // check for edge manifoldness
//    if(!isManifoldEdge(Dt,e))
//        return 0;

//    Cell_handle current_cell = e->first;

//    Vertex_handle v1 = current_cell->vertex(e->second);
//    Vertex_handle v2 = current_cell->vertex(e->third);
//    if(!isManifoldPoint(Dt, v1) && !isManifoldPoint(Dt, v2))
//        return 0;

//    Delaunay::Cell_circulator cc = Dt.incident_cells(*e, e->first);
//    Cell_handle first_cell = cc;
//    Edge current_edge = *e;
//    do{

//        int i = current_edge.second;
//        int j = current_edge.third;
//        int k = Dt.next_around_edge(i,j);

//        // check point manifoldness
//        Vertex_handle v3 = cc->vertex(k);
//        if(!isManifoldPoint(Dt, v3))
//            return 0;

//        // check first edge
//        Edge e1 = CGAL::Triple<Cell_handle,int,int>(cc,i,k);
//        if(!isManifoldEdge(Dt,e1))
//            return 0;
//        // check second edge
//        Edge e2 = CGAL::Triple<Cell_handle,int,int>(cc,j,k);
//        if(!isManifoldEdge(Dt,e2))
//            return 0;
//        // check third edge
//        Edge e3 = CGAL::Triple<Cell_handle,int,int>(cc,k,6-(i+j+k));
//        if(!isManifoldEdge(Dt,e3))
//            return 0;

//        // go to next cell with its corresponding edge
//        cc++;
//        int new_i;
//        switch(i){
//        case 0: new_i = 1;
//            break;
//        case 1: new_i = 0;
//            break;
//        case 2: new_i = 2;
//            break;
//        case 3: new_i = 3;
//        }
//        int new_j;
//        switch(j){
//        case 0: new_j = 1;
//            break;
//        case 1: new_j = 0;
//            break;
//        case 2: new_j = 2;
//            break;
//        case 3: new_j = 3;
//        }
//        current_edge = CGAL::Triple<Cell_handle, int, int>(cc, new_j, new_i);
//    }
//    while(first_cell != cc);

//    return 1;

//}


//double nonManifoldCliqueEnergy(const Delaunay& Dt, std::vector<Cell_handle> all_cells, double reg_weight){

//    std::vector<Cell_handle>::iterator cit;
//    double unary = 0;
//    double binary = 0;
//    for(cit = all_cells.begin(); cit != all_cells.end(); cit++){

//        // first cell
//        int cc_label = getCellLabel(*cit);

//        // add the unary term of this cell
//        if(cc_label==1){unary+=(*cit)->info().inside_score;} // if label = 1 = outside, penalize with inside score
//        else{unary+=(*cit)->info().outside_score;} // if label = 0 = inside, penalize with outside score

//        int current_index = (*cit)->info().global_idx;
//        for(int i = 0; i < 4; i++){

//            Cell_handle neighbour_cell = (*cit)->neighbor(i);
//            int neighbour_index = neighbour_cell->info().global_idx;

//            // if current_cell and neighbour_cell are BOTH infinite, then continue
//            if(Dt.is_infinite(*cit) && Dt.is_infinite(neighbour_cell)){
//                continue;
//            }

//            // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
//            if(neighbour_index < current_index)
//                continue;

//            // since i is giving me the cell that is opposite of vertex i, as well as the facet that is opposite of vertex i, I can just use that same index
//            Triangle tri = Dt.triangle(*cit, i);
//            binary+=sqrt(tri.squared_area());

//        }

//    }
//    return unary+reg_weight*binary;
//}


//void getNonManifoldPoints(Delaunay& Dt, double regularization_weight, std::string path){

//    auto start = std::chrono::high_resolution_clock::now();

//    path+="_nmPoints.ply";

//    std::fstream fo;
//    fo.open(path, std::fstream::out);

//    int nv = Dt.number_of_vertices();

//    printPLYHeader(fo,
//                   nv, 0,
//                   false, true, false, false, false, false, 15);

//    // iterate over all edges of the Dt
//    Delaunay::Finite_vertices_iterator fvi;
//    int count = 0;
//    for(fvi = Dt.finite_vertices_begin(); fvi != Dt.finite_vertices_end(); fvi++){

//        if(!isManifoldPoint(Dt, fvi)){
//            fo << fvi->point() << " " << "0 255 0" << std::endl;
//            std::array<unsigned char, 3> col = {0,255,0};
//            fvi->info().color = col;
//            count++;
//        }
//        else{
//            fo << fvi->point() << " " << "0 0 0" << std::endl;
//        }
//    }
//    std::cout << "number of non-manifold points: " << count << std::endl;
//    fo.close();

//}


//void fixNonManifoldPoints(Delaunay& Dt, double regularization_weight){

//    auto start = std::chrono::high_resolution_clock::now();

//    // iterate over all points of the Dt
//    Delaunay::Finite_vertices_iterator fvi;
//    for(fvi = Dt.finite_vertices_begin(); fvi != Dt.finite_vertices_end(); fvi++){

//        std::vector<Cell_handle> incident_cells;
//        if(!isManifoldPoint(Dt, fvi, incident_cells))
//            continue;


//        // make the bitset with according size
//        std::vector<Combination_score> combinations;
//        int number_of_cells = incident_cells.size();
//        // NOTE: problem is that there are two many combinations (e.g. 2^30 possible) around a point
//        int number_of_possible_combinations = makeCombinations(number_of_cells, combinations);

//        // iterate over the container, relabel the cells with the bitset, and get their corresponding energy and save it in a vector
//        for(int c = 0; c < number_of_possible_combinations; c++){
//            for(int v = 0; v < number_of_cells; v++)
//                incident_cells[v]->info().manifold_label = combinations[c].second[v];
//            // calculate the energy of the current configuration
//            combinations[c].first = nonManifoldCliqueEnergy(Dt, incident_cells, regularization_weight);
//        }

//        // sort the container while keeping track of its original index, which gives you the corresponding configuration from the bitset index
//        std::sort(combinations.begin(), combinations.end(), sortCombinations);

//        // take the lowest energy and check if it is manifold
//        int solution_found = 0;
//        for(int sc = 0; sc < number_of_possible_combinations; sc++){

//            // relabel to the correct combination
//            for(int v = 0; v < number_of_cells; v++){
//                incident_cells[v]->info().manifold_label = combinations[sc].second[v];
//            }
//            if(isPointManifoldClique(Dt, fvi)){
//                for(int v = 0; v < number_of_cells; v++){
//                    incident_cells[v]->info().gc_label = combinations[sc].second[v];
//                }
//                solution_found = 1;
//                break;
//            }
//        } // end of trying all possible solutions
//        if(!solution_found){
//            std::cout << "no manifold combination found for edge "
//                      << fvi->info().global_idx
//                      << std::endl;
//        }
//    }

//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    std::cout << "Fixed non-manifold points in " << duration.count() << "s" << std::endl;
//}
