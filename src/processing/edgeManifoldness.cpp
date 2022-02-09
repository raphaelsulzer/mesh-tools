#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <processing/rayTracingTet.h>
#include <processing/tetIntersection.h>
#include <processing/edgeManifoldness.h>
#include <unordered_set>

// idea: take a clique, starting with a non manifold edge, ALL its sourounding facets and their facet neighbours.
// now optimize a binary linear program by enforcing that an edge of this clique either has 2 or 0 facet neighbours.
// the idea is that the neighbouring facets of the neighbours of the non manifold edge, which a lot of them did not get selected, now
// contribute with their data attachment (low difference between inside and outside label -> take the inverse)
// to the non selection of the obvious case where the
// 2 facet around the non-manifold edge with the highest difference of inside outside get selected.


// fix the region by applying face selection in the way that
// PolyFit does it. Use their Gurobi energy minimization solver for it (their energy cannot be minimized with a graph-cut).
// Problem will be that it might be to slow, and that it is not clear how to extract region with problem.
// problem: need to incrementally build the polygon surfaces for that to get a "controlled" facet ID with which I can query the neighbouring
// cells to get their label

////////////////////////////////////
///////// Edge manifoldness ////////
////////////////////////////////////

int isManifoldEdge(Delaunay& Dt, Delaunay::Finite_edges_iterator& e){

    Delaunay::Cell_circulator cc = Dt.incident_cells(*e);
    Cell_handle first_cell = cc;
    int borders = 0;
    do{
        // first cell
        int cc_label;
        if(cc->info().manifold_label < 2)
            cc_label = cc->info().manifold_label;
        else
            cc_label = cc->info().gc_label;
        // go to next cell
        int nc_label;
        cc++;
        if(cc->info().manifold_label < 2)
            nc_label = cc->info().manifold_label;
        else
            nc_label = cc->info().gc_label;
        // check if there is a border
        if(cc_label != nc_label){borders+=1;}
    }
    while(first_cell != cc && borders < 3);

    if(borders > 2)
        return 0;
    else
        return 1;
}

// get non manifold edges
int isManifoldEdge(Delaunay& Dt, Edge e){

    Delaunay::Cell_circulator cc = Dt.incident_cells(e);
    Cell_handle first_cell = cc;
    int borders = 0;
    do{
        // first cell
        int cc_label;
        if(cc->info().manifold_label < 2)
            cc_label = cc->info().manifold_label;
        else
            cc_label = cc->info().gc_label;
        // go to next cell
        int nc_label;
        cc++;
        if(cc->info().manifold_label < 2)
            nc_label = cc->info().manifold_label;
        else
            nc_label = cc->info().gc_label;
        // check if there is a border
        if(cc_label != nc_label){borders+=1;}
    }
    while(first_cell != cc && borders < 3);

    if(borders > 2)
        return 0;
    else
        return 1;
}


int isEdgeManifoldClique(Delaunay& Dt, Delaunay::Finite_edges_iterator& e){

    if(!isManifoldEdge(Dt,e))
        return 0;

    Delaunay::Cell_circulator cc = Dt.incident_cells(*e, e->first);
    Cell_handle first_cell = cc;
    Edge current_edge = *e;
    do{

        int i = current_edge.second;
        int j = current_edge.third;
        int k = Dt.next_around_edge(i,j);

        // check first edge
        Edge e1 = CGAL::Triple<Cell_handle,int,int>(cc,i,k);
        if(!isManifoldEdge(Dt,e1))
            return 0;
        // check second edge
        Edge e2 = CGAL::Triple<Cell_handle,int,int>(cc,j,k);
        if(!isManifoldEdge(Dt,e2))
            return 0;
        // check third edge
        Edge e3 = CGAL::Triple<Cell_handle,int,int>(cc,k,6-(i+j+k));
        if(!isManifoldEdge(Dt,e3))
            return 0;

        // go to next cell with its corresponding edge
        cc++;
        int new_i;
        switch(i){
        case 0: new_i = 1;
            break;
        case 1: new_i = 0;
            break;
        case 2: new_i = 2;
            break;
        case 3: new_i = 3;
        }
        int new_j;
        switch(j){
        case 0: new_j = 1;
            break;
        case 1: new_j = 0;
            break;
        case 2: new_j = 2;
            break;
        case 3: new_j = 3;
        }
        current_edge = CGAL::Triple<Cell_handle, int, int>(cc, new_j, new_i);
    }
    while(first_cell != cc);
    return 1;
}

int isEdgeManifoldClique(Delaunay& Dt, Edge edge, std::vector<Edge>& nm_edges){

    if(!isManifoldEdge(Dt,edge))
        return 0;

    Delaunay::Cell_circulator cc = Dt.incident_cells(edge, edge.first);
    Cell_handle first_cell = cc;
    do{

        int i = edge.second;
        int j = edge.third;
        int k = Dt.next_around_edge(i,j);

        // check first edge
        Edge e1 = CGAL::Triple<Cell_handle,int,int>(cc,i,k);
        if(std::find(nm_edges.begin(), nm_edges.end(), e1) == nm_edges.end()){
            if(!isManifoldEdge(Dt,e1))
                return 0;
        }
        // check second edge
        Edge e2 = CGAL::Triple<Cell_handle,int,int>(cc,j,k);
        if(std::find(nm_edges.begin(), nm_edges.end(), e2) == nm_edges.end()){
            if(!isManifoldEdge(Dt,e2))
                return 0;
        }
        // check third edge
        Edge e3 = CGAL::Triple<Cell_handle,int,int>(cc,k,6-(i+j+k));
        if(std::find(nm_edges.begin(), nm_edges.end(), e3) == nm_edges.end()){
            if(!isManifoldEdge(Dt,e3))
                return 0;
        }

        // go to next cell with its corresponding edge
        cc++;
        int new_i;
        switch(i){
        case 0: new_i = 1;
            break;
        case 1: new_i = 0;
            break;
        case 2: new_i = 2;
            break;
        case 3: new_i = 3;
        }
        int new_j;
        switch(j){
        case 0: new_j = 1;
            break;
        case 1: new_j = 0;
            break;
        case 2: new_j = 2;
            break;
        case 3: new_j = 3;
        }
        edge = CGAL::Triple<Cell_handle, int, int>(cc, new_j, new_i);
    }
    while(first_cell != cc);
    return 1;
}


int getCellLabel(Cell_handle& c, bool manifold=true){

    if(c->info().manifold_label < 2 && manifold)
        return c->info().manifold_label;
    else
        return c->info().gc_label;
}

double nonManifoldCliqueEnergy(const Delaunay& Dt, Edge e, double reg_weight){

    Delaunay::Cell_circulator cc = Dt.incident_cells(e, e.first);
    Cell_handle first_cell = cc;
    double unary = 0;
    double binary = 0;
    do{
        // first cell
        Cell_handle current_cell = cc;
        int cc_label = getCellLabel(current_cell);

        // add the unary term of this cell
        if(cc_label==1){unary+=cc->info().inside_score;} // if label = 1 = outside, penalize with inside score
        else{unary+=cc->info().outside_score;} // if label = 0 = inside, penalize with outside score

        // add the binary term of the other two facets that are not connected to this edge
        Facet f1 = std::make_pair(cc,e.second);
        Cell_handle mc1 = Dt.mirror_facet(f1).first;
        int mc1_label = getCellLabel(mc1);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(cc) && Dt.is_infinite(mc1)) && mc1_label != cc_label){
            Triangle tri = Dt.triangle(f1);
            binary+=sqrt(tri.squared_area());
        }
        Facet f2 = std::make_pair(cc,e.third);
        Cell_handle mc2 = Dt.mirror_facet(f2).first;
        int mc2_label = getCellLabel(mc2);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(cc) && Dt.is_infinite(mc2)) && mc2_label != cc_label){
            Triangle tri = Dt.triangle(f2);
            binary+=sqrt(tri.squared_area());
        }

        // go to next cell and add the binary term between previous cell and next cell
        cc++;
        int nidx = current_cell->index(cc); // the index of the new cell seen from the old cell = facet between previous and new cell
        Facet f3 = std::make_pair(current_cell,nidx);
        int mc3_label = getCellLabel(current_cell);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(current_cell) && Dt.is_infinite(cc)) && mc3_label != cc_label){
            Triangle tri = Dt.triangle(f3);
            binary+=sqrt(tri.squared_area());
        }
    }
    while(first_cell != cc);
    return unary+reg_weight*binary;
}

double nonManifoldCliqueEnergy(const Delaunay& Dt, Delaunay::Finite_edges_iterator& e, double reg_weight){

    Delaunay::Cell_circulator cc = Dt.incident_cells(*e, e->first);
    Cell_handle first_cell = cc;
    double unary = 0;
    double binary = 0;
    do{
        // first cell
        Cell_handle current_cell = cc;
        int cc_label = getCellLabel(current_cell);

        // add the unary term of this cell
        if(cc_label==1){unary+=cc->info().inside_score;} // if label = 1 = outside, penalize with inside score
        else{unary+=cc->info().outside_score;} // if label = 0 = inside, penalize with outside score

        // add the binary term of the other two facets that are not connected to this edge
        Facet f1 = std::make_pair(cc,e->second);
        Cell_handle mc1 = Dt.mirror_facet(f1).first;
        int mc1_label = getCellLabel(mc1);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(cc) && Dt.is_infinite(mc1)) && mc1_label != cc_label){
            Triangle tri = Dt.triangle(f1);
            binary+=sqrt(tri.squared_area());
        }
        Facet f2 = std::make_pair(cc,e->third);
        Cell_handle mc2 = Dt.mirror_facet(f2).first;
        int mc2_label = getCellLabel(mc2);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(cc) && Dt.is_infinite(mc2)) && mc2_label != cc_label){
            Triangle tri = Dt.triangle(f2);
            binary+=sqrt(tri.squared_area());
        }

        // go to next cell and add the binary term between previous cell and next cell
        cc++;
        int nidx = current_cell->index(cc); // the index of the new cell seen from the old cell = facet between previous and new cell
        Facet f3 = std::make_pair(current_cell,nidx);
        int mc3_label = getCellLabel(current_cell);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(current_cell) && Dt.is_infinite(cc)) && mc3_label != cc_label){
            Triangle tri = Dt.triangle(f3);
            binary+=sqrt(tri.squared_area());
        }
    }
    while(first_cell != cc);
    return unary+reg_weight*binary;
}


int makeCombinations(int n, std::vector<Combination_score>& c){
    const int x = pow(2,n);
    for(int i = 0; i < x; i++){
        const boost::dynamic_bitset<> b(n, i);
        c.push_back(std::make_pair(0.0,b));
    }
    return x;
}

bool sortCombinations(const Combination_score &a,
              const Combination_score &b)
{
    return (a.first < b.first);
}

void fixNonManifoldEdges(Delaunay& Dt, double regularization_weight){

    auto start = std::chrono::high_resolution_clock::now();


    // first collect all non-manifold cells
    Delaunay::Finite_edges_iterator fei;
    std::vector<Edge> nm_edges;
    for(fei = Dt.finite_edges_begin(); fei != Dt.finite_edges_end(); fei++){
        if(!isManifoldEdge(Dt, fei))
            nm_edges.push_back(*fei);
    }

    std::cout << "\nFound " << nm_edges.size() << " non-manifold edges. Trying to fix them..." << std::endl;

    int fixed_count=0;
    for(const auto& edge : nm_edges){
        // make a container with all incident cells
        std::vector<Cell_handle> cells_around_nmedge;
        Delaunay::Cell_circulator cc = Dt.incident_cells(edge);
        Cell_handle first_cell = cc;
        do{cells_around_nmedge.push_back(cc++);}
        while(first_cell != cc);

        // make the bitset with according size
        std::vector<Combination_score> combinations;
        int number_of_cells = cells_around_nmedge.size();
        int number_of_possible_combinations = makeCombinations(number_of_cells, combinations);

//        std::cout << "number of possible combinations n^2: " << number_of_cells << std::endl;

        // iterate over the container, relabel the cells with the bitset, and get their corresponding energy and save it in a vector
        for(int c = 0; c < number_of_possible_combinations; c++){
            for(int v = 0; v < number_of_cells; v++)
                cells_around_nmedge[v]->info().manifold_label = combinations[c].second[v];
            // calculate the energy of the current configuration
            combinations[c].first = nonManifoldCliqueEnergy(Dt, edge, regularization_weight);
        }

        // sort the container while keeping track of its original index, which gives you the corresponding configuration from the bitset index
        std::sort(combinations.begin(), combinations.end(), sortCombinations);

        // take the lowest energy and check if it is manifold
        int solution_found = 0;
        for(int sc = 0; sc < number_of_possible_combinations; sc++){

            // relabel to the correct combination
            for(int v = 0; v < number_of_cells; v++){
                cells_around_nmedge[v]->info().manifold_label = combinations[sc].second[v];
            }
            // NOTE: problem with point-edge-manifold (PEManifold) constraint is that it most of the time does not find valid configurations
            // solution could be to put a new point inside each cell in question, or increase the search space for manifold solutions
            // try first solution first!
//            if(isManifoldClique(Dt, fei)){
            if(isEdgeManifoldClique(Dt, edge, nm_edges)){          // has to be this
//            if(isManifoldEdge(Dt, fei)){              // not this
                // NOTE: in fact, in 3D, flickering one cell in an edge clique CAN change an adjacent
                // edge clique from manifold to non-manifold!!!
                // the furthest outside check is not allowed to change/flicker (more formal: it must keep constant unary weight!)
                // if I want it to change as well, I again need to increase my check for one more level of neighbourhood
                for(int v = 0; v < number_of_cells; v++){
                    cells_around_nmedge[v]->info().gc_label = combinations[sc].second[v];
                }
                solution_found = 1;
                break;
            }
        } // end of trying all possible solutions
        if(!solution_found){
            // std::cout << "no manifold combination found for edge ("
            //           << fei->first->info().global_idx << ", "
            //           << fei->second << ", "
            //           << fei->third << ")"
            //           << std::endl;
            for(int v = 0; v < number_of_cells; v++){
                cells_around_nmedge[v]->info().manifold_label = 2;
            }
        }
        else{
            fixed_count++;
        }

    } // end of iteration over all edges of the Delaunay

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Fixed  " << fixed_count << "/" << nm_edges.size() << " edges in "
                 << full_duration.count() << "s\n" << std::endl;

}


//void getNonManifoldCliques(Delaunay& Dt, std::vector<std::vector<std::vector<Cell_handle>>> nm_cliques){


//    // first collect all non-manifold cells
//    Delaunay::Finite_edges_iterator fei;
//    std::vector<std::vector<int>> v;
//    // TODO: use Delaunay::Simplex to make sets of edges.
//    // https://doc.cgal.org/latest/Triangulation_3/Triangulation_3_2simplex_8cpp-example.html
//    int nm_edge_id = 0;
//    std::vector<Edge> nm_edges;
//    for(fei = Dt.finite_edges_begin(); fei != Dt.finite_edges_end(); fei++){
//        if(isManifoldEdge(Dt, fei))
//            continue;
//        nm_edges.push_back(*fei);
//        std::unordered_set<Cell_handle> cells_around_nmedges;
//        Delaunay::Cell_circulator cc = Dt.incident_cells(*fei);
//        Cell_handle first_cell = cc;
//        do{
//            cc->info.nm_edge_ids.push_back(nm_edge_id);
//            cells_around_nmedges.insert(cc++);
//        }
//        while(first_cell != cc);
//        nm_edge_id++;
//    }

//    for(const auto& cell : cells_around_nmedges){



//    }





//}




