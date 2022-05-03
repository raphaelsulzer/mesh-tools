#include <base/cgal_typedefs.h>

// for GCoptimization
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "gco-v3.0/GCoptimization.h"
#include <processing/meshProcessing.h>
#include <util/geometricOperations.h>

////////////////////////////////////////////////////////////
////////////////////// Feature scaling /////////////////////
////////////////////////////////////////////////////////////


typedef float gtype;

struct ForSmoothFn{
    gtype *data;
};

gtype smoothCost(int s1, int s2, int l1, int l2, void *data){

    ForSmoothFn *myData = (ForSmoothFn *) data;

    // difference of outside observation of s1 and s2
    gtype g = myData->data[s1*2] - myData->data[s2*2];

    gtype smooth_cost;

    if(l1 == 0 && l2 == 1)
        smooth_cost = 0.5*(1+g);
    else if(l1 == 1 && l2 == 0)
        smooth_cost = 0.5*(1-g);
    else
        smooth_cost = 0.5*abs(g);
//        smooth_cost = (exp(-pow(g-1,2) * 2) + exp(-pow(g+1,2) * 2) - 2*exp(-2) ) / 2;

    return smooth_cost;
}


////////////////////////////////////////////////////////////
//////////////////////// Optimization //////////////////////
////////////////////////////////////////////////////////////
//// in this version, set data and smoothness terms using arrays
//// grid neighborhood is set up "manually". Uses spatially varying terms. Namely
//// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with
//// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1
//std::pair<std::map<Cell_handle, int>, std::vector<int>> GeneralGraph_DArraySArraySpatVarying(std::pair<Delaunay&, Cell_map&> dt_cells, std::map<Cell_handle, int>& cell_indexMap, std::vector<int> result, int num_iterations)
void graphCutTet(dataHolder& data, runningOptions options)
{

    gtype div = 100;

    if(options.area_reg_weight > 0.0)
        calculateMeanTriangleArea(data);
    if(options.vol_reg)
        calculateMeanCellVolume(data);

    cout << "\nOptimization..." << endl;
    if(options.area_reg_weight > 0.0){
            options.area_reg_weight/=data.mean_triangle_area;
            cout << "\t-Area regularization weight (/ mean_triangle area):  " << options.area_reg_weight  << endl;
//            cout << "\t-Area regularization weight:  " << options.area_reg_weight  << endl;
    }
    if(options.angle_reg_weight > 0.0)
            cout << "\t-Angle regularization weight:  " << options.angle_reg_weight  << endl;
    if(options.cc_reg_weight > 0.0)
            cout << "\t-Connected components regularization weight: " << options.cc_reg_weight << endl;
    if(options.sv_reg_weight > 0.0)
        cout << "\t-spatially varying weight according to observation with: " << options.sv_reg_weight << endl;
    if(options.ob_reg_weight > 0.0)
        cout << "\t-only binaries with regularization weight: " << options.ob_reg_weight << endl;


    auto start = chrono::high_resolution_clock::now();

    int num_cells = data.Dt.number_of_cells();
    int num_labels = 2;
    GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_cells,num_labels);

    // first set up the array for data costs and set the initial label in the same loop
    gtype *data_term = new gtype[num_cells*num_labels];
    gtype *gt_unaries = new gtype[num_cells*num_labels];
    int idx = 0;
    Delaunay::All_cells_iterator cft;
    // iterate over the all_cells map
    for(cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++)
    {   
        // set an index for each cell
        cft->info().global_idx = idx;

        // set the data term
        if(options.vol_reg){
            if(!data.Dt.is_infinite(cft)){
                double vol = data.Dt.tetrahedron(cft).volume();
                cft->info().outside_score*=(vol/data.mean_cell_volume);
                cft->info().inside_score*=(vol/data.mean_cell_volume);
            }

        }
        if(options.ob_reg_weight > 0.0){
            if(data.Dt.is_infinite(cft)){
                data_term[idx*2+0] = 10;
                data_term[idx*2+1] = 0;
            }
            else{
                data_term[idx*2+0] = 0;
                data_term[idx*2+1] = 0;
            }

        }
        else{
            data_term[idx*2+0] = cft->info().outside_score/div;
            data_term[idx*2+1] = cft->info().inside_score/div;

        }

//        gt_unaries[idx*2+0] = cft->info().gt_outside/div;
//        gt_unaries[idx*2+1] = cft->info().gt_inside/div;


        // set an initial label
        if(cft->info().outside_score > cft->info().inside_score)
            {gc->setLabel(idx, 1);}
        else
            {gc->setLabel(idx, 0);}

        // increase index for next cell
        idx++;
    }

    // next set up the array for smooth costs
    gtype *smooth_term = new gtype[num_labels*num_labels];
    for ( int l1 = 0; l1 < num_labels; l1++ )
        for (int l2 = 0; l2 < num_labels; l2++ )
            if(l1 == l2){smooth_term[l1+l2*num_labels] = 0.0;}
            else{smooth_term[l1+l2*num_labels] = 1.0;}

//    try{

//    ForSmoothFn toFn;
//    toFn.data = gt_unaries;

    gc->setDataCost(data_term);

//    if(options.sv_reg_weight > 0.0)
//        gc->setSmoothCost(&smoothCost,&toFn);
//    else
    gc->setSmoothCost(smooth_term);

    // set neighborhood:
    int current_index;
    int neighbour_index;
    // iterate over all cells
    for(cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++)
    {
        Cell_handle current_cell = cft;
        current_index = current_cell->info().global_idx;
        for(int i = 0; i < 4; i++){

            Cell_handle neighbour_cell = current_cell->neighbor(i);
            neighbour_index = neighbour_cell->info().global_idx;



            // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
            if(neighbour_index < current_index)
                continue;

            // calc smoothness term
            gtype total_local_weight,area_weight,angle_weight,ob_weight;
            Facet current_facet(current_cell, i);
            if(data.Dt.is_infinite(current_facet))
                area_weight = 1;
            else
                area_weight = sqrt(data.Dt.triangle(current_facet).squared_area());
            area_weight*=options.area_reg_weight;

            // this also returns 1 as infinite weight
            angle_weight = 1-std::min(computeCosFacetCellAngle(data.Dt, current_facet),
                                       computeCosFacetCellAngle(data.Dt, data.Dt.mirror_facet(current_facet)));
            angle_weight*=options.angle_reg_weight;


            // unary to binary weight
            ob_weight = 1-abs(current_cell->info().outside_score - neighbour_cell->info().outside_score);
            ob_weight*=options.ob_reg_weight;

            total_local_weight=area_weight+angle_weight+options.cc_reg_weight+options.sv_reg_weight+ob_weight;

            // call the neighbourhood function
            gc->setNeighbors(current_index, neighbour_index, total_local_weight/div);
//            gc->setNeighbors(current_index, neighbour_index);

        }
    }

    // Optimization
    cout << "\t-Before optimization data + smoothness energy: " << endl;
    cout << "\t\t" << gc->giveDataEnergy() << " + " << gc->giveSmoothEnergy() << " = "  << gc->compute_energy() << endl;

    // use swap because it is good when you have two labels
    gc->swap(options.gco_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
//        gc->expansion(num_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);

    // iterate over all cells and get their optimized label
//    for(cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++)
//        cft->info().gc_label = gc->whatLabel(cft->info().global_idx);
//    cout << "\nWARNING: infinite cells are hard coded to be outside\n" << endl;
    for(cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++){
//        if(data.Dt.is_infinite(cft)){
//           cft->info().gc_label = 1;
//        }
//        else{
            cft->info().gc_label = gc->whatLabel(cft->info().global_idx);
//        }
    }

    cout << "\t-After optimization energy is " << gc->compute_energy() << endl;
    delete gc;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;
}

void graphCutFacet(dataHolder& data, runningOptions options)
{
    gtype div = 100;


    // TODO:
    // introduce a mode without node weights (except inf for inf cells maybe) and set binary weights as done in sv


    if(options.area_reg_weight > 0.0)
        calculateMeanTriangleArea(data);
    if(options.vol_reg)
        calculateMeanCellVolume(data);

    cout << "\nOptimization..." << endl;
    if(options.vol_reg){
            cout << "\t-Volume regularization with mean cell volume:  " << data.mean_cell_volume  << endl;
    }
    if(options.area_reg_weight > 0.0){
            options.area_reg_weight/=data.mean_triangle_area;
            cout << "\t-Area regularization weight (/ mean triangle area):  " << options.area_reg_weight  << endl;
    }
    if(options.angle_reg_weight > 0.0)
            cout << "\t-Angle regularization weight:  " << options.angle_reg_weight  << endl;
    if(options.cc_reg_weight > 0.0)
            cout << "\t-Connected components regularization weight: " << options.cc_reg_weight << endl;

    auto start = chrono::high_resolution_clock::now();

    int num_cells = data.Dt.number_of_cells();
    int num_labels = 2;
    GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_cells,num_labels);

    // first set up the array for data costs and set the initial label in the same loop
    gtype *data_term = new gtype[num_cells*num_labels];
    int idx = 0;
    Delaunay::All_cells_iterator cft;
    // iterate over the all_cells map
    for(cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++)
    {

        if(options.vol_reg){
            if(!data.Dt.is_infinite(cft)){
                double vol = data.Dt.tetrahedron(cft).volume();
                cft->info().outside_score*=(vol/data.mean_cell_volume);
                cft->info().inside_score*=(vol/data.mean_cell_volume);
            }

        }

        // set an index for each cell
//        cft->info().global_idx = idx;
        idx = cft->info().global_idx;
        // I am initializing my label s.t. if the outside vote is bigger than the inside vote, than it should be labelled 1 (outside) - and vice versa (0 = inside)
        // that means the cost/penalty for labelling a cell that initially has label 1 with the opposite label, is the opposite vote (so the inside vote)
        // so labelling 0 costs outside votes
        data_term[idx*2+0] = (gtype)cft->info().outside_score/div;
        // and labelling 1 costs inside votes
        data_term[idx*2+1] = (gtype)cft->info().inside_score/div;

//        cout << " inside: " << cft->info().inside_score << " outside: " << cft->info().outside_score << endl;

        // set an initial label
        // this makes sense and should not be >=, because cells that are not crossed at all
        // are much more likely to be inside, and should thus be initialized as inside.
        if(cft->info().outside_score >= cft->info().inside_score)
            {gc->setLabel(idx, 1);}
        else
            {gc->setLabel(idx, 0);}
    }

    // next set up the array for smooth costs
    gtype *smooth_term = new gtype[num_labels*num_labels];
    for(int l1 = 0; l1 < num_labels; l1++)
        for(int l2 = 0; l2 < num_labels; l2++)
            if(l1 == l2){smooth_term[l1+l2*num_labels] = 0;}
            else{smooth_term[l1+l2*num_labels] = 1;}

    gc->setDataCost(data_term);
    gc->setSmoothCost(smooth_term);

    // set neighborhood:
    int current_index;
    int neighbour_index;
    // iterate over all cells
    for(cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++)
    {
        Cell_handle current_cell = cft;
        current_index = current_cell->info().global_idx;
        // iterate over all facets of the cellF
        for(int i = 0; i < 4; i++){

            Cell_handle neighbour_cell = current_cell->neighbor(i);
            neighbour_index = neighbour_cell->info().global_idx;

            // labatut version is not working properly if this is activated
            // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
//            if(neighbour_index < current_index)
//                continue;

            // calc smoothness term
            gtype area_weight,angle_weight;
            Facet current_facet(current_cell, i);
            if(data.Dt.is_infinite(current_facet))
                area_weight = 1;
            else
                area_weight = sqrt(data.Dt.triangle(current_facet).squared_area());
            area_weight*=options.area_reg_weight;

            // this also returns 1 if the first cell is infinite
            angle_weight = 1-std::min(computeCosFacetCellAngle(data.Dt, current_facet),
                                       computeCosFacetCellAngle(data.Dt, data.Dt.mirror_facet(current_facet)));
            // if both cells are infinite, still put an edge weight of 1 on there, to prevent making an interface between infinite cells
            if(data.Dt.is_infinite(current_cell) && data.Dt.is_infinite(neighbour_cell))
                angle_weight = 1;
            angle_weight*=options.angle_reg_weight;

            // call the neighbourhood function
            gc->setNeighbors(current_index, neighbour_index, (area_weight+angle_weight+current_cell->info().facet_weights[i])/div);


//            cout << "binary: " << shape_weight+current_cell->info().facet_weights[i] << endl;

        }
    }

    // Optimization
    cout << "\t-Before optimization data + smoothness energy: " << endl;
    cout << "\t\t" << gc->giveDataEnergy() << " + " << gc->giveSmoothEnergy() << " = "  << gc->compute_energy() << endl;

    // use swap because it is good when you have two labels
    gc->swap(options.gco_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
//        gc->expansion(num_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);

    // iterate over all cells and get their optimized label
//    cout << "\nWARNING: infinite cells are hard coded to be outside" << endl;
    for(cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++){
//        if(data.Dt.is_infinite(cft)){
//           cft->info().gc_label = 1;
//        }
//        else{
            cft->info().gc_label = gc->whatLabel(cft->info().global_idx);
//        }
    }

    cout << "\t-After optimization energy is " << gc->compute_energy() << endl;
    delete gc;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-done in " << duration.count() << "s" << endl;
}
