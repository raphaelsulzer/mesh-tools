#include "learning/learningIO.h"
#include <learning/learningMath.h>
#include <processing/graphCut.h>
#include <util/helper.h>
#include <cstdlib>
#include <ctime>
#include <CGAL/Random.h>
#include <boost/filesystem.hpp>

#include <xtensor/xbuilder.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor-io/xnpz.hpp>

namespace fs = boost::filesystem;

///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// EXPORT FROM POISSON INPUT ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////



void cellIndexAndLabel(ofstream& out, Delaunay& Dt, Delaunay::All_cells_iterator& fci, int idx, runningOptions options){

    if(options.ground_truth){
        //    il << "global_index inside_perc outside_perc inside_outside_max gc_label" << endl;
        out << idx << " ";
        out << fci->info().inside_score << " ";
        out << fci->info().outside_score << " ";
        // this means 50/50 cells will be labelled as inside
        int l = fci->info().outside_score > fci->info().inside_score ? 1 : 0;
        out << l << " ";
        out << fci->info().gc_label << " ";
        out << Dt.is_infinite(fci) << endl;
    }
    else
        out << "999 999 999 999 999 " << Dt.is_infinite(fci) << endl;



}
void cellBasedGeometricFeatures(ofstream& out, Delaunay& Dt, Delaunay::All_cells_iterator& fci){

//        cgeom << "radius vol longest_edge shortest_edge" << endl;
    if(Dt.is_infinite(fci))
        out << "0 0 0 0" << endl;
    else{
        auto tet = Dt.tetrahedron(fci);
        tetFeatures tf = calcTetFeatures(tet);
        out << tf.radius << " ";
        out << tf.vol << " ";
        out << tf.longest_edge << " ";
        out << tf.shortest_edge << " ";
        out << endl;
    }
}

void cellBasedVertexFeatures(ofstream& out, Delaunay& Dt, Delaunay::All_cells_iterator& fci){

//    cbvf << "cb_vertex_inside_count cb_vertex_inside_dist_min cb_vertex_inside_dist_max cb_vertex_inside_dist_sum ";
//    cbvf << "cb_vertex_outside_count cb_vertex_outside_dist_min cb_vertex_outside_dist_max cb_vertex_outside_dist_sum ";
//    cbvf << "cb_vertex_last_count cb_vertex_last_dist_min cb_vertex_last_dist_max cb_vertex_last_dist_sum" << endl;

    // inside
    int vdi = fci->info().cb_vertex_inside.size();
    out << vdi << " ";
    if(vdi > 0){
        out << *min_element(fci->info().cb_vertex_inside.begin(), fci->info().cb_vertex_inside.end()) << " ";
        out << *max_element(fci->info().cb_vertex_inside.begin(), fci->info().cb_vertex_inside.end()) << " ";
        out << accumulate(fci->info().cb_vertex_inside.begin(), fci->info().cb_vertex_inside.end(),0.0) << " ";

    }
    else
        out << "0 0 0 ";
    // outside
    int vdo = fci->info().cb_vertex_outside.size();
    out << vdo << " ";
    if(vdo > 0){
        out << *min_element(fci->info().cb_vertex_outside.begin(), fci->info().cb_vertex_outside.end()) << " ";
        out << *max_element(fci->info().cb_vertex_outside.begin(), fci->info().cb_vertex_outside.end()) << " ";
        out << accumulate(fci->info().cb_vertex_outside.begin(), fci->info().cb_vertex_outside.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // last
    int vdl = fci->info().cb_vertex_last.size();
    out << vdl << " ";
    if(vdl > 0){
        out << *min_element(fci->info().cb_vertex_last.begin(), fci->info().cb_vertex_last.end()) << " ";
        out << *max_element(fci->info().cb_vertex_last.begin(), fci->info().cb_vertex_last.end()) << " ";
        out << accumulate(fci->info().cb_vertex_last.begin(), fci->info().cb_vertex_last.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";

    out << endl;
}

void cellBasedFacetFeatures(ofstream& out, Delaunay& Dt, Delaunay::All_cells_iterator& fci){

//    cbff << "cb_facet_inside_first_count cb_facet_inside_first_dist_min cb_facet_inside_first_dist_max cb_facet_inside_first_dist_sum ";
//    cbff << "cb_facet_inside_second_count cb_facet_inside_second_dist_min cb_facet_inside_second_dist_max cb_facet_inside_second_dist_sum ";
//    cbff << "cb_facet_outside_first_count cb_facet_outside_first_dist_min cb_facet_outside_first_dist_max cb_facet_outside_first_dist_sum ";
//    cbff << "cb_facet_outside_second_count cb_facet_outside_second_dist_min cb_facet_outside_second_dist_max cb_facet_outside_second_dist_sum ";
//    cbff << "cb_facet_last_first_count cb_facet_last_first_dist_min cb_facet_last_first_dist_max cb_facet_last_first_dist_sum" << endl;
//    cbff << "cb_facet_last_second_count cb_facet_last_second_dist_min cb_facet_last_second_dist_max cb_facet_last_second_dist_sum" << endl;

    // inside first
    int fdif = fci->info().cb_facet_inside_first.size();
    out << fdif << " ";
    if(fdif > 0){
        out << *min_element(fci->info().cb_facet_inside_first.begin(), fci->info().cb_facet_inside_first.end()) << " ";
        out << *max_element(fci->info().cb_facet_inside_first.begin(), fci->info().cb_facet_inside_first.end()) << " ";
        out << accumulate(fci->info().cb_facet_inside_first.begin(), fci->info().cb_facet_inside_first.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // inside second
    int fdis = fci->info().cb_facet_inside_second.size();
    out << fdis << " ";
    if(fdis > 0){
        out << *min_element(fci->info().cb_facet_inside_second.begin(), fci->info().cb_facet_inside_second.end()) << " ";
        out << *max_element(fci->info().cb_facet_inside_second.begin(), fci->info().cb_facet_inside_second.end()) << " ";
        out << accumulate(fci->info().cb_facet_inside_second.begin(), fci->info().cb_facet_inside_second.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // outside first
    int fdof = fci->info().cb_facet_outside_first.size();
    out << fdof << " ";
    if(fdof > 0){
        out << *min_element(fci->info().cb_facet_outside_first.begin(), fci->info().cb_facet_outside_first.end()) << " ";
        out << *max_element(fci->info().cb_facet_outside_first.begin(), fci->info().cb_facet_outside_first.end()) << " ";
        out << accumulate(fci->info().cb_facet_outside_first.begin(), fci->info().cb_facet_outside_first.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // outside second
    int fdos = fci->info().cb_facet_outside_second.size();
    out << fdos << " ";
    if(fdos > 0){
        out << *min_element(fci->info().cb_facet_outside_second.begin(), fci->info().cb_facet_outside_second.end()) << " ";
        out << *max_element(fci->info().cb_facet_outside_second.begin(), fci->info().cb_facet_outside_second.end()) << " ";
        out << accumulate(fci->info().cb_facet_outside_second.begin(), fci->info().cb_facet_outside_second.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // last first
    int fdlf = fci->info().cb_facet_last_first.size();
    out << fdlf << " ";
    if(fdlf > 0){
        out << *min_element(fci->info().cb_facet_last_first.begin(), fci->info().cb_facet_last_first.end()) << " ";
        out << *max_element(fci->info().cb_facet_last_first.begin(), fci->info().cb_facet_last_first.end()) << " ";
        out << accumulate(fci->info().cb_facet_last_first.begin(), fci->info().cb_facet_last_first.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // last second
    int fdls = fci->info().cb_facet_last_second.size();
    out << fdls << " ";
    if(fdls > 0){
        out << *min_element(fci->info().cb_facet_last_second.begin(), fci->info().cb_facet_last_second.end()) << " ";
        out << *max_element(fci->info().cb_facet_last_second.begin(), fci->info().cb_facet_last_second.end()) << " ";
        out << accumulate(fci->info().cb_facet_last_second.begin(), fci->info().cb_facet_last_second.end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";

    out << endl;
}

void facetBasedGeometricFeatures(ofstream& out, Delaunay& Dt, Delaunay::All_cells_iterator& fci, int j){

//    fgeom << "area angle beta dist" << endl;
    // export the area and angle off the current facet as edge feature
    auto fac = make_pair(fci,j);
    out << computeFacetArea(Dt,fac) << " ";
    out << computeCosFacetCellAngle(Dt,fac) << " ";
    out << 1-min(computeCosFacetCellAngle(Dt, fac),
                             computeCosFacetCellAngle(Dt, Dt.mirror_facet(fac))) << " ";

    auto neighbor = fci->neighbor(j);
    if(!Dt.is_infinite(neighbor) && !Dt.is_infinite(fci)){
        auto cen1 = CGAL::centroid(fci->vertex(0)->point(),
                                   fci->vertex(1)->point(),
                                   fci->vertex(2)->point(),
                                   fci->vertex(3)->point());
        auto cen2 = CGAL::centroid(neighbor->vertex(0)->point(),
                       neighbor->vertex(1)->point(),
                       neighbor->vertex(2)->point(),
                       neighbor->vertex(3)->point());
        out << CGAL::squared_distance(cen1,cen2);
    }
    else
        out << "-1";
    out << endl;
}
void facetBasedVertexFeatures(ofstream& out, Delaunay& Dt, Delaunay::All_cells_iterator& fci, int j){

//    fbvf << "fb_vertex_inside_count fb_vertex_inside_dist_min fb_vertex_inside_dist_max fb_vertex_inside_dist_sum ";
//    fbvf << "fb_vertex_outside_count fb_vertex_outside_dist_min fb_vertex_outside_dist_max fb_vertex_outside_dist_sum" << endl;
//    fbvf << "fb_vertex_last_count fb_vertex_last_dist_min fb_vertex_last_dist_max fb_vertex_last_dist_sum" << endl;

    // inside
    int vdi = fci->info().fb_vertex_inside[j].size();
    out << vdi << " ";
    if(vdi > 0){
        out << *min_element(fci->info().fb_vertex_inside[j].begin(), fci->info().fb_vertex_inside[j].end()) << " ";
        out << *max_element(fci->info().fb_vertex_inside[j].begin(), fci->info().fb_vertex_inside[j].end()) << " ";
        out << accumulate(fci->info().fb_vertex_inside[j].begin(), fci->info().fb_vertex_inside[j].end(),0.0) << " ";

    }
    else
        out << "0 0 0 ";
    // outside
    int vdo = fci->info().fb_vertex_outside[j].size();
    out << vdo << " ";
    if(vdo > 0){
        out << *min_element(fci->info().fb_vertex_outside[j].begin(), fci->info().fb_vertex_outside[j].end()) << " ";
        out << *max_element(fci->info().fb_vertex_outside[j].begin(), fci->info().fb_vertex_outside[j].end()) << " ";
        out << accumulate(fci->info().fb_vertex_outside[j].begin(), fci->info().fb_vertex_outside[j].end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // last
    int vdl = fci->info().fb_vertex_last[j].size();
    out << vdl << " ";
    if(vdl > 0){
        out << *min_element(fci->info().fb_vertex_last[j].begin(), fci->info().fb_vertex_last[j].end()) << " ";
        out << *max_element(fci->info().fb_vertex_last[j].begin(), fci->info().fb_vertex_last[j].end()) << " ";
        out << accumulate(fci->info().fb_vertex_last[j].begin(), fci->info().fb_vertex_last[j].end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";

    out << endl;

}

void facetBasedFacetFeatures(ofstream& out, Delaunay& Dt, Delaunay::All_cells_iterator& fci, int j){

//    fbff << "fb_facet_inside_count fb_facet_inside_dist_min fb_facet_inside_dist_max fb_facet_inside_dist_sum ";
//    fbff << "fb_facet_outside_count fb_facet_outside_dist_min fb_facet_outside_dist_max fb_facet_outside_dist_sum" << endl;
//    fbff << "fb_facet_last_count fb_facet_last_dist_min fb_facet_last_dist_max fb_facet_last_dist_sum" << endl;

    // inside
    int vdi = fci->info().fb_facet_inside[j].size();
    out << vdi << " ";
    if(vdi > 0){
        out << *min_element(fci->info().fb_facet_inside[j].begin(), fci->info().fb_facet_inside[j].end()) << " ";
        out << *max_element(fci->info().fb_facet_inside[j].begin(), fci->info().fb_facet_inside[j].end()) << " ";
        out << accumulate(fci->info().fb_facet_inside[j].begin(), fci->info().fb_facet_inside[j].end(),0.0) << " ";

    }
    else
        out << "0 0 0 ";
    // outside
    int vdo = fci->info().fb_facet_outside[j].size();
    out << vdo << " ";
    if(vdo > 0){
        out << *min_element(fci->info().fb_facet_outside[j].begin(), fci->info().fb_facet_outside[j].end()) << " ";
        out << *max_element(fci->info().fb_facet_outside[j].begin(), fci->info().fb_facet_outside[j].end()) << " ";
        out << accumulate(fci->info().fb_facet_outside[j].begin(), fci->info().fb_facet_outside[j].end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";
    // last
    int vdl = fci->info().fb_facet_last[j].size();
    out << vdl << " ";
    if(vdl > 0){
        out << *min_element(fci->info().fb_facet_last[j].begin(), fci->info().fb_facet_last[j].end()) << " ";
        out << *max_element(fci->info().fb_facet_last[j].begin(), fci->info().fb_facet_last[j].end()) << " ";
        out << accumulate(fci->info().fb_facet_last[j].begin(), fci->info().fb_facet_last[j].end(),0.0) << " ";
    }
    else
        out << "0 0 0 ";

    out << endl;

}





void exportGraph(dirHolder& dir, runningOptions& options, Delaunay& Dt){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nExport graph..." << endl;
    boost::filesystem::path p(dir.write_file);
    string outfile = p.stem().string();
    cout << "\t-to " << "gt/"+outfile+"_X.txt" << endl;
    string rays = (options.export_rays) ? "\t-with rays" : "\t-without rays";
    cout << rays << endl;
    int print_precision = 8;


    ////////////////////////////////////////////
    //////////////    HEADERS    ///////////////
    ////////////////////////////////////////////

    /////////////////////////////////////////////
    //////////////////// CELLS //////////////////
    /////////////////////////////////////////////
    ////// info file
    ofstream info;
    info.open(dir.path+"gt/"+outfile+"_info.txt");
    info << "vertices " << Dt.number_of_vertices() << endl;
    info << "facets " << Dt.number_of_facets() << endl;
    info << "cells " << Dt.number_of_cells() << endl;
    info.close();

    ////// label and index
    // cell index and labels
    ofstream il;
    il.open(dir.path+"gt/"+outfile+"_labels.txt");
    il << "global_index inside_perc outside_perc inside_outside_max gc_label infinite" << endl;

    ////// geometric features
    ofstream cgeom;
    cgeom.open(dir.path+"gt/"+outfile+"_cgeom.txt");
    cgeom << "radius vol longest_edge shortest_edge" << endl;
    cgeom << setprecision(print_precision);

    ////// cell based vertex features
    ofstream cbvf;
    cbvf.open(dir.path+"gt/"+outfile+"_cbvf.txt");
    cbvf << "cb_vertex_inside_count cb_vertex_inside_dist_min cb_vertex_inside_dist_max cb_vertex_inside_dist_sum ";
    cbvf << "cb_vertex_outside_count cb_vertex_outside_dist_min cb_vertex_outside_dist_max cb_vertex_outside_dist_sum ";
    cbvf << "cb_vertex_last_count cb_vertex_last_dist_min cb_vertex_last_dist_max cb_vertex_last_dist_sum" << endl;
    cbvf << setprecision(print_precision);

    ////// cell based facet features
    ofstream cbff;
    cbff.open(dir.path+"gt/"+outfile+"_cbff.txt");
    cbff << "cb_facet_inside_first_count cb_facet_inside_first_dist_min cb_facet_inside_first_dist_max cb_facet_inside_first_dist_sum ";
    cbff << "cb_facet_inside_second_count cb_facet_inside_second_dist_min cb_facet_inside_second_dist_max cb_facet_inside_second_dist_sum ";
    cbff << "cb_facet_outside_first_count cb_facet_outside_first_dist_min cb_facet_outside_first_dist_max cb_facet_outside_first_dist_sum ";
    cbff << "cb_facet_outside_second_count cb_facet_outside_second_dist_min cb_facet_outside_second_dist_max cb_facet_outside_second_dist_sum ";
    cbff << "cb_facet_last_first_count cb_facet_last_first_dist_min cb_facet_last_first_dist_max cb_facet_last_first_dist_sum ";
    cbff << "cb_facet_last_second_count cb_facet_last_second_dist_min cb_facet_last_second_dist_max cb_facet_last_second_dist_sum" << endl;
    cbff << setprecision(print_precision);


    /////////////////////////////////////////////
    //////////////////// EDGES //////////////////
    /////////////////////////////////////////////
    // adjacency files
    ofstream adjacency_ij, adjacency_ji;
    adjacency_ij.open(dir.path+"gt/"+outfile+"_adjacency_ij.txt");
    adjacency_ij << "# list of edges" << endl;
    adjacency_ji.open(dir.path+"gt/"+outfile+"_adjacency_ji.txt");
    adjacency_ji << "# list of edges" << endl;

    // geometric features
    ofstream fgeom;
    fgeom.open(dir.path+"gt/"+outfile+"_fgeom.txt");
    fgeom << "area angle beta dist" << endl;
    // facet based vertex features
    ofstream fbvf;
    fbvf.open(dir.path+"gt/"+outfile+"_fbvf.txt");
    fbvf << "fb_vertex_inside_count fb_vertex_inside_dist_min fb_vertex_inside_dist_max fb_vertex_inside_dist_sum ";
    fbvf << "fb_vertex_outside_count fb_vertex_outside_dist_min fb_vertex_outside_dist_max fb_vertex_outside_dist_sum ";
    fbvf << "fb_vertex_last_count fb_vertex_last_dist_min fb_vertex_last_dist_max fb_vertex_last_dist_sum" << endl;

    // facet based facet features
    ofstream fbff;
    fbff.open(dir.path+"gt/"+outfile+"_fbff.txt");
    fbff << "fb_facet_inside_count fb_facet_inside_dist_min fb_facet_inside_dist_max fb_facet_inside_dist_sum ";
    fbff << "fb_facet_outside_count fb_facet_outside_dist_min fb_facet_outside_dist_max fb_facet_outside_dist_sum ";
    fbff << "fb_facet_last_count fb_facet_last_dist_min fb_facet_last_dist_max fb_facet_last_dist_sum" << endl;

    int edge_count = 0;
    int i, j;

    string npz_out = dir.path+"gt/"+outfile+"_npz.npz";
    Delaunay::All_cells_iterator fci;
    for(fci = Dt.all_cells_begin(); fci != Dt.all_cells_end(); fci++){

        // TODO: there is a problem here, not all cells are exported, running index i < nc at the end;
        i = fci->info().global_idx;

        /////////////////////////////////////////////
        //////////////     LABELS    ////////////////
        /////////////////////////////////////////////
        cellIndexAndLabel(il, Dt, fci, i, options);

        /////////////////////////////////////////////
        ////////////// CELL FEATURES ////////////////
        /////////////////////////////////////////////
        cellBasedGeometricFeatures(cgeom, Dt, fci);
        cellBasedVertexFeatures(cbvf,Dt,fci);
        cellBasedFacetFeatures(cbff,Dt,fci);

//        xt::xarray<double> radius;
//        xt::xarray<double> vol;
//        if(Dt.is_infinite(fci)){
//            radius = {0.0};
//            vol = {0.0};
//        }
//        else{
//            auto tet = Dt.tetrahedron(fci);
//            tetFeatures tf = calcTetFeatures(tet);
//            radius = {tf.radius};
//            vol = {tf.vol};
//        }

//        xt::dump_npz(npz_out,"radius",radius,true,true);
//        xt::dump_npz(npz_out,"volume",vol,true,true);


        /////////////////////////////////////////////
        ////////////// EDGE FEATURES ////////////////
        /////////////////////////////////////////////
        for(int c = 0; c < 4; c++){
            // set the edge to the mirror cell
            j = fci->neighbor(c)->info().global_idx;

            if(j == i) // skip self loops
                continue;
            adjacency_ij << i << " ";
            adjacency_ji << j << " ";

            //// geometric features
            facetBasedGeometricFeatures(fgeom,Dt,fci,c);
            facetBasedVertexFeatures(fbvf,Dt,fci,c);
            facetBasedFacetFeatures(fbff,Dt,fci,c);

            edge_count++;
        }
    }

    il.close();
    cgeom.close();
    cbvf.close();
    cbff.close();

    adjacency_ij.close();
    adjacency_ji.close();
    fgeom.close();
    fbvf.close();
    fbff.close();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Exported " << i << " nodes and " << edge_count << " edges" << endl;
    cout << "\t-in " << duration.count() << "s" << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// IMPORT ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
int loadPrediction(dirHolder dir, dataHolder& data, runningOptions options){


    fs::path path = fs::path(dir.path) / fs::path(dir.prediction_file);
    ifstream file(path.string());

    cout << "\nLoad prediction score..." << endl;
    cout << "\t-from " << path.string() << endl;

    if(!file){
        cout << "\nFILE DOES NOT EXIST OR IS EMPTY!" << endl;
        return 1;
    }
    auto npz_map = xt::load_npz(path.string());

    Delaunay::All_cells_iterator fci;
    int b = data.Dt.number_of_cells();
    auto nc = npz_map["number_of_cells"].cast<long>();
    const int a = nc(0,0);
    // check for same sampling of input and scan
    if(a != b){
        cout << "ERROR: SAMPLING IS NOT THE SAME!" << endl;
        cout << "\t-" << b << " Delaunay cells" << endl;
        cout << "\t-" << a << " imported cells" << endl;
        return 1;
    }

    int i = 0;
    array<size_t, 2> shape = { a, 2 };
    xt::xtensor<float, 2> pred(shape);
    if(options.prediction_type == "lo"){
        pred = npz_map["logits"].cast<float>();
    }
    else if(options.prediction_type == "sm"){
        pred = npz_map["softmax"].cast<float>();
    }
    else if(options.prediction_type == "si"){
        pred = npz_map["sigmoid"].cast<float>();
    }
    else{
        cout << options.prediction_type << " is not a valid prediction type, choose either lo or sm" << endl;
        return 1;
    }
    for(fci = data.Dt.all_cells_begin(); fci != data.Dt.all_cells_end(); fci++){

        if(data.Dt.is_infinite(fci)){
            fci->info().outside_score = 1.0;
            fci->info().inside_score = 0.0;
        }
        else{

            // check if there is a camera inside the cell
            for(auto s = data.sensor_map.begin(); s != data.sensor_map.end(); s++){
                if(data.Dt.tetrahedron(fci).has_on_positive_side(s->second)){
                    fci->info().outside_score = 1000;
                    fci->info().inside_score = -1000;
                    break;
                }
            }


//            xt::xtensor_fixed<double, xt::xshape<a, 2>> pred;
            if(pred.shape(1)==1){
//                cout << pred(i,0) << endl;
                if(pred(i,0)<=0.0){
                    fci->info().outside_score = 0;
                    fci->info().inside_score = 1;
                }
                else{
                    fci->info().outside_score = 1;
                    fci->info().inside_score = 0;
                }
            }
            else{
                fci->info().inside_score = pred(i,0);
                fci->info().outside_score = pred(i,1);
                // TODO: maybe it should be
                // fci->info().inside_score = pred(i,0) - pred(i,1);
                // fci->info().outside_score = pred(i,1) - pred(i,0);
                // because in the case of logits pred(i,0) = 1 - pred(i,1)
                // UPDATE: in fact it doesn't matter! what matters in  the graph cut is the difference between the two labelling options.
            }
        }
        i++;

    }
    cout<< "\t-read " << a << " cell scores. " << endl;

    return 0;
}

