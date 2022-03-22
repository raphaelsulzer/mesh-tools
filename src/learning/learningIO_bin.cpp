#include "learning/learningIO_bin.h"
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
void graphExporter::cellLabel(Delaunay::All_cells_iterator& fci){

    int i = fci->info().global_idx;
    infinite_[i] = Dt_.is_infinite(fci);

    if(options_.ground_truth){
        inside_perc_[i] = fci->info().inside_score;
        outside_perc_[i] = fci->info().outside_score;
    }

}
void graphExporter::cellGeometricFeatures(Delaunay::All_cells_iterator& fci){

    int i = fci->info().global_idx;
    if(!Dt_.is_infinite(fci)){
        auto tet = Dt_.tetrahedron(fci);
        tetFeatures tf = calcTetFeatures(tet);
        radius_[i] = tf.radius;
        vol_[i] = tf.vol;
        longest_edge_[i] = tf.longest_edge;
        shortest_edge_[i] = tf.shortest_edge;
    }
}

void graphExporter::cellVertexFeatures(Delaunay::All_cells_iterator& fci){

//    cbvf << "cb_vertex_inside_count cb_vertex_inside_dist_min cb_vertex_inside_dist_max cb_vertex_inside_dist_sum ";
//    cbvf << "cb_vertex_outside_count cb_vertex_outside_dist_min cb_vertex_outside_dist_max cb_vertex_outside_dist_sum ";
//    cbvf << "cb_vertex_last_count cb_vertex_last_dist_min cb_vertex_last_dist_max cb_vertex_last_dist_sum" << endl;

    // inside
    int i = fci->info().global_idx;
    int vdi = fci->info().cb_vertex_inside.size();
    cb_vertex_inside_count_[i] = vdi;
    if(vdi > 0){
        cb_vertex_inside_dist_min_[i] = *min_element(fci->info().cb_vertex_inside.begin(), fci->info().cb_vertex_inside.end());
        cb_vertex_inside_dist_max_[i] = *max_element(fci->info().cb_vertex_inside.begin(), fci->info().cb_vertex_inside.end());
        cb_vertex_inside_dist_sum_[i] = accumulate(fci->info().cb_vertex_inside.begin(), fci->info().cb_vertex_inside.end(),0.0);
    }

    // outside
    int vdo = fci->info().cb_vertex_outside.size();
    cb_vertex_outside_count_[i] = vdo;
    if(vdo > 0){
        cb_vertex_outside_dist_min_[i] = *min_element(fci->info().cb_vertex_outside.begin(), fci->info().cb_vertex_outside.end());
        cb_vertex_outside_dist_max_[i] = *max_element(fci->info().cb_vertex_outside.begin(), fci->info().cb_vertex_outside.end());
        cb_vertex_outside_dist_sum_[i] = accumulate(fci->info().cb_vertex_outside.begin(), fci->info().cb_vertex_outside.end(),0.0);
    }

    // last
    int vdl = fci->info().cb_vertex_last.size();
    cb_vertex_last_count_[i] = vdl;
    if(vdl > 0){
        cb_vertex_last_dist_min_[i] = *min_element(fci->info().cb_vertex_last.begin(), fci->info().cb_vertex_last.end());
        cb_vertex_last_dist_max_[i] = *max_element(fci->info().cb_vertex_last.begin(), fci->info().cb_vertex_last.end());
        cb_vertex_last_dist_sum_[i] = accumulate(fci->info().cb_vertex_last.begin(), fci->info().cb_vertex_last.end(),0.0);
    }
}

void graphExporter::cellFacetFeatures(Delaunay::All_cells_iterator& fci){

//    cbff << "cb_facet_inside_first_count cb_facet_inside_first_dist_min cb_facet_inside_first_dist_max cb_facet_inside_first_dist_sum ";
//    cbff << "cb_facet_inside_second_count cb_facet_inside_second_dist_min cb_facet_inside_second_dist_max cb_facet_inside_second_dist_sum ";
//    cbff << "cb_facet_outside_first_count cb_facet_outside_first_dist_min cb_facet_outside_first_dist_max cb_facet_outside_first_dist_sum ";
//    cbff << "cb_facet_outside_second_count cb_facet_outside_second_dist_min cb_facet_outside_second_dist_max cb_facet_outside_second_dist_sum ";
//    cbff << "cb_facet_last_first_count cb_facet_last_first_dist_min cb_facet_last_first_dist_max cb_facet_last_first_dist_sum" << endl;
//    cbff << "cb_facet_last_second_count cb_facet_last_second_dist_min cb_facet_last_second_dist_max cb_facet_last_second_dist_sum" << endl;

    // inside first
    int i = fci->info().global_idx;
    int fdif = fci->info().cb_facet_inside_first.size();
    cb_facet_inside_first_count_[i] = fdif;
    if(fdif > 0){
        cb_facet_inside_first_dist_min_[i] = *min_element(fci->info().cb_facet_inside_first.begin(), fci->info().cb_facet_inside_first.end());
        cb_facet_inside_first_dist_max_[i] = *max_element(fci->info().cb_facet_inside_first.begin(), fci->info().cb_facet_inside_first.end());
        cb_facet_inside_first_dist_sum_[i] = accumulate(fci->info().cb_facet_inside_first.begin(), fci->info().cb_facet_inside_first.end(),0.0);
    }
    // inside second
    int fdis = fci->info().cb_facet_inside_second.size();
    cb_facet_inside_second_count_[i] = fdis;
    if(fdis > 0){
        cb_facet_inside_second_dist_min_[i] = *min_element(fci->info().cb_facet_inside_second.begin(), fci->info().cb_facet_inside_second.end());
        cb_facet_inside_second_dist_max_[i] = *max_element(fci->info().cb_facet_inside_second.begin(), fci->info().cb_facet_inside_second.end());
        cb_facet_inside_second_dist_sum_[i] = accumulate(fci->info().cb_facet_inside_second.begin(), fci->info().cb_facet_inside_second.end(),0.0);
    }
    // outside first
    int fdof = fci->info().cb_facet_outside_first.size();
    cb_facet_outside_first_count_[i] = fdof;
    if(fdof > 0){
        cb_facet_outside_first_dist_min_[i] = *min_element(fci->info().cb_facet_outside_first.begin(), fci->info().cb_facet_outside_first.end());
        cb_facet_outside_first_dist_max_[i] = *max_element(fci->info().cb_facet_outside_first.begin(), fci->info().cb_facet_outside_first.end());
        cb_facet_outside_first_dist_sum_[i] = accumulate(fci->info().cb_facet_outside_first.begin(), fci->info().cb_facet_outside_first.end(),0.0);
    }
    // outside second
    int fdos = fci->info().cb_facet_outside_second.size();
    cb_facet_outside_second_count_[i] = fdos;
    if(fdos > 0){
        cb_facet_outside_second_dist_min_[i] = *min_element(fci->info().cb_facet_outside_second.begin(), fci->info().cb_facet_outside_second.end());
        cb_facet_outside_second_dist_max_[i] = *max_element(fci->info().cb_facet_outside_second.begin(), fci->info().cb_facet_outside_second.end());
        cb_facet_outside_second_dist_sum_[i] = accumulate(fci->info().cb_facet_outside_second.begin(), fci->info().cb_facet_outside_second.end(),0.0);
    }
    // last first
    int fdlf = fci->info().cb_facet_last_first.size();
    cb_facet_last_first_count_[i] = fdlf;
    if(fdlf > 0){
        cb_facet_last_first_dist_min_[i] = *min_element(fci->info().cb_facet_last_first.begin(), fci->info().cb_facet_last_first.end());
        cb_facet_last_first_dist_max_[i] = *max_element(fci->info().cb_facet_last_first.begin(), fci->info().cb_facet_last_first.end());
        cb_facet_last_first_dist_sum_[i] = accumulate(fci->info().cb_facet_last_first.begin(), fci->info().cb_facet_last_first.end(),0.0);
    }
    // last second
    int fdls = fci->info().cb_facet_last_second.size();
    cb_facet_last_second_count_[i] = fdls;
    if(fdls > 0){
        cb_facet_last_second_dist_min_[i] = *min_element(fci->info().cb_facet_last_second.begin(), fci->info().cb_facet_last_second.end());
        cb_facet_last_second_dist_max_[i] = *max_element(fci->info().cb_facet_last_second.begin(), fci->info().cb_facet_last_second.end());
        cb_facet_last_second_dist_sum_[i] = accumulate(fci->info().cb_facet_last_second.begin(), fci->info().cb_facet_last_second.end(),0.0);
    }
}

void graphExporter::facetGeometricFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count){

//    fgeom << "area angle beta dist" << endl;
    // export the area and angle off the current facet as edge feature
    auto fac = make_pair(fci,j);
    area_[edge_count] = computeFacetArea(Dt_,fac);
    angle_[edge_count] = computeCosFacetCellAngle(Dt_,fac);
    beta_[edge_count] = 1-min(computeCosFacetCellAngle(Dt_, fac),
                             computeCosFacetCellAngle(Dt_, Dt_.mirror_facet(fac)));

    auto neighbor = fci->neighbor(j);
    if(!Dt_.is_infinite(neighbor) && !Dt_.is_infinite(fci)){
        auto cen1 = CGAL::centroid(fci->vertex(0)->point(),
                                   fci->vertex(1)->point(),
                                   fci->vertex(2)->point(),
                                   fci->vertex(3)->point());
        auto cen2 = CGAL::centroid(neighbor->vertex(0)->point(),
                       neighbor->vertex(1)->point(),
                       neighbor->vertex(2)->point(),
                       neighbor->vertex(3)->point());
        dist_[edge_count] = CGAL::squared_distance(cen1,cen2);
    }
    else
        dist_[edge_count] = -1.0;
}
void graphExporter::facetVertexFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count){

//    fbvf << "fb_vertex_inside_count fb_vertex_inside_dist_min fb_vertex_inside_dist_max fb_vertex_inside_dist_sum ";
//    fbvf << "fb_vertex_outside_count fb_vertex_outside_dist_min fb_vertex_outside_dist_max fb_vertex_outside_dist_sum" << endl;
//    fbvf << "fb_vertex_last_count fb_vertex_last_dist_min fb_vertex_last_dist_max fb_vertex_last_dist_sum" << endl;

    // inside
    int vdi = fci->info().fb_vertex_inside[j].size();
    fb_vertex_inside_count_[edge_count] = vdi;
    if(vdi > 0){
        fb_vertex_inside_dist_min_[edge_count] = *min_element(fci->info().fb_vertex_inside[j].begin(), fci->info().fb_vertex_inside[j].end());
        fb_vertex_inside_dist_max_[edge_count] = *max_element(fci->info().fb_vertex_inside[j].begin(), fci->info().fb_vertex_inside[j].end());
        fb_vertex_inside_dist_sum_[edge_count] = accumulate(fci->info().fb_vertex_inside[j].begin(), fci->info().fb_vertex_inside[j].end(),0.0);

    }
    // outside
    int vdo = fci->info().fb_vertex_outside[j].size();
    fb_vertex_outside_count_[edge_count] = vdo;
    if(vdo > 0){
        fb_vertex_outside_dist_min_[edge_count] = *min_element(fci->info().fb_vertex_outside[j].begin(), fci->info().fb_vertex_outside[j].end());
        fb_vertex_outside_dist_max_[edge_count] = *max_element(fci->info().fb_vertex_outside[j].begin(), fci->info().fb_vertex_outside[j].end());
        fb_vertex_outside_dist_sum_[edge_count] = accumulate(fci->info().fb_vertex_outside[j].begin(), fci->info().fb_vertex_outside[j].end(),0.0);
    }
    // last
    int vdl = fci->info().fb_vertex_last[j].size();
    fb_vertex_last_count_[edge_count] = vdl;
    if(vdl > 0){
        fb_vertex_last_dist_min_[edge_count] = *min_element(fci->info().fb_vertex_last[j].begin(), fci->info().fb_vertex_last[j].end());
        fb_vertex_last_dist_max_[edge_count] = *max_element(fci->info().fb_vertex_last[j].begin(), fci->info().fb_vertex_last[j].end());
        fb_vertex_last_dist_sum_[edge_count] = accumulate(fci->info().fb_vertex_last[j].begin(), fci->info().fb_vertex_last[j].end(),0.0);
    }
}

void graphExporter::facetFacetFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count){

//    fbff << "fb_facet_inside_count fb_facet_inside_dist_min fb_facet_inside_dist_max fb_facet_inside_dist_sum ";
//    fbff << "fb_facet_outside_count fb_facet_outside_dist_min fb_facet_outside_dist_max fb_facet_outside_dist_sum" << endl;
//    fbff << "fb_facet_last_count fb_facet_last_dist_min fb_facet_last_dist_max fb_facet_last_dist_sum" << endl;

    // inside
    int vdi = fci->info().fb_facet_inside[j].size();
    fb_facet_inside_count_[edge_count] = vdi;
    if(vdi > 0){
        fb_facet_inside_dist_min_[edge_count] = *min_element(fci->info().fb_facet_inside[j].begin(), fci->info().fb_facet_inside[j].end());
        fb_facet_inside_dist_max_[edge_count] = *max_element(fci->info().fb_facet_inside[j].begin(), fci->info().fb_facet_inside[j].end());
        fb_facet_inside_dist_sum_[edge_count] = accumulate(fci->info().fb_facet_inside[j].begin(), fci->info().fb_facet_inside[j].end(),0.0);
    }
    // outside
    int vdo = fci->info().fb_facet_outside[j].size();
    fb_facet_outside_count_[edge_count] = vdo;
    if(vdo > 0){
        fb_facet_outside_dist_min_[edge_count] = *min_element(fci->info().fb_facet_outside[j].begin(), fci->info().fb_facet_outside[j].end());
        fb_facet_outside_dist_max_[edge_count] = *max_element(fci->info().fb_facet_outside[j].begin(), fci->info().fb_facet_outside[j].end());
        fb_facet_outside_dist_sum_[edge_count] = accumulate(fci->info().fb_facet_outside[j].begin(), fci->info().fb_facet_outside[j].end(),0.0);
    }
    // last
    int vdl = fci->info().fb_facet_last[j].size();
    fb_facet_last_count_[edge_count] = vdl;
    if(vdl > 0){
        fb_facet_last_dist_min_[edge_count] = *min_element(fci->info().fb_facet_last[j].begin(), fci->info().fb_facet_last[j].end());
        fb_facet_last_dist_max_[edge_count] = *max_element(fci->info().fb_facet_last[j].begin(), fci->info().fb_facet_last[j].end());
        fb_facet_last_dist_sum_[edge_count] = accumulate(fci->info().fb_facet_last[j].begin(), fci->info().fb_facet_last[j].end(),0.0);
    }
}





void graphExporter::run(bool compress){

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nExport graph..." << endl;
    cout << "\t-to " << "dgnn/"+outfile_+"_X.npz" << endl;
//    string rays = (options_.export_rays) ? "\t-with rays" : "\t-without rays";
//    cout << rays << endl;


    /////////////////////////////////////////////
    //////////////////// CELLS //////////////////
    /////////////////////////////////////////////
    ////// info file
    ofstream info;
    info.open(dir_.path+"dgnn/"+outfile_+"_info.txt");
    info << "vertices " << nv_ << endl;
    info << "facets " << nf_ << endl;
    info << "cells " << nc_ << endl;
    info.close();


    int j, i, edge_count=0;


    Delaunay::All_cells_iterator fci;
    for(fci = Dt_.all_cells_begin(); fci != Dt_.all_cells_end(); fci++){

        i = fci->info().global_idx;

        /////////////////////////////////////////////
        //////////////     LABELS    ////////////////
        /////////////////////////////////////////////
        cellLabel(fci);

        /////////////////////////////////////////////
        ////////////// CELL FEATURES ////////////////
        /////////////////////////////////////////////
        cellGeometricFeatures(fci);
        cellVertexFeatures(fci);
        cellFacetFeatures(fci);



        /////////////////////////////////////////////
        ////////////// EDGE FEATURES ////////////////
        /////////////////////////////////////////////
        for(int c = 0; c < 4; c++){
            // set the edge to the mirror cell
            j = fci->neighbor(c)->info().global_idx;

            if(j == i) // skip self loops
                continue;

            adjacencies_(edge_count,0) = i;
            adjacencies_(edge_count,1) = j;

//            adjacency_ij << i << " ";
//            adjacency_ji << j << " ";

            //// geometric features
            facetGeometricFeatures(fci,c,edge_count);
            facetVertexFeatures(fci,c,edge_count);
            facetFacetFeatures(fci,c,edge_count);

            edge_count++;
        }
    }

    string ilf = dir_.path+"dgnn/"+outfile_+"_labels.npz";
    xt::dump_npz(ilf,"infinite",infinite_,compress,false);
    if(options_.ground_truth){
        xt::dump_npz(ilf,"inside_perc",inside_perc_,compress,true);
        xt::dump_npz(ilf,"outside_perc",outside_perc_,compress,true);
    }

    string cgeom = dir_.path+"dgnn/"+outfile_+"_cgeom.npz";
    xt::dump_npz(cgeom,"radius",radius_,compress,false);
    xt::dump_npz(cgeom,"vol",vol_,compress,true);
    xt::dump_npz(cgeom,"longest_edge",longest_edge_,compress,true);
    xt::dump_npz(cgeom,"shortest_edge",shortest_edge_,compress,true);

    ////// cell based vertex features
    string cvff = dir_.path+"dgnn/"+outfile_+"_cbvf.npz";
    xt::dump_npz(cvff,"cb_vertex_inside_count",cb_vertex_inside_count_,compress,false);
    xt::dump_npz(cvff,"cb_vertex_inside_dist_min",cb_vertex_inside_dist_min_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_inside_dist_max",cb_vertex_inside_dist_max_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_inside_dist_sum",cb_vertex_inside_dist_sum_,compress,true);

    xt::dump_npz(cvff,"cb_vertex_outside_count",cb_vertex_outside_count_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_outside_dist_min",cb_vertex_outside_dist_min_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_outside_dist_max",cb_vertex_outside_dist_max_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_outside_dist_sum",cb_vertex_outside_dist_sum_,compress,true);

    xt::dump_npz(cvff,"cb_vertex_last_count",cb_vertex_last_count_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_last_dist_min",cb_vertex_last_dist_min_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_last_dist_max",cb_vertex_last_dist_max_,compress,true);
    xt::dump_npz(cvff,"cb_vertex_last_dist_sum",cb_vertex_last_dist_sum_,compress,true);

    ////// cell based facet features
    string cbff = dir_.path+"dgnn/"+outfile_+"_cbff.npz";
    xt::dump_npz(cbff,"cb_facet_inside_first_count",cb_facet_inside_first_count_,compress,false);
    xt::dump_npz(cbff,"cb_facet_inside_first_dist_min",cb_facet_inside_first_dist_min_,compress,true);
    xt::dump_npz(cbff,"cb_facet_inside_first_dist_max",cb_facet_inside_first_dist_max_,compress,true);
    xt::dump_npz(cbff,"cb_facet_inside_first_dist_sum",cb_facet_inside_first_dist_sum_,compress,true);

    xt::dump_npz(cbff,"cb_facet_inside_second_count",cb_facet_inside_second_count_,compress,true);
    xt::dump_npz(cbff,"cb_facet_inside_second_dist_min",cb_facet_inside_second_dist_min_,compress,true);
    xt::dump_npz(cbff,"cb_facet_inside_second_dist_max",cb_facet_inside_second_dist_max_,compress,true);
    xt::dump_npz(cbff,"cb_facet_inside_second_dist_sum",cb_facet_inside_second_dist_sum_,compress,true);

    xt::dump_npz(cbff,"cb_facet_outside_first_count",cb_facet_outside_first_count_,compress,true);
    xt::dump_npz(cbff,"cb_facet_outside_first_dist_min",cb_facet_outside_first_dist_min_,compress,true);
    xt::dump_npz(cbff,"cb_facet_outside_first_dist_max",cb_facet_outside_first_dist_max_,compress,true);
    xt::dump_npz(cbff,"cb_facet_outside_first_dist_sum",cb_facet_outside_first_dist_sum_,compress,true);

    xt::dump_npz(cbff,"cb_facet_outside_second_count",cb_facet_outside_second_count_,compress,true);
    xt::dump_npz(cbff,"cb_facet_outside_second_dist_min",cb_facet_outside_second_dist_min_,compress,true);
    xt::dump_npz(cbff,"cb_facet_outside_second_dist_max",cb_facet_outside_second_dist_max_,compress,true);
    xt::dump_npz(cbff,"cb_facet_outside_second_dist_sum",cb_facet_outside_second_dist_sum_,compress,true);

    xt::dump_npz(cbff,"cb_facet_last_first_count",cb_facet_last_first_count_,compress,true);
    xt::dump_npz(cbff,"cb_facet_last_first_dist_min",cb_facet_last_first_dist_min_,compress,true);
    xt::dump_npz(cbff,"cb_facet_last_first_dist_max",cb_facet_last_first_dist_max_,compress,true);
    xt::dump_npz(cbff,"cb_facet_last_first_dist_sum",cb_facet_last_first_dist_sum_,compress,true);

    xt::dump_npz(cbff,"cb_facet_last_second_count",cb_facet_last_second_count_,compress,true);
    xt::dump_npz(cbff,"cb_facet_last_second_dist_min",cb_facet_last_second_dist_min_,compress,true);
    xt::dump_npz(cbff,"cb_facet_last_second_dist_max",cb_facet_last_second_dist_max_,compress,true);
    xt::dump_npz(cbff,"cb_facet_last_second_dist_sum",cb_facet_last_second_dist_sum_,compress,true);


    /////////////////////////////////////////////
    //////////////////// EDGES //////////////////
    /////////////////////////////////////////////
    string adjf = dir_.path+"dgnn/"+outfile_+"_adjacencies.npz";
    xt::dump_npz(adjf,"adjacencies",adjacencies_,compress,false);


    // geometric features
    string fgeom = dir_.path+"dgnn/"+outfile_+"_fgeom.npz";
    xt::dump_npz(fgeom,"area",area_,compress,false);
    xt::dump_npz(fgeom,"angle",angle_,compress,true);
    xt::dump_npz(fgeom,"beta",beta_,compress,true);
    xt::dump_npz(fgeom,"dist",dist_,compress,true);


    // facet based vertex features
    string fbvf = dir_.path+"dgnn/"+outfile_+"_fbvf.npz";
    xt::dump_npz(fbvf,"fb_vertex_inside_count",fb_vertex_inside_count_,compress,false);
    xt::dump_npz(fbvf,"fb_vertex_inside_dist_min",fb_vertex_inside_dist_min_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_inside_dist_max",fb_vertex_inside_dist_max_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_inside_dist_sum",fb_vertex_inside_dist_sum_,compress,true);

    xt::dump_npz(fbvf,"fb_vertex_outside_count",fb_vertex_outside_count_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_outside_dist_min",fb_vertex_outside_dist_min_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_outside_dist_max",fb_vertex_outside_dist_max_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_outside_dist_sum",fb_vertex_outside_dist_sum_,compress,true);

    xt::dump_npz(fbvf,"fb_vertex_last_count",fb_vertex_last_count_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_last_dist_min",fb_vertex_last_dist_min_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_last_dist_max",fb_vertex_last_dist_max_,compress,true);
    xt::dump_npz(fbvf,"fb_vertex_last_dist_sum",fb_vertex_last_dist_sum_,compress,true);


    // facet based facet features
    string fbff = dir_.path+"dgnn/"+outfile_+"_fbff.npz";
    xt::dump_npz(fbff,"fb_facet_inside_count",fb_facet_inside_count_,compress,false);
    xt::dump_npz(fbff,"fb_facet_inside_dist_min",fb_facet_inside_dist_min_,compress,true);
    xt::dump_npz(fbff,"fb_facet_inside_dist_max",fb_facet_inside_dist_max_,compress,true);
    xt::dump_npz(fbff,"fb_facet_inside_dist_sum",fb_facet_inside_dist_sum_,compress,true);

    xt::dump_npz(fbff,"fb_facet_outside_count",fb_facet_outside_count_,compress,true);
    xt::dump_npz(fbff,"fb_facet_outside_dist_min",fb_facet_outside_dist_min_,compress,true);
    xt::dump_npz(fbff,"fb_facet_outside_dist_max",fb_facet_outside_dist_max_,compress,true);
    xt::dump_npz(fbff,"fb_facet_outside_dist_sum",fb_facet_outside_dist_sum_,compress,true);

    xt::dump_npz(fbff,"fb_facet_last_count",fb_facet_last_count_,compress,true);
    xt::dump_npz(fbff,"fb_facet_last_dist_min",fb_facet_last_dist_min_,compress,true);
    xt::dump_npz(fbff,"fb_facet_last_dist_max",fb_facet_last_dist_max_,compress,true);
    xt::dump_npz(fbff,"fb_facet_last_dist_sum",fb_facet_last_dist_sum_,compress,true);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-Exported " << i+1 << " nodes and " << edge_count << " edges" << endl;
    cout << "\t-in " << duration.count() << "s" << endl;
}


bool importOccPoints(dirHolder dir, dataHolder& data){

    fs::path path = fs::path(dir.path) / fs::path(dir.occ_file).replace_extension(".npz");
    path = path.lexically_normal();
    std::ifstream file(path.string());

    cout << "\nLoad occupancy points..." << endl;
    cout << "\t-from " << path.string() << endl;

    // load the occupancy points here and put them in a new vector
    // called data.gt_points

    if(!file){
        cout << "\nFILE DOES NOT EXIST OR IS EMPTY!" << endl;
        return 1;
    }
    auto a = xt::load_npz(path.string());
    if(a.find("points") == a.end()){
        cout << "\nERROR: No points found in NPZ file!" << endl;
        return 1;
    }
    if(a.find("occupancies") == a.end()){
        cout << "\nERROR: No occupancies found in NPZ file!" << endl;
        return 1;
    }

    data.xgt_points = a["points"].cast<float>();
    data.xgt_occupancies = a["occupancies"].cast<uint8_t>();

//    data.gt_points.clear();
//    data.gt_infos.clear();
//    Point p;
//    for(int i = 0; i < points.shape()[0]; i++){
//        p = Point(points(i,0),points(i,1),points(i,2));
//        data.gt_points.push_back(p);
//    }


}
bool point2TetraIndex(dataHolder& data){

    Point p;
    Cell_handle c;
    data.xgt_point2tet.clear();
//    assert(data.xgt_point2tet.size()==0);
    for(int i = 0; i < data.xgt_points.shape()[0]; i++){
        p = Point(data.xgt_points(i,0),data.xgt_points(i,1),data.xgt_points(i,2));
        c = data.Dt.locate(p);
        data.xgt_point2tet.push_back(c->info().global_idx);
    }

};
bool exportOccPoints(dirHolder dir, dataHolder& data){

    // make my own
    fs::path p(dir.write_file);
    string outfile = p.stem().string()+"_eval";

    fs::path path = fs::path(dir.path) / fs::path("dgnn") / fs::path(outfile).replace_extension(".npz");
    path = path.lexically_normal();

    cout << "\nExport GT points..." << endl;
//    cout << "\t-to " << path.string() << endl;
    cout << "\t-to " << "dgnn/"+outfile+".npz" << endl;

    bool compress=true;

    xt::dump_npz(path.string(),"points",data.xgt_points,compress,false);
    xt::dump_npz(path.string(),"occupancies",data.xgt_occupancies,compress,true);
    xt::dump_npz(path.string(),"tetIndex",xt::adapt(data.xgt_point2tet),compress,true);

    return 0;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// IMPORT ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool importPrediction(dirHolder dir, dataHolder& data, runningOptions options){


    fs::path path = fs::path(dir.path) / fs::path(dir.prediction_file).replace_extension(".npz");
    path = path.lexically_normal();
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
    if(options.occupancy_type == "lo"){
        pred = npz_map["logits"].cast<float>();
    }
    else if(options.occupancy_type == "so"){
        pred = npz_map["softmax"].cast<float>();
    }
    else if(options.occupancy_type == "si"){
        pred = npz_map["sigmoid"].cast<float>();
    }
    else{
        cout << options.occupancy_type << " is not a valid occupancy type, choose either lo or sm" << endl;
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



