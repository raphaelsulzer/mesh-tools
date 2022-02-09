#ifndef LEARNINGIO_H
#define LEARNINGIO_H

#include <base/cgal_typedefs.h>
#include <learning/learning.h>

#include <xtensor/xbuilder.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor-io/xnpz.hpp>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
using namespace std;


class graphExporter{

private:
    size_t nc_;
    size_t nf_;
    size_t nv_;

    dirHolder& dir_;
    runningOptions& options_;
    Delaunay& Dt_;

    string outfile_;

    /////////////////////////////////////////////
    //////////////////// CELLS //////////////////
    /////////////////////////////////////////////
    // cell index and labels
    xt::xarray<double> inside_perc_;
    xt::xarray<double> outside_perc_;
    xt::xarray<int> infinite_;

    xt::xarray<double> radius_;
    xt::xarray<double> vol_;
    xt::xarray<double> longest_edge_;
    xt::xarray<double> shortest_edge_;

    ////// cell based vertex features
    xt::xarray<double> cb_vertex_inside_count_;
    xt::xarray<double> cb_vertex_inside_dist_min_;
    xt::xarray<double> cb_vertex_inside_dist_max_;
    xt::xarray<double> cb_vertex_inside_dist_sum_;

    xt::xarray<double> cb_vertex_outside_count_;
    xt::xarray<double> cb_vertex_outside_dist_min_;
    xt::xarray<double> cb_vertex_outside_dist_max_;
    xt::xarray<double> cb_vertex_outside_dist_sum_;

    xt::xarray<double> cb_vertex_last_count_;
    xt::xarray<double> cb_vertex_last_dist_min_;
    xt::xarray<double> cb_vertex_last_dist_max_;
    xt::xarray<double> cb_vertex_last_dist_sum_;

    ////// cell based facet features
    xt::xarray<double> cb_facet_inside_first_count_ ;
    xt::xarray<double> cb_facet_inside_first_dist_min_;
    xt::xarray<double> cb_facet_inside_first_dist_max_;
    xt::xarray<double> cb_facet_inside_first_dist_sum_;

    xt::xarray<double> cb_facet_inside_second_count_;
    xt::xarray<double> cb_facet_inside_second_dist_min_;
    xt::xarray<double> cb_facet_inside_second_dist_max_;
    xt::xarray<double> cb_facet_inside_second_dist_sum_;

    xt::xarray<double> cb_facet_outside_first_count_;
    xt::xarray<double> cb_facet_outside_first_dist_min_;
    xt::xarray<double> cb_facet_outside_first_dist_max_;
    xt::xarray<double> cb_facet_outside_first_dist_sum_ ;

    xt::xarray<double> cb_facet_outside_second_count_;
    xt::xarray<double> cb_facet_outside_second_dist_min_;
    xt::xarray<double> cb_facet_outside_second_dist_max_;
    xt::xarray<double> cb_facet_outside_second_dist_sum_;

    xt::xarray<double> cb_facet_last_first_count_;
    xt::xarray<double> cb_facet_last_first_dist_min_;
    xt::xarray<double> cb_facet_last_first_dist_max_;
    xt::xarray<double> cb_facet_last_first_dist_sum_;

    xt::xarray<double> cb_facet_last_second_count_;
    xt::xarray<double> cb_facet_last_second_dist_min_;
    xt::xarray<double> cb_facet_last_second_dist_max_;
    xt::xarray<double> cb_facet_last_second_dist_sum_;


    /////////////////////////////////////////////
    //////////////////// EDGES //////////////////
    /////////////////////////////////////////////
    xt::xarray<int> adjacencies_;

    // geometric features
    xt::xarray<double> area_;
    xt::xarray<double> angle_;
    xt::xarray<double> beta_;
    xt::xarray<double> dist_;

    // facet based vertex features
    xt::xarray<double> fb_vertex_inside_count_;
    xt::xarray<double> fb_vertex_inside_dist_min_;
    xt::xarray<double> fb_vertex_inside_dist_max_;
    xt::xarray<double> fb_vertex_inside_dist_sum_;

    xt::xarray<double> fb_vertex_outside_count_;
    xt::xarray<double> fb_vertex_outside_dist_min_;
    xt::xarray<double> fb_vertex_outside_dist_max_;
    xt::xarray<double> fb_vertex_outside_dist_sum_;

    xt::xarray<double> fb_vertex_last_count_;
    xt::xarray<double> fb_vertex_last_dist_min_;
    xt::xarray<double> fb_vertex_last_dist_max_;
    xt::xarray<double> fb_vertex_last_dist_sum_;

    // facet based facet features
    xt::xarray<double> fb_facet_inside_count_;
    xt::xarray<double> fb_facet_inside_dist_min_;
    xt::xarray<double> fb_facet_inside_dist_max_;
    xt::xarray<double> fb_facet_inside_dist_sum_;

    xt::xarray<double> fb_facet_outside_count_;
    xt::xarray<double> fb_facet_outside_dist_min_;
    xt::xarray<double> fb_facet_outside_dist_max_;
    xt::xarray<double> fb_facet_outside_dist_sum_;

    xt::xarray<double> fb_facet_last_count_;
    xt::xarray<double> fb_facet_last_dist_min_;
    xt::xarray<double> fb_facet_last_dist_max_;
    xt::xarray<double> fb_facet_last_dist_sum_;



public:
    graphExporter(dirHolder& dir, Delaunay& Dt, runningOptions& options):
        dir_(dir),
        Dt_(Dt),
        options_(options)
    {

        nv_ = Dt_.number_of_vertices();
        nf_ = Dt_.number_of_facets();
        nc_ = Dt_.number_of_cells();

        boost::filesystem::path p(dir.write_file);
        outfile_ = p.stem().string();

        infinite_ = xt::xarray<int>::from_shape({nc_});
        inside_perc_ = xt::xarray<double>::from_shape({nc_});
        outside_perc_ = xt::xarray<double>::from_shape({nc_});

        radius_ = xt::zeros<double>({nc_});
        vol_ = xt::zeros<double>({nc_});
        longest_edge_ = xt::zeros<double>({nc_});
        shortest_edge_ = xt::zeros<double>({nc_});

        ////// cell based vertex features
        cb_vertex_inside_count_ = xt::zeros<double>({nc_});
        cb_vertex_inside_dist_min_ = xt::zeros<double>({nc_});
        cb_vertex_inside_dist_max_ = xt::zeros<double>({nc_});
        cb_vertex_inside_dist_sum_ = xt::zeros<double>({nc_});

        cb_vertex_outside_count_ = xt::zeros<double>({nc_});
        cb_vertex_outside_dist_min_ = xt::zeros<double>({nc_});
        cb_vertex_outside_dist_max_ = xt::zeros<double>({nc_});
        cb_vertex_outside_dist_sum_ = xt::zeros<double>({nc_});

        cb_vertex_last_count_ = xt::zeros<double>({nc_});
        cb_vertex_last_dist_min_ = xt::zeros<double>({nc_});
        cb_vertex_last_dist_max_ = xt::zeros<double>({nc_});
        cb_vertex_last_dist_sum_ = xt::zeros<double>({nc_});

        ////// cell based facet features
        cb_facet_inside_first_count_  = xt::zeros<double>({nc_});
        cb_facet_inside_first_dist_min_ = xt::zeros<double>({nc_});
        cb_facet_inside_first_dist_max_ = xt::zeros<double>({nc_});
        cb_facet_inside_first_dist_sum_ = xt::zeros<double>({nc_});

        cb_facet_inside_second_count_ = xt::zeros<double>({nc_});
        cb_facet_inside_second_dist_min_ = xt::zeros<double>({nc_});
        cb_facet_inside_second_dist_max_ = xt::zeros<double>({nc_});
        cb_facet_inside_second_dist_sum_ = xt::zeros<double>({nc_});

        cb_facet_outside_first_count_ = xt::zeros<double>({nc_});
        cb_facet_outside_first_dist_min_ = xt::zeros<double>({nc_});
        cb_facet_outside_first_dist_max_ = xt::zeros<double>({nc_});
        cb_facet_outside_first_dist_sum_  = xt::zeros<double>({nc_});

        cb_facet_outside_second_count_ = xt::zeros<double>({nc_});
        cb_facet_outside_second_dist_min_ = xt::zeros<double>({nc_});
        cb_facet_outside_second_dist_max_ = xt::zeros<double>({nc_});
        cb_facet_outside_second_dist_sum_ = xt::zeros<double>({nc_});

        cb_facet_last_first_count_ = xt::zeros<double>({nc_});
        cb_facet_last_first_dist_min_ = xt::zeros<double>({nc_});
        cb_facet_last_first_dist_max_ = xt::zeros<double>({nc_});
        cb_facet_last_first_dist_sum_ = xt::zeros<double>({nc_});

        cb_facet_last_second_count_ = xt::zeros<double>({nc_});
        cb_facet_last_second_dist_min_ = xt::zeros<double>({nc_});
        cb_facet_last_second_dist_max_ = xt::zeros<double>({nc_});
        cb_facet_last_second_dist_sum_ = xt::zeros<double>({nc_});


        /////////////////////////////////////////////
        //////////////////// EDGES //////////////////
        /////////////////////////////////////////////
        adjacencies_ = xt::xarray<int>::from_shape({nf_*2,2});

        // geometric facet based features
        area_ = xt::zeros<double>({nf_*2});
        angle_ = xt::zeros<double>({nf_*2});
        beta_ = xt::zeros<double>({nf_*2});
        dist_ = xt::zeros<double>({nf_*2});

        // facet based vertex features
        fb_vertex_inside_count_ = xt::zeros<double>({nf_*2});
        fb_vertex_inside_dist_min_ = xt::zeros<double>({nf_*2});
        fb_vertex_inside_dist_max_ = xt::zeros<double>({nf_*2});
        fb_vertex_inside_dist_sum_ = xt::zeros<double>({nf_*2});

        fb_vertex_outside_count_ = xt::zeros<double>({nf_*2});
        fb_vertex_outside_dist_min_ = xt::zeros<double>({nf_*2});
        fb_vertex_outside_dist_max_ = xt::zeros<double>({nf_*2});
        fb_vertex_outside_dist_sum_ = xt::zeros<double>({nf_*2});

        fb_vertex_last_count_ = xt::zeros<double>({nf_*2});
        fb_vertex_last_dist_min_ = xt::zeros<double>({nf_*2});
        fb_vertex_last_dist_max_ = xt::zeros<double>({nf_*2});
        fb_vertex_last_dist_sum_ = xt::zeros<double>({nf_*2});

        // facet based facet features
        fb_facet_inside_count_ = xt::zeros<double>({nf_*2});
        fb_facet_inside_dist_min_ = xt::zeros<double>({nf_*2});
        fb_facet_inside_dist_max_ = xt::zeros<double>({nf_*2});
        fb_facet_inside_dist_sum_ = xt::zeros<double>({nf_*2});

        fb_facet_outside_count_ = xt::zeros<double>({nf_*2});
        fb_facet_outside_dist_min_ = xt::zeros<double>({nf_*2});
        fb_facet_outside_dist_max_ = xt::zeros<double>({nf_*2});
        fb_facet_outside_dist_sum_ = xt::zeros<double>({nf_*2});

        fb_facet_last_count_ = xt::zeros<double>({nf_*2});
        fb_facet_last_dist_min_ = xt::zeros<double>({nf_*2});
        fb_facet_last_dist_max_ = xt::zeros<double>({nf_*2});
        fb_facet_last_dist_sum_ = xt::zeros<double>({nf_*2});

    }

    void run(bool compress);
    void cellLabel(Delaunay::All_cells_iterator& fci);

    void cellGeometricFeatures(Delaunay::All_cells_iterator& fci);
    void cellVertexFeatures(Delaunay::All_cells_iterator& fci);
    void cellFacetFeatures(Delaunay::All_cells_iterator& fci);

    void facetGeometricFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count);
    void facetVertexFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count);
    void facetFacetFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count);






//    ////// geometric features
//    string cgeom = dir.path+"gt/"+outfile+"_cgeom.npz";
//    = xt::xarray<double> radius = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> vol = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> longest_edge = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> shortest_edge = xt::xarray<double>::from_shape({nc});

//    ////// cell based vertex features
//    string cbvf = dir.path+"gt/"+outfile+"_cbvf.npz";
//    = xt::xarray<double> cb_vertex_inside_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_inside_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_inside_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_inside_dist_sum = xt::xarray<double>::from_shape({nc});

//    = xt::xarray<double> cb_vertex_outside_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_outside_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_outside_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_outside_dist_sum = xt::xarray<double>::from_shape({nc});

//    = xt::xarray<double> cb_vertex_last_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_last_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_last_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_vertex_last_dist_sum = xt::xarray<double>::from_shape({nc});

//    ////// cell based facet features
//    string cbff = dir.path+"gt/"+outfile+"_cbff.npz";
//    = xt::xarray<double> cb_facet_inside_first_count  = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_inside_first_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_inside_first_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_inside_first_dist_sum = xt::xarray<double>::from_shape({nc});

//    = xt::xarray<double> cb_facet_inside_second_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_inside_second_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_inside_second_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_inside_second_dist_sum  = xt::xarray<double>::from_shape({nc});

//    = xt::xarray<double> cb_facet_outside_first_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_outside_first_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_outside_first_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_outside_first_dist_sum  = xt::xarray<double>::from_shape({nc});

//    = xt::xarray<double> cb_facet_outside_second_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_outside_second_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_outside_second_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_outside_second_dist_sum = xt::xarray<double>::from_shape({nc});

//    = xt::xarray<double> cb_facet_last_first_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_last_first_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_last_first_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_last_first_dist_sum = xt::xarray<double>::from_shape({nc});

//    = xt::xarray<double> cb_facet_last_second_count = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_last_second_dist_min = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_last_second_dist_max = xt::xarray<double>::from_shape({nc});
//    = xt::xarray<double> cb_facet_last_second_dist_sum = xt::xarray<double>::from_shape({nc});


//    /////////////////////////////////////////////
//    //////////////////// EDGES //////////////////
//    /////////////////////////////////////////////
//    // adjacency files
//    string adj = dir.path+"gt/"+outfile+"_adjacencies.npz";
//    = xt::xarray<double> adjacencies = xt::xarray<double>::from_shape({nf*2,2});

//    // geometric features
//    string fgeom = dir.path+"gt/"+outfile+"_fgeom.npz";
//    = xt::xarray<double> area = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> angle = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> beta = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> dist = xt::xarray<double>::from_shape({nf*2});

//    // facet based vertex features
//    string fbvf = dir.path+"gt/"+outfile+"_fbvf.npz";
//    = xt::xarray<double> fb_vertex_inside_count = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_inside_dist_min = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_inside_dist_max = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_inside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//    = xt::xarray<double> fb_vertex_outside_count = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_outside_dist_min = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_outside_dist_max = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_outside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//    = xt::xarray<double> fb_vertex_last_count = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_last_dist_min = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_last_dist_max = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_vertex_last_dist_sum = xt::xarray<double>::from_shape({nf*2});

//    // facet based facet features
//    string fbff = dir.path+"gt/"+outfile+"_fbff.npz";
//    = xt::xarray<double> fb_facet_inside_count = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_inside_dist_min = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_inside_dist_max = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_inside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//    = xt::xarray<double> fb_facet_outside_count = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_outside_dist_min = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_outside_dist_max = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_outside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//    = xt::xarray<double> fb_facet_last_count = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_last_dist_min = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_last_dist_max = xt::xarray<double>::from_shape({nf*2});
//    = xt::xarray<double> fb_facet_last_dist_sum = xt::xarray<double>::from_shape({nf*2});



};


//////// label and index
//// cell index and labels
//string il = dir.path+"gt/"+outfile+"_labels.npz";

//xt::xarray<double> inside_perc;

////    static const = xt::xarray<double> inside_perc = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> outside_perc = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> infinite = xt::xarray<double>::from_shape({nc});


//////// geometric features
//string cgeom = dir.path+"gt/"+outfile+"_cgeom.npz";
//= xt::xarray<double> radius = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> vol = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> longest_edge = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> shortest_edge = xt::xarray<double>::from_shape({nc});

//////// cell based vertex features
//string cbvf = dir.path+"gt/"+outfile+"_cbvf.npz";
//= xt::xarray<double> cb_vertex_inside_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_inside_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_inside_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_inside_dist_sum = xt::xarray<double>::from_shape({nc});

//= xt::xarray<double> cb_vertex_outside_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_outside_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_outside_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_outside_dist_sum = xt::xarray<double>::from_shape({nc});

//= xt::xarray<double> cb_vertex_last_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_last_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_last_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_vertex_last_dist_sum = xt::xarray<double>::from_shape({nc});

//////// cell based facet features
//string cbff = dir.path+"gt/"+outfile+"_cbff.npz";
//= xt::xarray<double> cb_facet_inside_first_count  = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_inside_first_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_inside_first_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_inside_first_dist_sum = xt::xarray<double>::from_shape({nc});

//= xt::xarray<double> cb_facet_inside_second_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_inside_second_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_inside_second_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_inside_second_dist_sum  = xt::xarray<double>::from_shape({nc});

//= xt::xarray<double> cb_facet_outside_first_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_outside_first_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_outside_first_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_outside_first_dist_sum  = xt::xarray<double>::from_shape({nc});

//= xt::xarray<double> cb_facet_outside_second_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_outside_second_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_outside_second_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_outside_second_dist_sum = xt::xarray<double>::from_shape({nc});

//= xt::xarray<double> cb_facet_last_first_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_last_first_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_last_first_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_last_first_dist_sum = xt::xarray<double>::from_shape({nc});

//= xt::xarray<double> cb_facet_last_second_count = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_last_second_dist_min = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_last_second_dist_max = xt::xarray<double>::from_shape({nc});
//= xt::xarray<double> cb_facet_last_second_dist_sum = xt::xarray<double>::from_shape({nc});


///////////////////////////////////////////////
////////////////////// EDGES //////////////////
///////////////////////////////////////////////
//// adjacency files
//string adj = dir.path+"gt/"+outfile+"_adjacencies.npz";
//= xt::xarray<double> adjacencies = xt::xarray<double>::from_shape({nf*2,2});

//// geometric features
//string fgeom = dir.path+"gt/"+outfile+"_fgeom.npz";
//= xt::xarray<double> area = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> angle = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> beta = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> dist = xt::xarray<double>::from_shape({nf*2});

//// facet based vertex features
//string fbvf = dir.path+"gt/"+outfile+"_fbvf.npz";
//= xt::xarray<double> fb_vertex_inside_count = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_inside_dist_min = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_inside_dist_max = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_inside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//= xt::xarray<double> fb_vertex_outside_count = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_outside_dist_min = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_outside_dist_max = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_outside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//= xt::xarray<double> fb_vertex_last_count = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_last_dist_min = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_last_dist_max = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_vertex_last_dist_sum = xt::xarray<double>::from_shape({nf*2});

//// facet based facet features
//string fbff = dir.path+"gt/"+outfile+"_fbff.npz";
//= xt::xarray<double> fb_facet_inside_count = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_inside_dist_min = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_inside_dist_max = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_inside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//= xt::xarray<double> fb_facet_outside_count = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_outside_dist_min = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_outside_dist_max = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_outside_dist_sum = xt::xarray<double>::from_shape({nf*2});

//= xt::xarray<double> fb_facet_last_count = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_last_dist_min = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_last_dist_max = xt::xarray<double>::from_shape({nf*2});
//= xt::xarray<double> fb_facet_last_dist_sum = xt::xarray<double>::from_shape({nf*2});

// import functions
int loadPrediction(dirHolder dir, dataHolder& data, runningOptions options);


#endif // LEARNINGIO_H
