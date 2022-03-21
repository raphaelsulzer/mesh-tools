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

    typedef float xtype;

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
    xt::xarray<xtype> inside_perc_;
    xt::xarray<xtype> outside_perc_;
    xt::xarray<int> infinite_;

    xt::xarray<xtype> radius_;
    xt::xarray<xtype> vol_;
    xt::xarray<xtype> longest_edge_;
    xt::xarray<xtype> shortest_edge_;

    ////// cell based vertex features
    xt::xarray<xtype> cb_vertex_inside_count_;
    xt::xarray<xtype> cb_vertex_inside_dist_min_;
    xt::xarray<xtype> cb_vertex_inside_dist_max_;
    xt::xarray<xtype> cb_vertex_inside_dist_sum_;

    xt::xarray<xtype> cb_vertex_outside_count_;
    xt::xarray<xtype> cb_vertex_outside_dist_min_;
    xt::xarray<xtype> cb_vertex_outside_dist_max_;
    xt::xarray<xtype> cb_vertex_outside_dist_sum_;

    xt::xarray<xtype> cb_vertex_last_count_;
    xt::xarray<xtype> cb_vertex_last_dist_min_;
    xt::xarray<xtype> cb_vertex_last_dist_max_;
    xt::xarray<xtype> cb_vertex_last_dist_sum_;

    ////// cell based facet features
    xt::xarray<xtype> cb_facet_inside_first_count_ ;
    xt::xarray<xtype> cb_facet_inside_first_dist_min_;
    xt::xarray<xtype> cb_facet_inside_first_dist_max_;
    xt::xarray<xtype> cb_facet_inside_first_dist_sum_;

    xt::xarray<xtype> cb_facet_inside_second_count_;
    xt::xarray<xtype> cb_facet_inside_second_dist_min_;
    xt::xarray<xtype> cb_facet_inside_second_dist_max_;
    xt::xarray<xtype> cb_facet_inside_second_dist_sum_;

    xt::xarray<xtype> cb_facet_outside_first_count_;
    xt::xarray<xtype> cb_facet_outside_first_dist_min_;
    xt::xarray<xtype> cb_facet_outside_first_dist_max_;
    xt::xarray<xtype> cb_facet_outside_first_dist_sum_ ;

    xt::xarray<xtype> cb_facet_outside_second_count_;
    xt::xarray<xtype> cb_facet_outside_second_dist_min_;
    xt::xarray<xtype> cb_facet_outside_second_dist_max_;
    xt::xarray<xtype> cb_facet_outside_second_dist_sum_;

    xt::xarray<xtype> cb_facet_last_first_count_;
    xt::xarray<xtype> cb_facet_last_first_dist_min_;
    xt::xarray<xtype> cb_facet_last_first_dist_max_;
    xt::xarray<xtype> cb_facet_last_first_dist_sum_;

    xt::xarray<xtype> cb_facet_last_second_count_;
    xt::xarray<xtype> cb_facet_last_second_dist_min_;
    xt::xarray<xtype> cb_facet_last_second_dist_max_;
    xt::xarray<xtype> cb_facet_last_second_dist_sum_;


    /////////////////////////////////////////////
    //////////////////// EDGES //////////////////
    /////////////////////////////////////////////
    xt::xarray<int> adjacencies_;

    // geometric features
    xt::xarray<xtype> area_;
    xt::xarray<xtype> angle_;
    xt::xarray<xtype> beta_;
    xt::xarray<xtype> dist_;

    // facet based vertex features
    xt::xarray<xtype> fb_vertex_inside_count_;
    xt::xarray<xtype> fb_vertex_inside_dist_min_;
    xt::xarray<xtype> fb_vertex_inside_dist_max_;
    xt::xarray<xtype> fb_vertex_inside_dist_sum_;

    xt::xarray<xtype> fb_vertex_outside_count_;
    xt::xarray<xtype> fb_vertex_outside_dist_min_;
    xt::xarray<xtype> fb_vertex_outside_dist_max_;
    xt::xarray<xtype> fb_vertex_outside_dist_sum_;

    xt::xarray<xtype> fb_vertex_last_count_;
    xt::xarray<xtype> fb_vertex_last_dist_min_;
    xt::xarray<xtype> fb_vertex_last_dist_max_;
    xt::xarray<xtype> fb_vertex_last_dist_sum_;

    // facet based facet features
    xt::xarray<xtype> fb_facet_inside_count_;
    xt::xarray<xtype> fb_facet_inside_dist_min_;
    xt::xarray<xtype> fb_facet_inside_dist_max_;
    xt::xarray<xtype> fb_facet_inside_dist_sum_;

    xt::xarray<xtype> fb_facet_outside_count_;
    xt::xarray<xtype> fb_facet_outside_dist_min_;
    xt::xarray<xtype> fb_facet_outside_dist_max_;
    xt::xarray<xtype> fb_facet_outside_dist_sum_;

    xt::xarray<xtype> fb_facet_last_count_;
    xt::xarray<xtype> fb_facet_last_dist_min_;
    xt::xarray<xtype> fb_facet_last_dist_max_;
    xt::xarray<xtype> fb_facet_last_dist_sum_;



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
        inside_perc_ = xt::xarray<xtype>::from_shape({nc_});
        outside_perc_ = xt::xarray<xtype>::from_shape({nc_});

        radius_ = xt::zeros<xtype>({nc_});
        vol_ = xt::zeros<xtype>({nc_});
        longest_edge_ = xt::zeros<xtype>({nc_});
        shortest_edge_ = xt::zeros<xtype>({nc_});

        ////// cell based vertex features
        cb_vertex_inside_count_ = xt::zeros<xtype>({nc_});
        cb_vertex_inside_dist_min_ = xt::zeros<xtype>({nc_});
        cb_vertex_inside_dist_max_ = xt::zeros<xtype>({nc_});
        cb_vertex_inside_dist_sum_ = xt::zeros<xtype>({nc_});

        cb_vertex_outside_count_ = xt::zeros<xtype>({nc_});
        cb_vertex_outside_dist_min_ = xt::zeros<xtype>({nc_});
        cb_vertex_outside_dist_max_ = xt::zeros<xtype>({nc_});
        cb_vertex_outside_dist_sum_ = xt::zeros<xtype>({nc_});

        cb_vertex_last_count_ = xt::zeros<xtype>({nc_});
        cb_vertex_last_dist_min_ = xt::zeros<xtype>({nc_});
        cb_vertex_last_dist_max_ = xt::zeros<xtype>({nc_});
        cb_vertex_last_dist_sum_ = xt::zeros<xtype>({nc_});

        ////// cell based facet features
        cb_facet_inside_first_count_  = xt::zeros<xtype>({nc_});
        cb_facet_inside_first_dist_min_ = xt::zeros<xtype>({nc_});
        cb_facet_inside_first_dist_max_ = xt::zeros<xtype>({nc_});
        cb_facet_inside_first_dist_sum_ = xt::zeros<xtype>({nc_});

        cb_facet_inside_second_count_ = xt::zeros<xtype>({nc_});
        cb_facet_inside_second_dist_min_ = xt::zeros<xtype>({nc_});
        cb_facet_inside_second_dist_max_ = xt::zeros<xtype>({nc_});
        cb_facet_inside_second_dist_sum_ = xt::zeros<xtype>({nc_});

        cb_facet_outside_first_count_ = xt::zeros<xtype>({nc_});
        cb_facet_outside_first_dist_min_ = xt::zeros<xtype>({nc_});
        cb_facet_outside_first_dist_max_ = xt::zeros<xtype>({nc_});
        cb_facet_outside_first_dist_sum_  = xt::zeros<xtype>({nc_});

        cb_facet_outside_second_count_ = xt::zeros<xtype>({nc_});
        cb_facet_outside_second_dist_min_ = xt::zeros<xtype>({nc_});
        cb_facet_outside_second_dist_max_ = xt::zeros<xtype>({nc_});
        cb_facet_outside_second_dist_sum_ = xt::zeros<xtype>({nc_});

        cb_facet_last_first_count_ = xt::zeros<xtype>({nc_});
        cb_facet_last_first_dist_min_ = xt::zeros<xtype>({nc_});
        cb_facet_last_first_dist_max_ = xt::zeros<xtype>({nc_});
        cb_facet_last_first_dist_sum_ = xt::zeros<xtype>({nc_});

        cb_facet_last_second_count_ = xt::zeros<xtype>({nc_});
        cb_facet_last_second_dist_min_ = xt::zeros<xtype>({nc_});
        cb_facet_last_second_dist_max_ = xt::zeros<xtype>({nc_});
        cb_facet_last_second_dist_sum_ = xt::zeros<xtype>({nc_});


        /////////////////////////////////////////////
        //////////////////// EDGES //////////////////
        /////////////////////////////////////////////
        adjacencies_ = xt::xarray<int>::from_shape({nf_*2,2});

        // geometric facet based features
        area_ = xt::zeros<xtype>({nf_*2});
        angle_ = xt::zeros<xtype>({nf_*2});
        beta_ = xt::zeros<xtype>({nf_*2});
        dist_ = xt::zeros<xtype>({nf_*2});

        // facet based vertex features
        fb_vertex_inside_count_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_inside_dist_min_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_inside_dist_max_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_inside_dist_sum_ = xt::zeros<xtype>({nf_*2});

        fb_vertex_outside_count_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_outside_dist_min_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_outside_dist_max_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_outside_dist_sum_ = xt::zeros<xtype>({nf_*2});

        fb_vertex_last_count_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_last_dist_min_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_last_dist_max_ = xt::zeros<xtype>({nf_*2});
        fb_vertex_last_dist_sum_ = xt::zeros<xtype>({nf_*2});

        // facet based facet features
        fb_facet_inside_count_ = xt::zeros<xtype>({nf_*2});
        fb_facet_inside_dist_min_ = xt::zeros<xtype>({nf_*2});
        fb_facet_inside_dist_max_ = xt::zeros<xtype>({nf_*2});
        fb_facet_inside_dist_sum_ = xt::zeros<xtype>({nf_*2});

        fb_facet_outside_count_ = xt::zeros<xtype>({nf_*2});
        fb_facet_outside_dist_min_ = xt::zeros<xtype>({nf_*2});
        fb_facet_outside_dist_max_ = xt::zeros<xtype>({nf_*2});
        fb_facet_outside_dist_sum_ = xt::zeros<xtype>({nf_*2});

        fb_facet_last_count_ = xt::zeros<xtype>({nf_*2});
        fb_facet_last_dist_min_ = xt::zeros<xtype>({nf_*2});
        fb_facet_last_dist_max_ = xt::zeros<xtype>({nf_*2});
        fb_facet_last_dist_sum_ = xt::zeros<xtype>({nf_*2});

    }

    void run(bool compress);



    void cellLabel(Delaunay::All_cells_iterator& fci);

    void cellGeometricFeatures(Delaunay::All_cells_iterator& fci);
    void cellVertexFeatures(Delaunay::All_cells_iterator& fci);
    void cellFacetFeatures(Delaunay::All_cells_iterator& fci);

    void facetGeometricFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count);
    void facetVertexFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count);
    void facetFacetFeatures(Delaunay::All_cells_iterator& fci, int j, int edge_count);

};

// import functions
bool importPrediction(dirHolder dir, dataHolder& data, runningOptions options);
bool importOccPoints(dirHolder dir, dataHolder& data);
bool point2TetraIndex(dataHolder& data);
bool exportOccPoints(dirHolder dir, dataHolder& data);


#endif // LEARNINGIO_H
