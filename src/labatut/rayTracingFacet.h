#ifndef RAYTRACINGFACET_H
#define RAYTRACINGFACET_H

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <util/geometricOperations.h>

using namespace std;


//// partially taken from COLMAP c++ code


namespace processing{

class RayCaster{

    private:          // Access specifier
//        double alpha = 1;  // Attribute
        dataHolder& data_;
        runningOptions options_;

        int set_sink_weight_;

        // edge weight computer
        double visibility_threshold_;
        double visibility_normalization_;
        double distance_sigma_;
        double distance_threshold_;
        double distance_normalization_;

        // The standard deviation wrt. the number of images seen by each point.
        // Increasing this value decreases the influence of points seen in few images.
        double visibility_sigma=3.0;
        // The factor that is applied to the computed distance sigma, which is
        // automatically computed as the 25th percentile of edge lengths. A higher
        // value will increase the smoothness of the surface.
        double distance_sigma_factor=1.0;

        // for Labatu 2009
        double alpha_vis;
        double sigma; // should be adjusted to noise


    public:
        RayCaster(dataHolder& data, runningOptions options):
            data_(data),
            options_(options),
            sigma(options.labatut_sigma),
            alpha_vis(options.labatut_alpha),
            visibility_threshold_(5 * visibility_sigma),
            set_sink_weight_(1),
            visibility_normalization_(-0.5 / (visibility_sigma * visibility_sigma)){
                std::vector<float> edge_lengths;
                edge_lengths.reserve(data_.Dt.number_of_finite_edges());

                for (auto it = data_.Dt.finite_edges_begin(); it != data_.Dt.finite_edges_end(); ++it){
                edge_lengths.push_back((it->first->vertex(it->second)->point() -
                                       it->first->vertex(it->third)->point()).squared_length());
                }

                distance_sigma_ = distance_sigma_factor *
                               std::max(std::sqrt(Percentile(edge_lengths, 25)), 1e-7f);
                distance_threshold_ = 5 * distance_sigma_;
                distance_normalization_ = -0.5 / (distance_sigma_ * distance_sigma_);

                // make the global index of the cells
                int idx=0;
                for(auto cft = data.Dt.all_cells_begin(); cft!=data.Dt.all_cells_end(); cft++)
                    cft->info().global_idx = idx++;

            }

        void run();

        void traverseOutside(Vertex_handle vit, Cell_handle& current_cell,  int oppositeVertex, int sensor_idx, double tdist2);
        void traverseInside(Vertex_handle vit, Cell_handle& current_cell,  int oppositeVertex, int sensor_idx);

        void outside(Vertex_handle vit, int sensor_idx);
        void inside(Vertex_handle vit, int sensor_idx);

        double DistanceSigma() const { return distance_sigma_; }

        double ComputeVisibilityProb(const double distance_squared) const {

            return alpha_vis*(1-std::exp(-distance_squared/(2*sigma*sigma)));
        }

        // TODO: adjust this to the parameters given in experimental results in Labatu 2009

        double ComputeDistanceProb(const double distance_squared) const {

            return alpha_vis*(1-std::exp(-distance_squared/(2*sigma*sigma)));
        }

//        double ComputeVisibilityProb(const double visibility_squared) const {
////          if (visibility_squared < visibility_threshold_) {
////            return std::max(
////                0.0, 1.0 - std::exp(visibility_squared * visibility_normalization_));
////          } else {
//            return 1.0;
////          }
//        }

//        // TODO: adjust this to the parameters given in experimental results in Labatu 2009

//        double ComputeDistanceProb(const double distance_squared) const {
////          if (distance_squared < distance_threshold_) {
////            return std::max(
////                0.0, 1.0 - std::exp(distance_squared * distance_normalization_));
////          } else {
//            return 1.0;
////          }
//        }
};


}



#endif // RAYTRACINGFACET_H
