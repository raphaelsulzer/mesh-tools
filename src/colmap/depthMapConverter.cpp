// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

//#include "mvs/fusion.h"
#include "IO/depthMapConverter.h"
#include "util/misc.h"

#include <base/cgal_typedefs.h>
#include <math/vectorArithmetic.h>


namespace colmap {
namespace mvs {

namespace myInternal {

template <typename T>
float Median(std::vector<T>* elems) {
  CHECK(!elems->empty());
  const size_t mid_idx = elems->size() / 2;
  std::nth_element(elems->begin(), elems->begin() + mid_idx, elems->end());
  if (elems->size() % 2 == 0) {
    const float mid_element1 = static_cast<float>((*elems)[mid_idx]);
    const float mid_element2 = static_cast<float>(
        *std::max_element(elems->begin(), elems->begin() + mid_idx));
    return (mid_element1 + mid_element2) / 2.0f;
  } else {
    return static_cast<float>((*elems)[mid_idx]);
  }
}

// Use the sparse model to find most connected image that has not yet been
// fused. This is used as a heuristic to ensure that the workspace cache reuses
// already cached images as efficient as possible.
int FindNextImage(const std::vector<std::vector<int>>& overlapping_images,
                  const std::vector<char>& used_images,
                  const std::vector<char>& fused_images,
                  const int prev_image_idx) {
  CHECK_EQ(used_images.size(), fused_images.size());

  for (const auto image_idx : overlapping_images.at(prev_image_idx)) {
    if (used_images.at(image_idx) && !fused_images.at(image_idx)) {
      return image_idx;
    }
  }

  // If none of the overlapping images are not yet fused, simply return the
  // first image that has not yet been fused.
  for (size_t image_idx = 0; image_idx < fused_images.size(); ++image_idx) {
    if (used_images[image_idx] && !fused_images[image_idx]) {
      return image_idx;
    }
  }

  return -1;
}

}  // namespace myInternal






DepthMapConverter::DepthMapConverter(const StereoFusionOptions& options,
                           const std::string& workspace_path,
                           const std::string& workspace_format,
                           const std::string& pmvs_option_name,
                           const std::string& input_type)
    : StereoFusion(options, workspace_path, workspace_format, pmvs_option_name, input_type){}



void DepthMapConverter::run(){

  fused_points_.clear();
  fused_points_visibility_.clear();

  options_.Print();
  std::cout << std::endl;

  std::cout << "Reading workspace..." << std::endl;

  Workspace::Options workspace_options;

  auto workspace_format_lower_case = workspace_format_;
  StringToLower(&workspace_format_lower_case);
  if (workspace_format_lower_case == "pmvs") {
    workspace_options.stereo_folder =
        StringPrintf("stereo-%s", pmvs_option_name_.c_str());
  }

  workspace_options.max_image_size = options_.max_image_size;
  workspace_options.image_as_rgb = true;
  workspace_options.cache_size = options_.cache_size;
  workspace_options.workspace_path = workspace_path_;
  workspace_options.workspace_format = workspace_format_;
  workspace_options.input_type = input_type_;

  workspace_.reset(new Workspace(workspace_options));

  if (IsStopped()) {
    GetTimer().PrintMinutes();
    return;
  }

  std::cout << "Reading configuration..." << std::endl;

  const auto& model = workspace_->GetModel();

  const double kMinTriangulationAngle = 0;
  if (model.GetMaxOverlappingImagesFromPMVS().empty()) {
    overlapping_images_ = model.GetMaxOverlappingImages(
        options_.check_num_images, kMinTriangulationAngle);
  } else {
    overlapping_images_ = model.GetMaxOverlappingImagesFromPMVS();
  }

  used_images_.resize(model.images.size(), false);
  fused_images_.resize(model.images.size(), false);
  fused_pixel_masks_.resize(model.images.size());
  depth_map_sizes_.resize(model.images.size());
  bitmap_scales_.resize(model.images.size());
  P_.resize(model.images.size());
  inv_P_.resize(model.images.size());
  inv_R_.resize(model.images.size());

  const auto image_names = ReadTextFileLines(JoinPaths(
      workspace_path_, workspace_options.stereo_folder, "fusion.cfg"));

  // iterate over every image in the model/workspace
  // this loop just inits the fused_pixel_masks, depth_maps, etc.
  // it is like a getDataFromModelFunction
  for (const auto& image_name : image_names){
    const int image_idx = model.GetImageIdx(image_name);

    if (!workspace_->HasBitmap(image_idx) ||
        !workspace_->HasDepthMap(image_idx) ||
        !workspace_->HasNormalMap(image_idx)) {
      std::cout
          << StringPrintf(
                 "WARNING: Ignoring image %s, because input does not exist.",
                 image_name.c_str())
          << std::endl;
      continue;
    }

    // get the image and the depth map
    const auto& image = model.images.at(image_idx);
    const auto& depth_map = workspace_->GetDepthMap(image_idx);

    // set image as already used
    used_images_.at(image_idx) = true;

    // init the fused_pixel_mask, which will later be used to the if a pixel has already been visited or not
    fused_pixel_masks_.at(image_idx) =
        Mat<bool>(depth_map.GetWidth(), depth_map.GetHeight(), 1);
    fused_pixel_masks_.at(image_idx).Fill(false);

    depth_map_sizes_.at(image_idx) =
        std::make_pair(depth_map.GetWidth(), depth_map.GetHeight());

    bitmap_scales_.at(image_idx) = std::make_pair(
        static_cast<float>(depth_map.GetWidth()) / image.GetWidth(),
        static_cast<float>(depth_map.GetHeight()) / image.GetHeight());

    Eigen::Matrix<float, 3, 3, Eigen::RowMajor> K =
        Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(
            image.GetK());
    K(0, 0) *= bitmap_scales_.at(image_idx).first;
    K(0, 2) *= bitmap_scales_.at(image_idx).first;
    K(1, 1) *= bitmap_scales_.at(image_idx).second;
    K(1, 2) *= bitmap_scales_.at(image_idx).second;


    ComposeProjectionMatrix(K.data(), image.GetR(), image.GetT(),
                            P_.at(image_idx).data());
    ComposeInverseProjectionMatrix(K.data(), image.GetR(), image.GetT(),
                                   inv_P_.at(image_idx).data());
    inv_R_.at(image_idx) =
        Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(
            image.GetR())
            .transpose();
  } // iterate over all images closed

  // now iterate over the images to actually fuse them
  // so this here can be seen as the main loop
  size_t num_fused_images = 0;
  for (int image_idx = 0; image_idx >= 0;
       // start with an image that has very big overlap to another one
       image_idx = myInternal::FindNextImage(overlapping_images_, used_images_,
                                           fused_images_, image_idx)) {
    if (IsStopped()) {
      break;
    }

    Timer timer;
    timer.Start();

    std::cout << StringPrintf("Fusing image [%d/%d]", num_fused_images + 1,
                              model.images.size())
              << std::flush;

    // get width, height and the fused_pixel_mask of this image
    const int width = depth_map_sizes_.at(image_idx).first;
    const int height = depth_map_sizes_.at(image_idx).second;
    const auto& fused_pixel_mask = fused_pixel_masks_.at(image_idx);

    // init a data holder for the current image
    FusionData data;
    data.image_idx = image_idx;
    data.traversal_depth = 0;

    // iterate over all the pixels of the current image
    for (data.row = 0; data.row < height; ++data.row) {
      for (data.col = 0; data.col < width; ++data.col) {
        if (fused_pixel_mask.Get(data.row, data.col)) {
          continue;
        }

        fusion_queue_.push_back(data);

        fuse();
      }
    }

    num_fused_images += 1;
    // set the current image as fused
    fused_images_.at(image_idx) = true;

    std::cout << StringPrintf(" in %.3fs (%d points)", timer.ElapsedSeconds(),
                              fused_points_.size())
              << std::endl;
  } // closing main image iteration loop

  fused_points_.shrink_to_fit();
  fused_points_visibility_.shrink_to_fit();

  if (fused_points_.empty()) {
    std::cout << "WARNING: Could not fuse any points. This is likely caused by "
                 "incorrect settings - filtering must be enabled for the last "
                 "call to patch match stereo."
              << std::endl;
  }

  std::cout << "Number of fused points: " << fused_points_.size() << std::endl;
  GetTimer().PrintMinutes();

}

void DepthMapConverter::fuse() {
//    void DepthMapConverter::fuse(std::vector<Point>& points, std::vector<vertex_info>& points_info, std::vector<std::vector<image_info>>) {
  CHECK_EQ(fusion_queue_.size(), 1);

  Eigen::Vector4f fused_ref_point = Eigen::Vector4f::Zero();
  Eigen::Vector3f fused_ref_normal = Eigen::Vector3f::Zero();

  fused_point_x_.clear();
  fused_point_y_.clear();
  fused_point_z_.clear();
  fused_point_nx_.clear();
  fused_point_ny_.clear();
  fused_point_nz_.clear();
  fused_point_r_.clear();
  fused_point_g_.clear();
  fused_point_b_.clear();
  fused_point_visibility_.clear();
//  fused_point_col_.clear();
//  fused_point_row_.clear();

  while (!fusion_queue_.empty()) {
    // this information gets updated inside the while loop to always reflect the current image
    const auto data = fusion_queue_.back();
    const int image_idx = data.image_idx;
    const int row = data.row;
    const int col = data.col;
    // in the first run, coming from the row/col iterator this is 0
    const int traversal_depth = data.traversal_depth;

    fusion_queue_.pop_back();

    // Check if pixel (in this current image, not the global 3d coordinate) already processed.
    auto& fused_pixel_mask = fused_pixel_masks_.at(image_idx);
    if (fused_pixel_mask.Get(row, col)) {
      continue;
    }

    const auto& depth_map = workspace_->GetDepthMap(image_idx);
    const float depth = depth_map.Get(row, col);

    // Pixels with negative depth are filtered.
    if (depth <= 0.0f) {
      continue;
    }

    // If the traversal depth is greater than zero, the initial reference
    // pixel has already been added and we need to check for consistency.
    if (traversal_depth > 0) {
      // Project reference point into current view.
      const Eigen::Vector3f proj = P_.at(image_idx) * fused_ref_point;

      // Depth error of reference depth with current depth.
      const float depth_error = std::abs((proj(2) - depth) / depth);
      if (depth_error > options_.max_depth_error) {
        continue;
      }

      // Reprojection error reference point in the current view.
      const float col_diff = proj(0) / proj(2) - col;
      const float row_diff = proj(1) / proj(2) - row;
      const float squared_reproj_error =
          col_diff * col_diff + row_diff * row_diff;
      if (squared_reproj_error > max_squared_reproj_error_) {
        continue;
      }
    }

    // Determine normal direction in global reference frame.
    const auto& normal_map = workspace_->GetNormalMap(image_idx);
    const Eigen::Vector3f normal =
        inv_R_.at(image_idx) * Eigen::Vector3f(normal_map.Get(row, col, 0),
                                               normal_map.Get(row, col, 1),
                                               normal_map.Get(row, col, 2));

    // Check for consistent normal direction with reference normal.
    if (traversal_depth > 0) {
      const float cos_normal_error = fused_ref_normal.dot(normal);
      if (cos_normal_error < min_cos_normal_error_) {
        continue;
      }
    }

    // Determine 3D location of current depth value.
    const Eigen::Vector3f xyz =
        inv_P_.at(image_idx) *
        Eigen::Vector4f(col * depth, row * depth, depth, 1.0f);

    // Read the color of the pixel.
    BitmapColor<uint8_t> color;
    const auto& bitmap_scale = bitmap_scales_.at(image_idx);
    workspace_->GetBitmap(image_idx).InterpolateNearestNeighbor(
        col / bitmap_scale.first, row / bitmap_scale.second, &color);

    // Set the current pixel as visited.
    fused_pixel_mask.Set(row, col, true);

    // Accumulate statistics for fused point.
    fused_point_x_.push_back(xyz(0));
    fused_point_y_.push_back(xyz(1));
    fused_point_z_.push_back(xyz(2));
    fused_point_nx_.push_back(normal(0));
    fused_point_ny_.push_back(normal(1));
    fused_point_nz_.push_back(normal(2));
    fused_point_r_.push_back(color.r);
    fused_point_g_.push_back(color.g);
    fused_point_b_.push_back(color.b);
    fused_point_visibility_.insert(image_idx);
//    fused_point_visibility_.push_back(image_idx);
//    fused_point_col_.push_back(col);
//    fused_point_row_.push_back(row);

    // Remember the first pixel as the reference.
    if (traversal_depth == 0) {
      fused_ref_point = Eigen::Vector4f(xyz(0), xyz(1), xyz(2), 1.0f);
      fused_ref_normal = normal;
    }

    if (fused_point_x_.size() >= static_cast<size_t>(options_.max_num_pixels)) {
      break;
    }

    // now go to the next image
    FusionData next_data;
    // so traversal depth means how many images are traversed per pixel
    next_data.traversal_depth = traversal_depth + 1;

    if (next_data.traversal_depth >= options_.max_traversal_depth) {
      continue;
    }

    // now for the calculated 3d coordinate,
    // project it back to all images that are overlapping with the current image, and push it back
    // into the fusion queue
    // this will not be done again for every reprojected pixel, because the while loop will be continue
    // before it reaches here, if the pixel was already processed
    for (const auto next_image_idx : overlapping_images_.at(image_idx)) {
      if (!used_images_.at(next_image_idx) ||
          fused_images_.at(next_image_idx)) {
        continue;
      }

      next_data.image_idx = next_image_idx;

      const Eigen::Vector3f next_proj =
          P_.at(next_image_idx) * xyz.homogeneous();
      next_data.col = static_cast<int>(std::round(next_proj(0) / next_proj(2)));
      next_data.row = static_cast<int>(std::round(next_proj(1) / next_proj(2)));

      const auto& depth_map_size = depth_map_sizes_.at(next_image_idx);
      if (next_data.col < 0 || next_data.row < 0 ||
          next_data.col >= depth_map_size.first ||
          next_data.row >= depth_map_size.second) {
        continue;
      }

      fusion_queue_.push_back(next_data);
    }
  } // end of while loop over one pixel/point in several images

  fusion_queue_.clear();

  const size_t num_pixels = fused_point_x_.size();
  if (num_pixels >= static_cast<size_t>(options_.min_num_pixels)){

    PlyPoint fused_point;

    Eigen::Vector3f fused_normal;
    fused_normal.x() = myInternal::Median(&fused_point_nx_);
    fused_normal.y() = myInternal::Median(&fused_point_ny_);
    fused_normal.z() = myInternal::Median(&fused_point_nz_);
    const float fused_normal_norm = fused_normal.norm();
    if (fused_normal_norm < std::numeric_limits<float>::epsilon()) {
      return;
    }

    //    Point p(myInternal::Median(&fused_point_x_),
    //            myInternal::Median(&fused_point_y_),
    //            myInternal::Median(&fused_point_z_));

    //    Vector n(fused_normal.x() / fused_normal_norm,
    //            fused_normal.y() / fused_normal_norm,
    //            fused_normal.z() / fused_normal_norm);

    //    vertex_info v;
    //    unsigned char red = TruncateCast<float, uint8_t>(
    //        std::round(myInternal::Median(&fused_point_r_)));
    //    unsigned char green = TruncateCast<float, uint8_t>(
    //        std::round(myInternal::Median(&fused_point_g_)));
    //    unsigned char blue = TruncateCast<float, uint8_t>(
    //        std::round(myInternal::Median(&fused_point_b_)));
    //    Color col = CGAL::make_array(red, green, blue);
    //    v.color = col;
    //    v.normal = n;


    //    points.push_back(p);
    //    fused_points_visibility_.emplace_back(fused_point_visibility_.begin(),
    //                                          fused_point_visibility_.end());



    fused_point.x = myInternal::Median(&fused_point_x_);
    fused_point.y = myInternal::Median(&fused_point_y_);
    fused_point.z = myInternal::Median(&fused_point_z_);

    fused_point.nx = fused_normal.x() / fused_normal_norm;
    fused_point.ny = fused_normal.y() / fused_normal_norm;
    fused_point.nz = fused_normal.z() / fused_normal_norm;

    fused_point.r = TruncateCast<float, uint8_t>(
        std::round(myInternal::Median(&fused_point_r_)));
    fused_point.g = TruncateCast<float, uint8_t>(
        std::round(myInternal::Median(&fused_point_g_)));
    fused_point.b = TruncateCast<float, uint8_t>(
        std::round(myInternal::Median(&fused_point_b_)));


    fused_points_.push_back(fused_point);
    fused_points_visibility_.emplace_back(fused_point_visibility_.begin(),
                                          fused_point_visibility_.end());
  }
}

void DepthMapConverter::getFusedPoints(std::vector<Point>& points, std::vector<vertex_info>& infos, std::map<int, Point>& sensorMap){

    std::vector<int>::size_type i=0;
    for(PlyPoint fused_point : fused_points_){

      Point p(fused_point.x,
              fused_point.y,
              fused_point.z);

      Vector n(fused_point.nx,
              fused_point.ny,
              fused_point.nz);

      vertex_info v;
      unsigned char red = fused_point.r;
      unsigned char green = fused_point.g;
      unsigned char blue = fused_point.b;
      Color col = CGAL::make_array(red, green, blue);
      v.color = col;
      v.normal = n;
//      v.images = fused_points_visibility_[i];

      Point sensor_point = sensorMap.find(fused_points_visibility_[i][0])->second;
      //v.sensor_pos = sensor_point;
      v.sensor_vec = sensor_point - p;

      i++;

      points.push_back(p);
      infos.push_back(v);
  }


}

void DepthMapConverter::getSensorMap(std::map<int, Point>& sensorMap){

    Workspace::Options workspace_options;

    auto workspace_format_lower_case = workspace_format_;
    StringToLower(&workspace_format_lower_case);
    if (workspace_format_lower_case == "pmvs") {
      workspace_options.stereo_folder =
          StringPrintf("stereo-%s", pmvs_option_name_.c_str());
    }

    workspace_options.max_image_size = options_.max_image_size;
    workspace_options.image_as_rgb = true;
    workspace_options.cache_size = options_.cache_size;
    workspace_options.workspace_path = workspace_path_;
    workspace_options.workspace_format = workspace_format_;
    workspace_options.input_type = input_type_;

    workspace_.reset(new Workspace(workspace_options));

    const auto& model = workspace_->GetModel();

    const auto image_names = ReadTextFileLines(JoinPaths(
        workspace_path_, workspace_options.stereo_folder, "fusion.cfg"));


    for (const auto& image_name : image_names) {
        const int image_idx = model.GetImageIdx(image_name);
        const auto& image = model.images.at(image_idx);

        // get the sensor position
        auto R = Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(image.GetR());
        auto T = Eigen::Map<const Eigen::Vector3f>(image.GetT());
        Point sensor_position = eigen2Cgal(-R.transpose()*T);

        std::cout << image_idx << " " << image_name << " " << sensor_position << std::endl;
//        std::cout << sensor_position << std::endl;

        sensorMap[image_idx] = sensor_position;

    }
}



void DepthMapConverter::getFusedPointsVisibility(std::vector<std::vector<int>>& visibility)
    const {
  visibility = fused_points_visibility_;
}



}  // namespace mvs
}  // namespace colmap
