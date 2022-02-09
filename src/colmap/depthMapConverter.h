//// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
//// All rights reserved.
////
//// Redistribution and use in source and binary forms, with or without
//// modification, are permitted provided that the following conditions are met:
////
////     * Redistributions of source code must retain the above copyright
////       notice, this list of conditions and the following disclaimer.
////
////     * Redistributions in binary form must reproduce the above copyright
////       notice, this list of conditions and the following disclaimer in the
////       documentation and/or other materials provided with the distribution.
////
////     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
////       its contributors may be used to endorse or promote products derived
////       from this software without specific prior written permission.
////
//// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//// POSSIBILITY OF SUCH DAMAGE.
////
//// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

//#ifndef COLMAP_SRC_MVS_DEPTHMAPCONVERTER_H_
//#define COLMAP_SRC_MVS_DEPTHMAPCONVERTER_H_

//#include <unordered_set>
//#include <vector>

//#include <Eigen/Core>

//#include "mvs/fusion.h"
//#include "mvs/depth_map.h"
//#include "mvs/image.h"
//#include "mvs/mat.h"
//#include "mvs/model.h"
//#include "mvs/normal_map.h"
//#include "mvs/workspace.h"
//#include "util/alignment.h"
//#include "util/cache.h"
//#include "util/math.h"
//#include "util/ply.h"
//#include "util/threading.h"

//#include <base/cgal_typedefs.h>

//namespace colmap{
//namespace mvs{

////class colmapImage{

//// public:

////   colmapImage();
////   colmapImage(const int width,
////               const int height);

////   colmapImage(const int width,
////               const int height,
////               std::vector<int> data);

//// private:

////   const int width;
////   const int height;

////};

////struct image_info{
////    int idx;
////    int nrow;
////    int ncol;
////    std::unordered_set<int> image_point_idx;
////};
//// TODO:
//// the question is, do I form tetrahedrons, even of points that are not used
//// there are basically 3 options,
//// 1. I only form sensor tetrahedrons per image of points that are used in the final point cloud / after the fusion step
//// 2. I form sensor tetrahedrons of all points (that have a depth) in all the depth maps
//        // papers that do ray tracing usually choose this option.
//// 3. I form sensor tetrahedrons of fused points, even between images
//        // 3. is most complex to implement, as I need to know what is a neighbouring point, even when it was not selected, and what is
//        // this neighbouring points depth neighbour. and what to do if it has multiple depth neighbours??


//class DepthMapConverter : protected StereoFusion{

//public:

//  DepthMapConverter(const StereoFusionOptions& options,
//                const std::string& workspace_path,
//                const std::string& workspace_format,
//                const std::string& pmvs_option_name,
//                const std::string& input_type);

////  std::vector<int> fused_point_row_;
////  std::vector<int> fused_point_col_;
////  std::vector<int> fused_point_visibility_;

//  void run();
//  void fuse();

//  void getFusedPoints(std::vector<Point>& points, std::vector<vertex_info>& infos, std::map<int, Point>& sensorMap);
//  void getFusedPointsVisibility(std::vector<std::vector<int>>&) const;
//  void getSensorMap(std::map<int, Point>& sensorMap);


//}; // end of DepthMapConverter class definition

//}  // namespace mvs
//}  // namespace colmap

//#endif  // COLMAP_SRC_MVS_FUSION_H_
