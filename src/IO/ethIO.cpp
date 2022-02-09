#include <base/cgal_typedefs.h>

#include <IO/fileIO.h>
#include <IO/ttIO.h>
#include <IO/ethIO.h>
#include <util/vectorArithmetic.h>
#include <util/helper.h>
#include <learning/learning.h>
#include <IO/meshlab_project.h>
#include <rPLY/ply.h>
#include <processing/pointSetProcessing.h>

#include <Eigen/Core>
#include <boost/filesystem.hpp>

#include <pcl/console/parse.h>
#include <pcl/io/ply_io.h>
#include <pcl/common/transforms.h>
#include <pcl/search/kdtree.h>

#include "IO/meshlab_project.h"
#include "IO/util.h"

using namespace std;

////////////////////////////////////////////////////////////
/////////////////////// ETH3D I/O ///////////////////////////
////////////////////////////////////////////////////////////

int importETH3D(dirHolder dir, dataHolder& data){

    cout << "\nLoad ETH3D meshlab project..." << endl;


    pcl::console::setVerbosityLevel(pcl::console::L_ALWAYS);

    // Load the ground truth point cloud poses from the MeshLab project file.
    MeshLabMeshInfoVector scan_infos;

    string ground_truth_mlp_path = dir.path + "scan_clean/scan_alignment.mlp";

    if (!ReadMeshLabProject(ground_truth_mlp_path, &scan_infos)){
      std::cerr << "Cannot read scan poses from " << ground_truth_mlp_path
                << std::endl;
      return 1;
    }

    // Validate the scan transformations (partly: not checking that the top-left
    // 3x3 block is a rotation matrix).
    for (size_t scan_index = 0; scan_index < scan_infos.size(); ++scan_index) {
      const MeshLabProjectMeshInfo& info = scan_infos[scan_index];
      if (info.global_T_mesh(3, 0) != 0 || info.global_T_mesh(3, 1) != 0 ||
          info.global_T_mesh(3, 2) != 0 || info.global_T_mesh(3, 3) != 1) {
        std::cerr << "Error: Last row in a scan's transformation matrix is not"
                     " (0, 0, 0, 1)."
                  << std::endl;
        return 1;
      }
    }

    // Load the ground truth scan point clouds.
    std::vector<PointCloudPtr> scans;
    for (const MeshLabProjectMeshInfo& scan_info : scan_infos) {
      // Get absolute or compose relative path.
      std::string file_path =
          (scan_info.filename.empty() || scan_info.filename[0] == '/')
              ? scan_info.filename
              : (boost::filesystem::path(ground_truth_mlp_path).parent_path() /
                 scan_info.filename)
                    .string();

      std::cout << "\t-loading scan: " << scan_info.filename << std::endl;
      PointCloudPtr point_cloud(new PointCloud());
      if (pcl::io::loadPLYFile(file_path, *point_cloud) < 0) {
        std::cerr << "Cannot read scan file." << std::endl;
        return 1;
      }

      scans.push_back(point_cloud);
    }

//    PointCloudPtr scan(new PointCloud());
    for (size_t scan_index = 0; scan_index < scan_infos.size(); ++scan_index){
        PointCloud temp_cloud;
        pcl::transformPointCloud(*scans[scan_index], temp_cloud,
                               scan_infos[scan_index].global_T_mesh);
        for (const auto & pt : temp_cloud.points){
            Point p(pt.x, pt.y, pt.z );
            if(!pointInsideBbox(p,data.sm_obb))
                continue;
            data.points.push_back(p);
            Point s(scan_infos[scan_index].global_T_mesh(0,3),
                  scan_infos[scan_index].global_T_mesh(1,3),
                  scan_infos[scan_index].global_T_mesh(2,3));
            vertex_info vinf;
            vinf.sensor_positions.push_back(s);
            vinf.sensor_vec = s - p;
            vinf.sensor_vec = vinf.sensor_vec/std::sqrt(vinf.sensor_vec.squared_length());
            data.infos.push_back(vinf);
        }
    }

    data.has_sensor = true;
    data.has_normal = false;

    return 0;
}

//int pcl2cgal(dataHolder& data, PointCloudPtr& cloud)
//{

//    for (long long int scan_point_index = 0; scan_point_index < scan_point_size;
//         ++scan_point_index) {
//      const pcl::PointXYZ& scan_point = scan->at(scan_point_index);

//    vertex_info vinf;
//    for (const auto & pt : *cloud.points) {
//    Point p (pt.x, pt.y, pt.z );
//    data.points.push_back(p);
//    Vector n (pt.normal_x, pt.normal_y, pt.normal_z);
//    vinf.normal = n;
//    vinf.sensor_positions[0] = p;
//    data.infos.push_back(vinf);
//    }

//    return 0;
//}




