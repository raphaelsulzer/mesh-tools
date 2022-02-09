#include <base/cgal_typedefs.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <IO/fileIO.h>

#include <util/helper.h>
#include <util/vectorArithmetic.h>
#include <processing/pointSetProcessing.h>

#include <rPLY/ply.h>

// colmap util
#include <util/misc.h>
#include <util/endian.h>
//#include <util/ply.h>
#include <mvs/workspace.h>
//#include <mvs/fusion.h>

#ifdef Open3D
#include <Open3D/Open3D.h>
#include <Open3D/Utility/IJsonConvertible.h>
#endif

#include <CGAL/Bbox_3.h>


/////// read colmap reconstructions ///////////
void getSensorMap(const std::string path, std::map<int, Point>& sensorMap){

//    colmap::mvs::StereoFusionOptions options;

    colmap::mvs::Workspace::Options workspace_options;
    const std::string workspace_format = "COLMAP";
    const std::string input_type = "geometric";

//    workspace_options.max_image_size = options.max_image_size;
    workspace_options.max_image_size = -1;
    workspace_options.image_as_rgb = true;
//    workspace_options.cache_size = options.cache_size;
    workspace_options.cache_size = 32.0;
    workspace_options.workspace_path = path;
    workspace_options.workspace_format = workspace_format;
    workspace_options.input_type = input_type;

    auto ws = new colmap::mvs::Workspace(workspace_options);

    const auto& model = ws->GetModel();

    const auto image_names = colmap::ReadTextFileLines(colmap::JoinPaths(
        path, workspace_options.stereo_folder, "fusion.cfg"));

    for (const auto& image_name : image_names) {
        const int image_idx = model.GetImageIdx(image_name);
        const auto& image = model.images.at(image_idx);

        // get the sensor position
        auto R = Eigen::Map<const Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(image.GetR());
        auto T = Eigen::Map<const Eigen::Vector3f>(image.GetT());
        Point sensor_position = eigen2Cgal(-R.transpose()*T);
//        std::cout << image_idx << " " << image_name << " " << sensor_position << std::endl;
//        std::cout << sensor_position << std::endl;
        sensorMap[image_idx] = sensor_position;
    }
}

bool sortSensors(const std::pair<double, Point> &a, const std::pair<double, Point> &b)
{return (a.first < b.first);}

void readColmapRecon(dirHolder& dir, dataHolder& data){
    auto start = std::chrono::high_resolution_clock::now();
    cout << "\nRead colmap reconstruction: " << endl;
    getSensorMap(dir.path, data.sensor_map);

    // apply transformation to register tanks and temples colmap reconstruction with tanks and temples ground truth
    bool applyTrafo = true;
    cout << "\t-read transformation matrix " << dir.read_file+"_trans.txt" << endl;
    Eigen::Matrix4d transformation_matrix = importTranslationMatrix(dir.path+dir.read_file+"_trans.txt");
//    cout << "\t\t" << transformation_matrix << endl;

    open3d::geometry::PointCloud gt_pc;
    cout << "\t-read ground truth from " << dir.read_file + ".ply" << endl;
    open3d::io::ReadPointCloudFromPLY(dir.path+dir.read_file+".ply",gt_pc);
//    const auto bbox = gt_pc.GetOrientedBoundingBox();
    const auto bbox = gt_pc.GetAxisAlignedBoundingBox().Scale(1.5);

    // read in the points from the fused and fused.vis files
    const auto& ply_points = colmap::ReadPly(colmap::JoinPaths(dir.path, "fused.ply"));
    const std::string vis_path = colmap::JoinPaths(dir.path, "fused.ply.vis");
    fstream vis_file(vis_path, std::ios::in | std::ios::binary);
    CHECK(vis_file.is_open()) << vis_path;
    const size_t vis_num_points = colmap::ReadBinaryLittleEndian<uint64_t>(&vis_file);
    CHECK_EQ(vis_num_points, ply_points.size());

    data.points.reserve(ply_points.size());
    vector<Eigen::Vector3d> points_eigen;
    vector<Eigen::Vector3d> sensors_eigen;
//    vector<Eigen::Vector3d> colors_eigen;
    for (const auto& ply_point : ply_points){

        points_eigen.push_back(Eigen::Vector3d(ply_point.x, ply_point.y, ply_point.z));
        Point input_point(ply_point.x, ply_point.y, ply_point.z);

        vertex_info vinf;
        vinf.normal = Vector(ply_point.nx, ply_point.ny, ply_point.nz);
        // color
        unsigned char r,g,b;
        r = ply_point.r;
        g = ply_point.g;
        b = ply_point.b;
        std::array<unsigned char, 3> col = {r,g,b};
        vinf.color = col;

        // make a list of visible images, sorted by the angle between the sensor_vec and the normal of the point, smallest angles first
        int num_visible_images = colmap::ReadBinaryLittleEndian<uint32_t>(&vis_file);
        std::vector<std::pair<double,Point>> sensor_points;
        for (uint32_t i = 0; i < num_visible_images; ++i) {
            const int image_idx = colmap::ReadBinaryLittleEndian<uint32_t>(&vis_file);
            Point sensor_point = data.sensor_map.find(image_idx)->second;
            Vector v1 = sensor_point - input_point;
            Vector v2 = vinf.normal;
            double angle = std::acos( v1 * v2 / std::sqrt(v1.squared_length()*v2.squared_length()));
            if(CGAL::cross_product(v1,v2).z()<0)
                angle = 2*M_PI-angle;
            sensor_points.push_back(std::make_pair(angle, sensor_point));
        }
        sort(sensor_points.begin(), sensor_points.end(), sortSensors);

        for(auto& sensor_position : sensor_points){
            if(applyTrafo)
                applyTransformation(sensor_position.second, transformation_matrix);
            vinf.sensor_positions.push_back(sensor_position.second);
        }
        sensors_eigen.push_back(Eigen::Vector3d(vinf.sensor_positions[0].x(), vinf.sensor_positions[0].y(), vinf.sensor_positions[0].z()));
//        colors_eigen.push_back(Eigen::Vector3d(ply_point.x, ply_point.y, ply_point.z));
        if(applyTrafo)
            applyTransformation(input_point, transformation_matrix);
        if(pointInsideBbox(input_point,bbox)){
            data.points.push_back(input_point);
            vinf.sensor_vec = sensor_points[0].second-input_point;
            data.infos.push_back(vinf);
        }

    }


//    // load ground truth, make a bounding box and clip the reconstruction with this bounding box
//    const string cropfile = dir.path + dir.read_file + "_cropfile.json";
//    cout << "\t-read cropfile " << dir.read_file + "_cropfile.json" << endl;
//    open3d::visualization::SelectionPolygonVolume bounding_poly;
//    open3d::io::ReadIJsonConvertible(cropfile, bounding_poly);
//    open3d::geometry::PointCloud pc;
//    pc.normals_ = sensors_eigen;
//    pc.points_ = points_eigen;

//    auto pc_transformed = pc.Transform(transformation_matrix);

//    cout << "size before cropping: " << pc_transformed.points_.size() << endl;
//    auto pc_cropped = bounding_poly.CropPointCloud(pc_transformed);
//    cout << "size after cropping: " << pc_cropped->points_.size() << endl;
//    data.points.clear(); data.infos.clear();
//    for(int i = 0; i < pc_cropped->points_.size(); i++){
//        data.points.push_back(Point(pc_cropped->points_[i].x(), pc_cropped->points_[i].y(), pc_cropped->points_[i].z()));
//        vertex_info vinf;
//        vinf.sensor_positions.push_back(Point(pc_cropped->normals_[i].x(), pc_cropped->normals_[i].y(), pc_cropped->normals_[i].z()));
//        data.infos.push_back(vinf);
//    }


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\t-" << data.points.size() << " points read" << endl;
    cout << "\t-" << data.sensor_map.size() << " cameras read" << endl;
//    cout << "\t-camera positions of sensor_vec which is closest to normal will be exported to scalar_x0, scalar_y0, scalar_z0" << endl;
    cout << "\t-in " << duration.count() << "s" << endl;

}

//void readColmapRecon(dirHolder& dir, dataHolder& data){
//    auto start = std::chrono::high_resolution_clock::now();
//    cout << "\nRead colmap reconstruction: " << endl;
//    getSensorMap(dir.path, data.sensor_map);

//    // apply transformation to register tanks and temples colmap reconstruction with tanks and temples ground truth
//    bool applyTrafo = false;
//    Eigen::Matrix4d transformation_matrix;
//    if(applyTrafo){
//        // load trans file, e.g. Barn_trans

//        transformation_matrix << -3.024223114080324848e+00, 1.337200447181571550e-01, 2.822018693855588012e+00, -8.033102841327895760e+00,
//                                 -2.824574281070231230e+00, -2.290640822308912528e-01, -3.016107720144580728e+00, 6.856971949289079049e-01,
//                                 5.874257013764642293e-02, -4.130041718728730160e+00, 2.586517244893801193e-01, -2.187313805091164554e+01,
//                                 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00;
//    }
//    // load ground truth, make a bounding box and clip the reconstruction with this bounding box
//    const string cropfile = dir.path + dir.read_file + "_cropfile.json";
//    open3d::visualization::SelectionPolygonVolume bounding_poly;
//    open3d::io::ReadIJsonConvertible( cropfile, bounding_poly);
//    open3d::geometry::PointCloud pc;
//    vector<Eigen::Vector3d> normals;
//    pc.normals_ = normals;


//    // read in the points from the fused and fused.vis files
//    const auto& ply_points = colmap::ReadPly(colmap::JoinPaths(dir.path, "fused.ply"));
//    const std::string vis_path = colmap::JoinPaths(dir.path, "fused.ply.vis");
//    fstream vis_file(vis_path, std::ios::in | std::ios::binary);
//    CHECK(vis_file.is_open()) << vis_path;
//    const size_t vis_num_points = colmap::ReadBinaryLittleEndian<uint64_t>(&vis_file);
//    CHECK_EQ(vis_num_points, ply_points.size());

//    data.points.reserve(ply_points.size());
//    vector<Eigen::Vector3d> points;
//    vector<Eigen::Vector3d> sensors;
//    for (const auto& ply_point : ply_points){


////        Eigen::Vector3d input_point(ply_point.x, ply_point.y, ply_point.z);
////        points.push_back(input_point);

//        Point input_point(ply_point.x, ply_point.y, ply_point.z);

//        vertex_info vinf;
//        vinf.normal = Vector(ply_point.nx, ply_point.ny, ply_point.nz);
//        // color
//        unsigned char r,g,b;
//        r = ply_point.r;
//        g = ply_point.g;
//        b = ply_point.b;
//        std::array<unsigned char, 3> col = {r,g,b};
//        vinf.color = col;

//        // make a list of visible images, sorted by the angle between the sensor_vec and the normal of the point, smallest angles first
//        int num_visible_images = colmap::ReadBinaryLittleEndian<uint32_t>(&vis_file);
//        std::vector<std::pair<double,Point>> sensor_points;
//        for (uint32_t i = 0; i < num_visible_images; ++i) {
//            const int image_idx = colmap::ReadBinaryLittleEndian<uint32_t>(&vis_file);
//            Point sensor_point = data.sensor_map.find(image_idx)->second;
//            Vector v1 = sensor_point - input_point;
//            Vector v2 = vinf.normal;
//            double angle = std::acos( v1 * v2 / std::sqrt(v1.squared_length()*v2.squared_length()));
//            if(CGAL::cross_product(v1,v2).z()<0)
//                angle = 2*M_PI-angle;
//            sensor_points.push_back(std::make_pair(angle, sensor_point));
//        }
//        std::sort(sensor_points.begin(), sensor_points.end(), sortSensors);

//        for(auto& sensor_position : sensor_points){
//            if(applyTrafo)
//                applyTransformation(sensor_position.second, transformation_matrix);
//            vinf.sensor_positions.push_back(sensor_position.second);
//        }
//        data.points.push_back(input_point);
//        vinf.sensor_vec = sensor_points[0].second-input_point;
//        data.infos.push_back(vinf);

//    }

//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    cout << "\t-" << data.points.size() << " points read" << endl;
//    cout << "\t-" << data.sensor_map.size() << " cameras read" << endl;
//    cout << "\t-camera positions of sensor_vec which is closest to normal will be exported to scalar_x0, scalar_y0, scalar_z0" << endl;
//    cout << "\t-in " << duration.count() << "s" << endl;

//}

//void readColmapPoints(const dirHolder& dir, dataHolder& data)
//{
//    // this reads a file that has been processed with the readColmapRecon() function before.
//    // Meaning a points file with sensor data, saved as nx, ny, nz

//    auto start = std::chrono::high_resolution_clock::now();
//    cout << "\nRead colmap points from " << dir.read_file << endl;

//    std::string ifn = dir.path + dir.read_file + ".ply";

//    // read Binary PLY with sensor
//    Mesh_ply aMesh;
//    Import_PLY(ifn.c_str(), &aMesh);

//    for(int i = 0; i < aMesh.mVertices.size(); i++){
//        // save points
//        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
//        data.points.push_back(pt);
//        // sensor
//        vertex_info vec_info;
//        vec_info.sensor_positions.push_back(Point(aMesh.mNormals[i].x, aMesh.mNormals[i].y, aMesh.mNormals[i].z));
//        data.infos.push_back(vec_info);
//    }
//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    cout << "\t-" << data.points.size() << " points read" << endl;
//    cout << "\t-in " << duration.count() << "s" << endl;
//}





////////////////////////////////////////////////
//////////////// WU TENGS CODE /////////////////
////////////////////////////////////////////////
// for colmap fused.ply.vis(binary format)
bool loadVisibilityFile(const char* pszVisFile, std::vector<std::vector<int>>& vis_info)
{
    FILE * fp = fopen(pszVisFile, "rb");
    if (fp == NULL){
        printf("Can not open the file %s\n", pszVisFile);
        return false;
    }

    uint64_t vis_num_points = 0;
    fread(& vis_num_points, sizeof(uint64_t), 1, fp);

    if (vis_num_points < 0){
        printf("non visibilty information\n");
        fclose(fp);

        return true;
    }

    vis_info.reserve(vis_num_points);
    for (int i = 0; i < vis_num_points; ++i){
        uint32_t num_visible_images = 0;
        fread(& num_visible_images, sizeof(uint32_t), 1, fp);

        std::vector<int> image_view;
        image_view.reserve(num_visible_images);
        for (int j = 0; j < num_visible_images; ++j){
            uint32_t image_idx = 0;
            fread(& image_idx, sizeof(uint32_t), 1, fp);
            image_view.push_back(image_idx);
        }
        vis_info.push_back(image_view);
    }

    fclose(fp);
    return true;
}

// because colmap save qvec and tvec
// we need cvec
Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d& qvec){
  const double norm = qvec.norm();
  if (norm == 0) {
    // We do not just use (1, 0, 0, 0) because that is a constant and when used
    // for automatic differentiation that would lead to a zero derivative.
    return Eigen::Vector4d(1.0, qvec(1), qvec(2), qvec(3));
  } else {
    return qvec / norm;
  }
}

//Eigen::Vector3d ProjectionCenterFromPose(const Eigen::Vector4d& qvec,
//                                         const Eigen::Vector3d& tvec){
//  // Inverse rotation as conjugate quaternion.
//  const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
//  const Eigen::Quaterniond quat(normalized_qvec(0), -normalized_qvec(1),
//                                -normalized_qvec(2), -normalized_qvec(3));
//  return quat * -tvec;
//}


Eigen::Vector3d ProjectionCenterFromPose(const double* qvec,
                                         const Eigen::Vector3d& tvec){
  // Inverse rotation as conjugate quaternion.
//  const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R;
  ceres::QuaternionToRotation(qvec, R.data());

//  const Eigen::Quaterniond quat(normalized_qvec(0), -normalized_qvec(1),
//                                -normalized_qvec(2), -normalized_qvec(3));
  return -R.transpose() * tvec;
}


bool loadImageFile(const char* pszImageFile, std::map<int, Point>& image_proj_center)
{
    FILE * fp = fopen(pszImageFile, "rb");
    if (fp == NULL){
        printf("Can not open the file %s\n", pszImageFile);
        return false;
    }

    uint64_t num_reg_images = 0;

    fread(& num_reg_images, sizeof(uint64_t), 1, fp);

    if (num_reg_images < 0){
        printf("non image information\n");
        fclose(fp);

        return true;
    }

    for (int i = 0; i < num_reg_images; ++i){
        uint32_t num_image_index = 0;
        fread(& num_image_index, sizeof(uint32_t), 1, fp);

        double qvec_array[4] = { 0 };
        fread(qvec_array, sizeof(double), 4, fp);

        double tvec_array[3] = { 0 };
        fread(tvec_array, sizeof(double), 3, fp);

        Eigen::Vector4d qvec;
        for (int j = 0; j < 4; ++j){
            qvec(j) = qvec_array[j];
        }

        Eigen::Vector3d tvec;
        for (int j = 0; j < 3; ++j){
            tvec(j) = tvec_array[j];
        }

//        Eigen::Vector3d cvec = ProjectionCenterFromPose(qvec, tvec);
        Eigen::Vector3d cvec = ProjectionCenterFromPose(qvec_array, tvec);

        Point proj_center(cvec(0), cvec(1), cvec(2));

        uint32_t num_camera_index = 0;
        fread(& num_camera_index, sizeof(uint32_t), 1, fp);

        // read char in binary
        std::string img_name;

        char name_char;
        do {
            fread(& name_char, sizeof(char), 1, fp);
            if (name_char != '\0') {
                img_name += name_char;
            }
        } while (name_char != '\0');

        uint64_t num_points2D = 0;
        fread(& num_points2D, sizeof(uint64_t), 1, fp);

        for (int j = 0; j < num_points2D; ++j){
            double xy[2] = { 0 };
            fread(xy, sizeof(double), 2, fp);

            uint64_t num_point3d_index = 0;
            fread(& num_point3d_index, sizeof(uint64_t), 1, fp);
        }

        std::cout << num_camera_index << " " << img_name << " " << proj_center << std::endl;

        image_proj_center.insert(std::map<int, Point>::value_type(num_camera_index, proj_center));
    }
    fclose(fp);

    return true;
}




