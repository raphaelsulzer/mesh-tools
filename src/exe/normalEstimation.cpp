#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#include <exe/inputParser.h>

#include <processing/normalAndSensorProcessing.h>

#include <boost/filesystem.hpp>
using namespace std;

//int normalEstimationTanksAndTemples(dirHolder& dir, runningOptions options){

//    string fp = "scanner_pos.txt";
//    dir.path+="is/";
//    ifstream file(dir.path+fp);

//    cout << "\nLoad sensor position files..." << endl;
//    cout << "\t-from " << fp << endl;

//    if(!file){
//        cout << "\nERROR: FILE DOES NOT EXIST OR IS EMPTY!" << endl;
//        return 1;
//    }

//    string line;
//    vector<string> lineAsVector;
//    vector<dataHolder> datas;
//    for(int i = 0; i < options.number_of_scans; i++){
//        if(!getline(file, line))
//            break;
//        dataHolder data;
//        splitString(line, lineAsVector, ' ');
//        dir.read_file = lineAsVector[0].c_str();
//        importTanksAndTemplesIs(dir, data, Point(stod(lineAsVector[1]),stod(lineAsVector[2]),stod(lineAsVector[3])));
////        estimateNormals(data,options.neighborhood_size,options.orientation);
//        datas.push_back(data);
//    };

//    vector<dataHolder>::iterator it = datas.begin();
//    it++;
//    while(it != datas.end()){
//        datas[0].points.insert(datas[0].points.end(),it->points.begin(),it->points.end());
//        datas[0].infos.insert(datas[0].infos.end(),it->infos.begin(),it->infos.end());
//        it=datas.erase(it);
//    }



//    if(options.subsample_grid_spacing > 0)
//        gridSimplify(datas[0], options.subsample_grid_spacing);

//    estimateNormals(datas[0],options.neighborhood_size,options.orientation);


//    exportOptions eo;
//    eo.sensor_vec = false;
//    eo.color = false;
//    eo.normals = true;
//    dir.suffix = "_normals";
//    if(dir.write_file.empty())
//        dir.write_file = dir.read_file;

//    exportPoints(dir,datas[0],eo);
//}


//int normalEstimationETH3D(dirHolder& dir, runningOptions options){

//    dataHolder mvs_data;
//#ifdef OpenMVS
//    if(loadOMVSScene(dir, mvs_data))
//        return 1;
//#endif
//    string input = dir.read_file;

//    dataHolder lidar_data;
//    if(input!="mvs"){
//        dir.read_file = "dslr_scan_eval/scan12";
//        cout << "\nLoad scanner files..." << endl;
//        cout << "\t-from " << dir.path+dir.read_file << endl;

//        importLidarPoints(dir,lidar_data);
//        int orient=0;
//        estimateNormals(lidar_data,options.neighborhood_size,orient);
//    }

//    if(input=="lidar"){
//        orientNormalsWithSecondCloud(mvs_data,lidar_data);
//        mvs_data.points.swap(lidar_data.points), lidar_data.points.clear();
//        mvs_data.infos.swap(lidar_data.infos), lidar_data.infos.clear();
//        dir.write_file=dir.read_file;
//    }
//    else if(input=="mvs"){
//        dir.write_file= "dslr_scan_eval/densify_file";
//    }
//    else if(input =="all"){
//        orientNormalsWithSecondCloud(mvs_data,lidar_data);
//        mvs_data.points.insert(mvs_data.points.end(), lidar_data.points.begin(), lidar_data.points.end());
//        lidar_data.points.clear();
//        mvs_data.infos.insert(mvs_data.infos.end(), lidar_data.infos.begin(), lidar_data.infos.end());
//        lidar_data.infos.clear();
//        dir.write_file= "dslr_scan_eval/all";
//    }
//    else{
//        cout << "ERROR: not a valid input. Choose either lidar, mvs or all." << endl;
//        return 1;
//    }

//    if(options.subsample_grid_spacing > 0)
//        gridSimplify(mvs_data, options.subsample_grid_spacing);

//    exportOptions eo;
//    eo.sensor_vec = false;
//    eo.color = false;
//    eo.normals = true;
//    dir.suffix = "_normals";
////    if(dir.write_file.empty())
////        dir.write_file = dir.read_file;

//    exportPoints(dir,mvs_data,eo);
//}


int main(int argc, char const *argv[]){


    cliParser ip("normals");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getNormals())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----ESTIMATE NORMALS-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;

    importPLYPoints(ip.dh,data);

    if(ip.ro.normal_overwrite){
        data.has_normal = false;
        cout << "\nOverwriting existing normals..." << endl;
    }

    if(estimateNormals(data,ip.ro))
        return 1;

    exportOptions eo;
    eo.normals = true;
    eo.color = data.has_color;
    exportPLY(ip.dh,data.points,data.infos,eo);





    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----ESTIMATE NORMALS FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}


