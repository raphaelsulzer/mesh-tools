#include <base/cgal_typedefs.h>
#include <CGAL/IO/read_ply_points.h> // CGAL file that is not present in CGAL version < 4.11, so copied it to my IO folder

#include <IO/fileIO.h>
#include <IO/ttIO.h>
#include <util/geometricOperations.h>
#include <util/helper.h>
//#include <learning/learning.h>
#include <random>

#include <rPLY/ply.h>

using namespace std;

//////////////////////////////////////////////////////////
///////////////////// TANKS AND TEMPLES I/O ///////////////////////////
//////////////////////////////////////////////////////////

int importTanksAndTemplesIs(dirHolder dir, dataHolder& data, Point sensor_position)
{
    cout << "\t-read LiDAR points from " << dir.read_file << endl;

    std::string ifn = dir.path + dir.read_file;

    // read Binary PLY with sensor
    Mesh_ply aMesh;
    Import_PLY(ifn.c_str(), &aMesh);

    assert(aMesh.mVertices.size() > 0);

    data.has_sensor = true;
    data.has_normal = false;
    data.has_color = false;

    if(!data.has_sensor && !data.has_normal){
        cout << "\nWARNING: NEITHER NORMALS NOR SENSOR INFORMATION FOUND IN THE INPUT FILE!\n" << endl;
    }
    for(int i = 0; i < aMesh.mVertices.size(); i++){
        // parse points
        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
        data.points.push_back(pt);
        vertex_info vinf;

        // sensor
        if(data.has_sensor){ // if there is a sensor location saved under scalar_x0
            vinf.sensor_positions.push_back(sensor_position);
            Vector vec(vinf.sensor_positions[0].x() - pt.x(), vinf.sensor_positions[0].y() - pt.y(), vinf.sensor_positions[0].z() - pt.z());
            vinf.sensor_vec = vec;

        }

        // index, e.g. used in the concatenate data function
        vinf.global_idx = i;
        // save vertex info
        data.infos.push_back(vinf);
    }
    cout << "\t-" << data.points.size() << " points read" << endl;

    return 0;
}

int importTanksAndTemplesScannerLocations(dirHolder dir, dataHolder& data, runningOptions options){


    dir.path+="is/";
    dir.read_file = "scanner_pos.txt";
    ifstream file(dir.path+dir.read_file);

    cout << "\nLoad sensor position file..." << endl;
    cout << "\t-from " << dir.path+dir.read_file << endl;

    if(!file){
        cout << "\nERROR: FILE DOES NOT EXIST OR IS EMPTY!" << endl;
        return 1;
    }

    // load all cameras
    string line;
    vector<string> lineAsVector;
    vector<Point> ss;
    while(getline(file, line)){
        splitString(line, lineAsVector, ' ');
        ss.push_back(Point(stod(lineAsVector[1]),stod(lineAsVector[2]),stod(lineAsVector[3])));
    };

    // shuffle the cameras
    std::random_device rd;
    std::mt19937 g(rd());
    g.seed(42);

    std::shuffle(ss.begin(), ss.end(), g);

    for(int i = 0; i < options.number_of_cameras; i++)
        data.sampled_sensors.push_back(ss[i]);

    return 0;
}


int importTanksAndTemples(dirHolder dir, dataHolder& data, runningOptions options){

    auto start = std::chrono::high_resolution_clock::now();

    dir.path+="is/";
    dir.read_file = "scanner_pos.txt";
    ifstream file(dir.path+dir.read_file);

    cout << "\nLoad sensor position file..." << endl;
    cout << "\t-from " << dir.path+dir.read_file << endl;

    if(!file){
        cout << "\nERROR: FILE DOES NOT EXIST OR IS EMPTY!" << endl;
        return 1;
    }

    string line;
    vector<string> lineAsVector;
    for(int i = 0; i < options.number_of_scans; i++){
        if(!getline(file, line))
            break;
        dataHolder d;
        splitString(line, lineAsVector, ' ');
        dir.read_file = lineAsVector[0].c_str();
        importTanksAndTemplesIs(dir, d, Point(stod(lineAsVector[1]),stod(lineAsVector[2]),stod(lineAsVector[3])));
        data.points.insert(data.points.end(),d.points.begin(),d.points.end());
        data.infos.insert(data.infos.end(),d.infos.begin(),d.infos.end());
    };

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "\t-" << data.points.size() << " points read" << endl;
    cout << "\t-in " << duration.count() << "s" << endl;

    return 0;
}









