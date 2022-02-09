#include <exe/sure.h>
#include <exe/inputParser.h>

#include <base/cgal_typedefs.h>
#include <IO/fileIO.h>
#ifdef OpenMVS
#include <IO/mvIO.h>
#endif

#include <boost/filesystem.hpp>
using namespace boost::filesystem;


int main(int argc, char const *argv[]){


    cliParser ip("feat");
    if(ip.parse(argc, argv))
        return 1;
    if(ip.getInput())
        return 1;
    if(ip.getOutput())
        return 1;
    if(ip.getFeat())
        return 1;

    auto start = std::chrono::high_resolution_clock::now();
    cout << "\n-----FEATURE EXTRACTION-----" << endl;
    cout << "\nWorking dir set to:\n\t-" << ip.dh.path << endl;

    dataHolder data;
//    if(extractFeatures(ip.dh, data, ip.ro, ip.eo))
//        return 1;

    loadOMVSScene(ip.dh,data);

    if(ip.eo.cameras){
        ip.dh.suffix = "_cameras";
        exportCameraCenter(ip.dh, data);
    }
    // export the scanned points
    if(ip.eo.toply){
        ip.eo.color = data.has_color;
        ip.eo.normals = data.has_normal;
        exportPLY(ip.dh, data.points, data.infos, ip.eo);
    }

//    if(exportO.tonpz){
    toXTensor(data);
    exportNPZ(ip.dh,data);
//    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    cout << "\n-----FEATURE EXTRACTION FINISHED in "<< duration.count() << "s -----\n" << endl;

    return 0;

}



