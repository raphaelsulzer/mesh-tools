#pragma once

#include "base/cgal_typedefs.h"
#include <IO/fileIO.h>
#include <util/args.hxx>
#include <util/helper.h>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid.hpp>

using namespace std;
namespace po = boost::program_options;

struct StringReader
{
    void operator()(const string& name, const string& value, vector<string>& valueAsVector){

        splitString(value, valueAsVector, ',');
    }
};

class cliParser{

public:

    po::variables_map vm; // this is a std::map which will be populate with all options available for the current mode
    exportOptions eo;
    runningOptions ro;
    dirHolder dh;
    string mode;
    int description_width = 160;

    cliParser(string mode_):
        mode(mode_)
    {}

    int parse(int argc, char const *argv[]); // parse the input args and populate vm map

    po::options_description initInput();
    po::options_description initOutput();

    po::options_description initLabatut();
    po::options_description initOcc2Mesh();
    po::options_description initNormals();
    po::options_description initCollapse();
    po::options_description initVsa();
    po::options_description initSimplify();
    po::options_description initEval();
    po::options_description initSample();
    po::options_description initETH();
    po::options_description initScan();
    po::options_description initIso();
    po::options_description initFeat();

    int getInput(); // populate dh (dirHolder)
    int getOutput(); // populate eo (exportOptions)

    int getLabatut(); // populate ro (runningOptions)
    int getOcc2Mesh(); // populate ro (runningOptions)
    int getVsa(); // populate ro (runningOptions)
    int getSimplify(); // populate ro (runningOptions)
    int getCollapse(); // populate ro (runningOptions)
    int getEval(); // populate ro (runningOptions)
    int getSample(); // populate ro (runningOptions)
    int getNormals(); // populate ro (runningOptions)
    int getETH(); // populate ro (runningOptions)
    int getScan(); // populate ro (runningOptions)
    int getIso(); // populate ro (runningOptions)
    int getFeat(); // populate ro (runningOptions)

};
