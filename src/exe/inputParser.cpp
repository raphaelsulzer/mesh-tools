#include <util/args.hxx>
#include <util/helper.h>
#include <exe/inputParser.h>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;
#include <boost/filesystem.hpp>

int cliParser::parse(int argc, char const* argv[]){


    try{

        po::options_description all("");

        auto ioptions = initInput();
        all.add(ioptions);
        auto ooptions = initOutput();
        all.add(ooptions);

        if(mode == "labatut"){
            auto options = initLabatut();
            all.add(options);
        }
        if(mode == "occ2mesh"){
            auto options = initOcc2Mesh();
            all.add(options);
        }
        if(mode == "vsa"){
            auto options = initVsa();
            all.add(options);
        }
        if(mode == "collapse"){
            auto options = initCollapse();
            all.add(options);
        }
        if(mode == "simplify"){
            auto options = initSimplify();
            all.add(options);
        }
        if(mode == "eval"){
            auto options = initEval();
            all.add(options);
        }
        if(mode == "sample"){
            auto options = initSample();
            all.add(options);
        }
        if(mode == "normals"){
            auto options = initNormals();
            all.add(options);
        }
        if(mode == "eth"){
            auto options = initETH();
            all.add(options);
        }
        if(mode == "scan"){
            auto options = initScan();
            all.add(options);
        }
        if(mode == "iso"){
            auto options = initIso();
            all.add(options);
        }
        if(mode == "feat"){
            auto options = initFeat();
            all.add(options);
        }
        // parse all options to a single options map called vm
        po::store(po::parse_command_line(argc, argv, all), vm);

        if (vm.count("help")){
            cout << all << "\n";
            return 1;
        }

        // There must be an easy way to handle the relationship between the
        // option "help" and "host"-"port"-"config"
        // Yes, the magic is putting the po::notify after "help" option check
        po::notify(vm);
    }
    catch(std::exception& e){
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    catch(...){
        std::cerr << "Unknown error when parsing command line options!" << "\n";
        return 1;
    }

    return 0;
}


po::options_description cliParser::initInput(){

    po::options_description input_options("\nINPUT OPTIONS",description_width);
    input_options.add_options()
            ("help,h", "Help message")
            ("working_dir,w", po::value<string>()->required(), "Working directory.\nAll paths will be treated relative to this directory.")
            ("input_file,i", po::value<string>()->required(), "Input file")
            ("output_file,o", po::value<string>(), "Output file")
//            ("source,s", po::value<string>()->default_value("ply"), "Data source and options:"
//                                                            "\n\t-ply"
//                                                            "\n\t-npz"
//                                                            "\n\t-colmap"
//                                                            "\n\t-omvs"
//                                                            "\n\t-scan,#points,#cameras,std_noise,%outliers"
//                                                            "\n\t-tt,#scans"
//                                                            "\n\t-eth")
            ("source,s", po::value<string>()->default_value("ply"), "Data source:"
                                                            "\n\t-ply"
                                                            "\n\t-npz"
                                                            "\n\t-omvs (an OpenMVS project file)")
            ("groundtruth_file,g", po::value<string>(), "Groundtruth file")
            ("occ", po::value<string>(), "Occupancy file")
            ("transformation_file,t", po::value<string>(), "Transformation file")
            ("crop_file,c", po::value<string>(), "Crop file")
            ("scale", po::value<double>()->default_value(0.0), "Scale mean edge length of Delaunay to this value.")
            ("adt", po::value<double>()->default_value(-1), "Epsilon for adaptive 3DT.")
            ("icomp", po::value<int>()->default_value(1), "Number of connected components of input to keep. -1 = all.")
            ("iclose", po::value<int>()->default_value(1), "Try to close input open meshes with hole filling.")
        ;
    return input_options;
}
int cliParser::getInput(){

    /////////////// PATH & FILE INPUT ARGS + OUTPUT ARGS ///////////////

    /// required
    dh.path =  vm["working_dir"].as<string>();
    if(dh.path[dh.path.length()-1] != '/')
        dh.path+="/";
    dh.read_file = vm["input_file"].as<string>();
    boost::filesystem::path rf(dh.read_file);
    if(rf.has_extension()){
        dh.read_file_type = rf.extension().string();
        dh.read_file = dh.read_file.substr(0,dh.read_file.length() - 4);
    }

    /// data source
    vector<string> data_source;
    splitString(vm["source"].as<string>(), data_source, ',');
    ro.data_source = data_source[0];
    if(ro.data_source == "scan"){
        if(data_source.size() < 5){
            cout << "scan requires options #points,#cameras,std_noise,%outliers" << endl;
            return 1;
        }
        ro.number_of_points_to_scan = stoi(data_source[1]);
        ro.number_of_cameras = stoi(data_source[2]);
        ro.noise_std = stod(data_source[3]);
        ro.percentage_of_outliers = stod(data_source[4]);
    }
    else if(ro.data_source == "ply"){;}
    else if(ro.data_source == "npz"){;}
    else if(ro.data_source == "colmap"){;} // not really implemented anymore
    else if(ro.data_source == "omvs"){;}
    else if(ro.data_source == "eth"){;} // not really implemented anymore
    else if(ro.data_source == "tt"){  // not really implemented anymore
        if(data_source.size() < 5){
            cout << "scan requires options #points,#cameras,std_noise,%outliers" << endl;
            return 1;
        }
        ro.number_of_points_to_scan = stoi(data_source[1]);
        ro.number_of_cameras = stoi(data_source[2]);
        ro.noise_std = stod(data_source[3]);
        ro.percentage_of_outliers = stod(data_source[4]);
    }
    else{
        cerr << "\nNOT A VALID DATA SOURCE (-s).\n" << endl;
        cerr << "\nto see available data sources type sure --help\n" << endl;
        return 1;
    }

    /// optional
    if(vm.count("output_file"))
        dh.write_file =  vm["output_file"].as<string>();
    else
        dh.write_file = dh.read_file;
    if(vm.count("transformation_file"))
        dh.transformation_file = vm["transformation_file"].as<string>();
    if(vm.count("crop_file"))
        dh.crop_file = vm["crop_file"].as<string>();
    if(vm.count("groundtruth_file"))
        dh.gt_poly_file =  vm["groundtruth_file"].as<string>();
    if(vm.count("occ"))
        dh.occ_file =  vm["occ"].as<string>();

    // Delaunay
    if(vm.count("adt"))
        ro.Dt_epsilon = vm["adt"].as<double>();

    if(vm.count("iclose"))
        ro.try_to_close = vm["iclose"].as<int>();
    if(vm.count("icomp"))
        ro.number_of_components_to_keep = vm["icomp"].as<int>();

    if(vm.count("scale"))
        ro.scale = vm["scale"].as<double>();

    return 0;

}

po::options_description cliParser::initOutput(){


    ////////////////// EXPORT OPTIONS //////////////////

    po::options_description options("\nOUTPUT OPTIONS",description_width);
    options.add_options()
//        ("help", "produce help message")
            ("output_options,e", po::value<string>()->default_value("i"),"Specify export options as list of letters (without seperator):"
                                                                          "\n\t-n = normals"
                                                                          "\n\t-r = rgb color"
                                                                          "\n\t-v = sensor vector"
                                                                          "\n\t-p = sensor position"
                                                                          "\n-k = cameras"
                                                                          "\n-i = interface"
                                                                          "\n-z = isosurface"
                                                                          "\n-x = scan"
                                                                          "\n-m = mesh"
                                                                          "\n-s = sampling"
                                                                          "\n-h = convex Hull"
                                                                          "\n-f = colored Facets"
                                                                          "\n-c = cell Score")
            ("output_sampling", po::value<string>(), "Sample points on output mesh"
                                                                           "\n-gs,X = grid sampling + grid size"
                                                                           "\n-ps,X = point sampling + #points"
                                                                           "\n-as,X = point sampling + #points per area")

        ;


    return options;

}

int cliParser::getOutput(){
    vector<string> sampling;
    if(vm.count("output_sampling")){
        splitString(vm["output_sampling"].as<string>(), sampling, ',');
        eo.sampling_method = sampling[0];
        eo.sampling_method_option = stod(sampling[1]);
    }
    if(vm["output_options"].as<string>() == "all"){
        eo.normals = false;
        eo.color = true;
        eo.sensor_vec = true;
        eo.sensor_position = false;
        eo.convexHull = true;
        eo.sampling = true;
        eo.scan = true;
        eo.cameras = true;
        eo.mesh = true;
        eo.coloredFacets = true;
        eo.cellScore = true;
        eo.interface = true;
        eo.isosurface = true;
    }
    else{
        if (vm["output_options"].as<string>().find('n') != std::string::npos)
            eo.normals = true;
        if (vm["output_options"].as<string>().find('r') != std::string::npos)
            eo.color = true;
        if (vm["output_options"].as<string>().find('v') != std::string::npos)
            eo.sensor_vec = true;
        if (vm["output_options"].as<string>().find('p') != std::string::npos)
            eo.sensor_position = true;

        if (vm["output_options"].as<string>().find('k') != std::string::npos)
            eo.cameras = true;
        if (vm["output_options"].as<string>().find('i') == std::string::npos)
            eo.interface = false;
        if (vm["output_options"].as<string>().find('z') != std::string::npos)
            eo.isosurface = true;
        if (vm["output_options"].as<string>().find('x') != std::string::npos)
            eo.scan = true;
        if (vm["output_options"].as<string>().find('m') != std::string::npos)
            eo.mesh = true;
        if (vm["output_options"].as<string>().find('s') != std::string::npos)
            eo.sampling = true;
        if (vm["output_options"].as<string>().find('h') != std::string::npos)
            eo.convexHull = true;
        if (vm["output_options"].as<string>().find('f') != std::string::npos)
            eo.coloredFacets = true;
        if (vm["output_options"].as<string>().find('c') != std::string::npos)
            eo.cellScore = true;
    }
    return 0;
}

po::options_description cliParser::initOcc2Mesh(){



    po::options_description options("\nOCC2MESH OPTIONS",description_width);
    options.add_options()
            ("prediction_file,p", po::value<string>(), "Prediction file")
            ("occupancy_type", po::value<string>()->default_value("lo"), "Prediction type := [so(ftmax), si(gmoid), lo(git)].")
            ("gco", po::value<string>(), "Graph-cut optimization:"
             "\nBinary Type-Weight[;Type2-Weight2;Type3-Weight3]"
             "\nwith Type=[area,angle,cc] and Weight=(0,inf)")

        ;

    return options;
}

int cliParser::getOcc2Mesh(){

    // reconstruction and optimization
    dh.prediction_file = vm["prediction_file"].as<string>();
    boost::filesystem::path rf(dh.prediction_file);
    if(rf.has_extension()){
        dh.read_file_type = rf.extension().string();
        dh.prediction_file = dh.prediction_file.substr(0,dh.prediction_file.length() - 4);
    }
    ro.occupancy_type = vm["occupancy_type"].as<string>();

    vector<string> optimization;
    if(vm.count("gco")){
        splitString(vm["gco"].as<string>(), optimization, ',');
        ro.optimization = 1;
        for(int i = 0; i < optimization.size(); i++){
            vector<string> current_reg;
            splitString(optimization[i],current_reg,'-');
            if(current_reg.size()!=2){
                cout<<"\nERROR: not a valid regularization. Enter e.g. --gco area-0.1 ."<<endl;
                return 1;
            }
            if(current_reg[0] == "area")
                ro.area_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "angle")
                ro.angle_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "cc")
                ro.cc_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "vol")
                    ro.vol_reg = 1;
            else if(current_reg[0] == "sv")
                    ro.sv_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "ob")
                    ro.ob_reg_weight = stod(current_reg[1]);
            else{
                cout << "ERROR: " << current_reg[0] << " is not a valid regularization term!" << endl;
                return 1;
            }

        }
    }
    return 0;
}



po::options_description cliParser::initLabatut(){



    po::options_description options("\nLABATUT OPTIONS",description_width);
    options.add_options()
            ("cameras,c", po::value<int>()->default_value(1), "cameras=[1,inf). Number of cameras for ray tracing.")
            ("sigma", po::value<double>()->default_value(-1), "sigma=[-1,inf). -1 => mean noise with PCA.")
            ("alpha", po::value<double>()->default_value(32), "alpha=(0,inf).")
            ("gco", po::value<string>(), "Graph-cut optimization:"
             "\nBinary Type-Weight[;Type2-Weight2;Type3-Weight3]"
             "\nwith Type=[area,angle,cc] and Weight=(0,inf)")

            //// OLD:
//            ("smooth", po::value<int>()->default_value(0), "Smoothend the piecewise constant function.")
//            ("iso", po::value<string>()->default_value("30.0,0.1,0.1"), "Number of connected components of output to keep. Default = all.")
//            ("ocomp", po::value<int>()->default_value(0), "Number of connected components of output to keep. Default = all.")
//            ("omanifold", po::value<int>()->default_value(0), "Try to make output mesh manifold.")
//            ("oclose", po::value<int>()->default_value(0), "Try to close output open meshes with hole filling.")
//            ("clean", po::value<int>()->default_value(0), "Apply OpenMVS mesh cleaning.")
//            ("simplify", po::value<int>()->default_value(0), "Simplify mesh")
//            ("eval", po::value<int>()->default_value(0), "Evaluate mesh.")
        ;

    return options;
}

int cliParser::getLabatut(){



    // reconstruction and optimization
    ro.scoring = "rt";
    ro.occupancy_type = "labatut";
    ro.number_of_rays =  vm["cameras"].as<int>();
    ro.labatut_sigma = vm["sigma"].as<double>();
    ro.labatut_alpha = vm["alpha"].as<double>();


//    }
//    else if(ro.scoring == "cl"){
//        if(scoring.size() < 2){
//            cout << "ERROR: missing option in scoring" << endl;
//            cout << "put e.g. cl,sm" << endl;
//            return 1;

//        }
//    }
//    else if(ro.scoring == "clcs"){
//        if(scoring.size() < 3){
//            cout << "ERROR: missing option in scoring" << endl;
//            cout << "put e.g. clcs,500,sm" << endl;
//            return 1;

//        }
//        ro.number_of_points_per_cell = stoi(scoring[2]);
//    }
//    else if(ro.scoring == "cs")
//        ro.number_of_points_per_cell = stoi(scoring[1]);
//    else if(ro.scoring == "csrt"){
//        ro.number_of_points_per_cell = stoi(scoring[1]);
//        if(scoring.size() < 6){
//            cout << "ERROR: scoring method rt requires number_of_rays, type and labatut_sigma" << endl;
//            return 1;
//        }
//        ro.number_of_points_per_cell = stoi(scoring[1]);
//        ro.score_type = scoring[2];
//        ro.number_of_rays = stoi(scoring[3]);
//        ro.labatut_sigma = stod(scoring[4]);
//        ro.labatut_alpha = stod(scoring[5]);
//    }
//    else if(ro.scoring == "lrtcs")
//        ro.number_of_points_per_cell = stoi(scoring[1]);
//    else if(ro.scoring == "lrt")
//        ro.number_of_rays = stoi(scoring[1]);
//    else if(ro.scoring == "rt"){
//        if(scoring.size() < 5){
//            cout << "ERROR: scoring method rt requires type (=labatut), number_of_rays, labatut_sigma, labatut_alpha" << endl;
//            return 1;
//        }
//        ro.score_type = scoring[1];
//        ro.number_of_rays = stoi(scoring[2]);
//        ro.labatut_sigma = stod(scoring[3]);
//        ro.labatut_alpha = stod(scoring[4]);
//    }
//    else{
//        cerr << "\nNOT A VALID SCORING TYPE.\n" << endl;
//        cerr << "\nto see available scoring types type sure --help\n" << endl;
//        return 1;
//    }

    vector<string> optimization;
    if(vm.count("gco")){
        splitString(vm["gco"].as<string>(), optimization, ',');
        ro.optimization = 1;
        for(int i = 0; i < optimization.size(); i++){
            vector<string> current_reg;
            splitString(optimization[i],current_reg,'-');
            if(current_reg.size()!=2){
                cout<<"\nERROR: not a valid regularization. Enter e.g. --gco area-0.1 ."<<endl;
                return 1;
            }
            if(current_reg[0] == "area")
                ro.area_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "angle")
                ro.angle_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "cc")
                ro.cc_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "vol")
                    ro.vol_reg = 1;
            else if(current_reg[0] == "sv")
                    ro.sv_reg_weight = stod(current_reg[1]);
            else if(current_reg[0] == "ob")
                    ro.ob_reg_weight = stod(current_reg[1]);
            else{
                cout << "ERROR: " << current_reg[0] << " is not a valid regularization term!" << endl;
                return 1;
            }

        }
    }
//    vector<string> iso_options;
//    splitString(vm["iso"].as<string>(), iso_options, ',');
//    if(iso_options.size() < 3){
//        cout << "ERROR: need all three iso options (angle,radius,distance) in comma seperated form." << endl;
//        return 1;
//    }
//    ro.options[0] = stod(iso_options[0]);
//    ro.options[1] = stod(iso_options[1]);
//    ro.options[2] = stod(iso_options[2]);

//    if(vm.count("smooth"))
//        ro.smooth_field = vm["smooth"].as<int>();
//    if(vm.count("clean"))
//        ro.clean_mesh = vm["clean"].as<int>();
//    if(vm.count("ocomp"))
//        ro.number_of_components_to_keep = vm["ocomp"].as<int>();
//    if(vm.count("oclose"))
//        ro.try_to_close = vm["oclose"].as<int>();
//    if(vm.count("omanifold"))
//        ro.try_to_make_manifold = vm["omanifold"].as<int>();
//    if(vm.count("simplify"))
//        ro.simplify = vm["simplify"].as<int>();
//    if(vm.count("eval"))
//        ro.evaluate_mesh = vm["eval"].as<int>();
//    if(vm.count("gclosed"))
//        ro.gt_isclosed = vm["gclosed"].as<int>();
    return 0;
}


po::options_description cliParser::initFeat(){
    ////////////////// FEATURE EXTRACTION OPTIONS /////////////////
    po::options_description options("\nFEATURE EXTRACTION OPTIONS",description_width);
    options.add_options()
            ("gclosed", po::value<int>()->default_value(1), "Is ground truth closed?")
            ("rays", po::value<int>()->default_value(1), "Number of rays to trace.")
            ("pcp", po::value<int>()->default_value(100), "Number of points to sample per cell.")
            ("export", po::value<string>()->default_value(""), "Export the point cloud to [ply,npz,all].")
        ;
    return options;
}
int cliParser::getFeat(){

    // Delaunay
    ro.number_of_rays = vm["rays"].as<int>();
    ro.number_of_points_per_cell = vm["pcp"].as<int>();
    if(vm.count("gclosed"))
        ro.gt_isclosed = vm["gclosed"].as<int>();

    if(vm.count("export")){
        string etype = vm["export"].as<string>();
        if(etype == "all"){
            eo.tonpz = true;
            eo.toply = true;
        }
        else if(etype == "ply")
            eo.toply = true;
        else if(etype == "npz")
            eo.tonpz = true;
    }
    return 0;

}


po::options_description cliParser::initNormals(){
//    ///////////////// CALCULATE NORMALS OPTIONS /////////////////
    ////////////////// COLLAPSE OPTIONS /////////////////
    po::options_description normal_options("\nNORMAL ESTIMATION OPTIONS",description_width);
    normal_options.add_options()
            ("method", po::value<string>()->default_value("jet"), "[pca, jet, vsm]")
            ("neighborhood", po::value<int>()->default_value(0), "Neighborhood size to consider. 0 = automatic")
            ("orient", po::value<int>()->default_value(2), "Orientation. 0 = no, 1 = sensor, 2 = mst.")
            ("overwrite", po::value<int>()->default_value(0), "Overwrite existing normals.")

        ;
    return normal_options;
}
int cliParser::getNormals(){

    ro.normal_neighborhood = vm["neighborhood"].as<int>();
    ro.normal_orientation = vm["orient"].as<int>();
    ro.normal_method = vm["method"].as<string>();
    ro.normal_overwrite = vm["overwrite"].as<int>();

    return 0;
}


po::options_description cliParser::initIso(){
    po::options_description options("\nISOSURFACE EXTRACTION OPTIONS",description_width);
    options.add_options()
            ("method", po::value<string>()->default_value("mit"), "[mit, boissonnat]")
            ("field", po::value<string>()->default_value("occ"), "[occ, sdf]")
            ("value", po::value<double>()->default_value(0), "Neighborhood size to consider. 0 = automatic")
            ("options", po::value<string>()->default_value("30.0,0.1,0.1"), "Number of connected components of output to keep.")
            ("add_points", po::value<int>()->default_value(0), "Number of additional points to add for MIT.")
        ;
    return options;
}
int cliParser::getIso(){

    ro.field = vm["field"].as<string>();
    ro.method = vm["method"].as<string>();
    ro.value = vm["value"].as<double>();
    ro.add_points = vm["add_points"].as<int>();

    vector<string> iso_options;
    splitString(vm["options"].as<string>(), iso_options, ',');
    if(iso_options.size() < 3){
        cout << "ERROR: need all three iso options (angle,radius,distance) in comma seperated form." << endl;
        return 1;
    }
    ro.options[0] = stod(iso_options[0]);
    ro.options[1] = stod(iso_options[1]);
    ro.options[2] = stod(iso_options[2]);

    return 0;
}


po::options_description cliParser::initCollapse(){

    ////////////////// COLLAPSE OPTIONS /////////////////
    po::options_description collapse_options("\nEDGE COLLAPSe OPTIONS",description_width);
    collapse_options.add_options()
            ("edges", po::value<double>()->required(), "Percentage of edges to keep")
            ("clean", po::value<int>(), "Apply OpenMVS mesh cleaning")
            ("eval", po::value<int>(), "Evaluate mesh")

        ;
    return collapse_options;
}
int cliParser::getCollapse(){

    ro.keep_edges = vm["edges"].as<double>();
    if(vm.count("clean"))
        ro.clean_mesh = vm["clean"].as<int>();
    if(vm.count("eval"))
        ro.evaluate_mesh = vm["eval"].as<int>();

    return 0;
}

po::options_description cliParser::initVsa(){
    ////////////////// VSA OPTIONS /////////////////
    po::options_description vsa_options("\nVSA OPTIONS",description_width);
    vsa_options.add_options()
            ("proxies", po::value<int>()->required(), "Max number of proxies")
            ("comp", po::value<int>(), "Number of connected components to keep.")
            ("close", po::value<int>(), "Try to close open meshes with hole filling.")
            ("clean", po::value<int>(), "Apply OpenMVS mesh cleaning")
            ("eval", po::value<int>(), "Evaluate mesh")

        ;
    return vsa_options;
}
int cliParser::getVsa(){
    ro.max_number_of_proxies = vm["proxies"].as<int>();
    if(vm.count("close"))
        ro.try_to_close = vm["close"].as<int>();
    if(vm.count("comp"))
        ro.number_of_components_to_keep = vm["comp"].as<int>();
    if(vm.count("clean"))
        ro.clean_mesh = vm["clean"].as<int>();
    if(vm.count("eval"))
        ro.evaluate_mesh = vm["eval"].as<int>();

    return 0;
}

po::options_description cliParser::initSimplify(){
    ////////////////// VSA OPTIONS /////////////////
    po::options_description ransac_options("\nSIMPLIFY OPTIONS",description_width);
    ransac_options.add_options()
            ("repsilon", po::value<double>()->default_value(0.0), "Max distance point - plane [0, +inf)")
            ("rnormal", po::value<double>()->default_value(0.9, "0.9"), "Max normal variation point - plane [0,1]")
            ("reg", po::value<int>()->default_value(1), "Optional plane regularization")
            ("strepsilon", po::value<double>()->default_value(0.0), "Epsilon for Point set structuring [0, +inf)")
            ("snepsilon", po::value<double>()->default_value(0.0), "Epsilon for Point snapping [0, +inf)")
            ("gco", po::value<string>(), "Graph-cut optimization."
             "\nSpecify Binary: Type1-Weight1[;Type2-Weight2;Type3-Weight3]"
             "\nwith Type=[area,angle,cc] and Weight=(0,inf)")
            ("simplify", po::value<int>()->default_value(0), "Simplify mesh")
            ("eval", po::value<int>(), "Evaluate mesh")
        ;
    return ransac_options;
}
int cliParser::getSimplify(){

    ro.ransac_epsilon = vm["repsilon"].as<double>();
    ro.ransac_normal = vm["rnormal"].as<double>();
    if(vm.count("reg"))
        ro.ransac_regularization = vm["reg"].as<int>();
    ro.structure_epsilon = vm["strepsilon"].as<double>();
    ro.snap_epsilon = vm["snepsilon"].as<double>();
    if(vm.count("simplify"))
        ro.simplify = vm["simplify"].as<int>();
    if(vm.count("eval"))
        ro.evaluate_mesh = vm["eval"].as<int>();

    vector<string> optimization;
    if(vm.count("gco")){
        ro.optimization = 1;
        splitString(vm["gco"].as<string>(), optimization, ',');
        ro.free_weight = stod(optimization[0]);
        ro.reg_weight = stod(optimization[1]);
    }

    return 0;
}


po::options_description cliParser::initEval(){
    po::options_description eval_options("\nEVAL OPTIONS",description_width);
    eval_options.add_options()
            ("ocomp", po::value<int>()->default_value(0), "Number of connected components of output to keep. Default = all.")
            ("omanifold", po::value<int>()->default_value(0), "Try to make output mesh manifold.")
            ("oclose", po::value<int>()->default_value(0), "Try to close output open meshes with hole filling.")
        ;
    return eval_options;
}
int cliParser::getEval(){
    if(vm.count("ocomp"))
        ro.number_of_components_to_keep = vm["ocomp"].as<int>();
    if(vm.count("oclose"))
        ro.try_to_close = vm["oclose"].as<int>();
    if(vm.count("omanifold"))
        ro.try_to_make_manifold = vm["omanifold"].as<int>();
    return 0;
}


po::options_description cliParser::initSample(){
    po::options_description sample_options("\nSAMPLE OPTIONS",description_width);
    sample_options.add_options()
            ("ocomp", po::value<int>()->default_value(0), "Number of connected components of output to keep. Default = all.")
            ("omanifold", po::value<int>()->default_value(0), "Try to make output mesh manifold.")
            ("oclose", po::value<int>()->default_value(0), "Try to close output open meshes with hole filling.")
        ;
    return sample_options;
}
int cliParser::getSample(){
    if(vm.count("ocomp"))
        ro.number_of_components_to_keep = vm["ocomp"].as<int>();
    if(vm.count("oclose"))
        ro.try_to_close = vm["oclose"].as<int>();
    if(vm.count("omanifold"))
        ro.try_to_make_manifold = vm["omanifold"].as<int>();
    return 0;
}

po::options_description cliParser::initScan(){
    po::options_description options("\nSCAN OPTIONS",description_width);
    options.add_options()
            ("gclosed", po::value<int>()->default_value(1), "Is ground truth closed?")
            ("points,p", po::value<int>()->default_value(500), "Number of scanning points.")
            ("noise,n", po::value<double>()->default_value(0.0), "Std of Gaussian noise to add.")
            ("outliers", po::value<double>()->default_value(0), "Percentage of outliers to add to scanning points.")
            ("cameras,c", po::value<int>()->default_value(3), "Number of cameras.")
            ("normal_method", po::value<string>(), "Method for normal estimation: [none (default), pca, jet, vsm]")
            ("normal_neighborhood", po::value<int>()->default_value(0), "Neighborhood size to consider. 0 = automatic (default)")
            ("normal_orient", po::value<int>()->default_value(2), "Orientation. 0 = no, 1 = sensor, 2 = mst (default).")
            ("export", po::value<string>()->default_value("ply"), "Export to [ply,npz,all].")
        ;
    return options;
}
int cliParser::getScan(){
    if(vm.count("gclosed"))
        ro.gt_isclosed = vm["gclosed"].as<int>();
    if(vm.count("cameras"))
        ro.number_of_cameras = vm["cameras"].as<int>();
    if(vm.count("points"))
        ro.number_of_points_to_scan = vm["points"].as<int>();
    if(vm.count("noise"))
        ro.noise_std = vm["noise"].as<double>();
    if(vm.count("outliers"))
        ro.percentage_of_outliers = vm["outliers"].as<double>();



    if(vm.count("normal_method")){
        ro.normal_method = vm["normal_method"].as<string>();
        ro.normal_neighborhood = vm["normal_neighborhood"].as<int>();
        ro.normal_orientation = vm["normal_orient"].as<int>();
    }

    if(vm.count("export")){
        string etype = vm["export"].as<string>();
        if(etype == "all"){
            eo.tonpz = true;
            eo.toply = true;
        }
        else if(etype == "ply")
            eo.toply = true;
        else if(etype == "npz")
            eo.tonpz = true;
        else{
            cerr << "\nERROR: not a valid exoport type, choose either." << endl;
            return 1;
        }
    }
    return 0;
}


po::options_description cliParser::initETH(){
    po::options_description options("\nETH OPTIONS",description_width);
    options.add_options()
            ("method", po::value<string>()->default_value("jet"), "[pca, jet (default), vsm]")
            ("neighborhood", po::value<int>()->default_value(0), "Neighborhood size to consider. 0 = automatic (default)")
            ("orient", po::value<int>()->default_value(2), "Orientation. 0 = no, 1 = sensor, 2 = mst (default).")
        ;
    return options;
}
int cliParser::getETH(){
    ro.normal_neighborhood = vm["neighborhood"].as<int>();
    ro.normal_orientation = vm["orient"].as<int>();
    ro.normal_method = vm["method"].as<string>();
    return 0;
}
