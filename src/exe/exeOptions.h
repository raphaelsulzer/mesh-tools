#ifndef EXEOPTIONS_H
#define EXEOPTIONS_H

struct runningOptions{



    string data_source;
    string gt_scan_source;
    string scoring;

    int number_of_scans; // how many individual scans for the TanksAndTemples input

    // preprocess input
    double scale = 0.0;
    int calculateNoise = 0;
    int orderSensors = 0;
    double subsample_n_points = 0.0;
    double subsample_grid_spacing = 0.0;
    int gt_isclosed = 1;

    // rayTracingTet
    int number_of_rays = 1;
    string score_type = "";
    int tetTracing = 0;
    double labatut_sigma = -1;
    double labatut_alpha = 32;

    // Delaunay tiangulation
    double Dt_epsilon; // for adaptive Delaunay triangulation
    int insert_sensor = 0;

    int make_global_cell_idx = 0;
    int make_finite_cell_idx = 0;

    int make_global_vertex_idx = 0;
    int make_finite_vertex_idx = 0;

    // optimization
    int optimization = 0;
    int vol_reg = 0;
    double sv_reg_weight = 0.0;
    double ob_reg_weight = 0.0;
    double angle_reg_weight = 0.0;
    double area_reg_weight = 0.0;
    double cc_reg_weight = 0.0;
    int gco_iterations;

    // iso options
    string method;
    int smooth_field = 0;
    vector<double> options = vector<double>(3);
    string field;
    double value;
    int add_points;

    double read_regularization_weight;
    string read_scoring;

    // cleaning
    int clean_mesh;
    int number_of_components_to_keep = 0;
    int try_to_close = 0;
    int try_to_make_manifold = 0;
    int factor_for_removing_large_faces;
    int fix_nm_edges = 0;

    // evaluation
    bool evaluate_mesh = false;

    // learning scanning
    int number_of_cameras;
    int number_of_points_to_scan;
    double percentage_of_outliers = 0.0;
    double noise_std = 0.0;
    int number_of_points_per_cell;
    // learning import
    string prediction_type;
    // learning export
    int export_rays = 0;
    int export_noise;
    int ground_truth = 0;

    // normal calculation
    int normal_neighborhood;
    int normal_orientation;
    int normal_overwrite;
    string normal_method;

    // variational shape approximation
    int max_number_of_proxies;

    // edge collapse
    double keep_edges;

    // ransac
    double ransac_epsilon;
    double ransac_normal;
    int ransac_regularization;

    //structure
    double structure_epsilon;

    // snapping
    double snap_epsilon;
    double free_weight;
    double reg_weight;

    // simplify
    int simplify;




};

#endif // SURFACERECON_H
