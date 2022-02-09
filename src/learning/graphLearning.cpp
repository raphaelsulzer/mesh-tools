//#include "learning/graphLearning.h"


//void GraphLearning::Export(){

//    auto start = std::chrono::high_resolution_clock::now();
//    ofstream labels;
//    string out_path = GraphLearning.GraphLearningdir.path;
//    labels.open (out_path+"_labels.txt");


//    // n_edges in the graph = n_finite_facets + n_finite_facets / 4
//    // n_finite_facets / 4 = n_finite_cells -> n_edges = n_finite_facets + n_finite_cells
//    // actually not n_finite_facets but n_triangles, which is

//    int nc = Dt.number_of_cells();
//    Delaunay::All_cells_iterator fci;
//    labels << "# " << nc << " global_index a b c d n_inside_rays n_outside_rays dist_inside dist_outside "
//                              "// label_inside_outside" << endl;
//    // the traversed label is decided as follows:"
//    // 0: no intersection of segment+triangle or ray+triangle
//    // 1: intersection of segment+triangle, i.e. an outside traversal
//    // 2: intersection of ray+triangle, i.e. an inside traversal
//    // the inside_outside label is decided as follows:"
//    // 0: triangle is inside
//    // 1: triangle is outside
//    labels << setprecision(5);

//    // make adjacency files
//    ofstream adjacency_ij;
//    ofstream adjacency_ji;
//    adjacency_ij.open(out_path+"_adjacency_ij.txt");
//    adjacency_ji.open(out_path+"_adjacency_ji.txt");
//    adjacency_ij << "list of edges ij" << endl;
//    adjacency_ji << "list of edges ji" << endl;
//    // make edge_feature file
//    ofstream edge_features;
//    edge_features.open(out_path+"_edge_features.txt");
//    edge_features << "list of edge features" << endl;

//    int edge_count = 0;
//    int i, j;
//    Facet fac;
//    for(fci = Dt.all_cells_begin(); fci != Dt.all_cells_end(); fci++){

//        i = fci->info().global_idx;
//        labels << i << " ";
//        // iterate over the 4 triangles of the tetrahedron
//        for(int c = 0; c < 4; c++){
//            /// print the cell coords
//            labels << fci->vertex(c)->point() << " ";

//            // set the edge to the mirror cell
//            j = fci->neighbor(c)->info().global_idx;
//            if(j <= i) // skip doubled sided (done just below) and self-loops (added in pytorch)
//                continue;
//            adjacency_ij << i << " " << j << " ";
//            adjacency_ji << j << " " << i << " ";
//            fac = make_pair(fci,c);
//            if(Dt.is_infinite(fac)){
//                edge_features << edge_count++ << " 0.0 " << endl;
//                edge_features << edge_count++ << " 0.0 " << endl;
//            }
//            else{ // write two times the same features, because it will be used for two half-edges
//                edge_features << edge_count++ << " " << sqrt(Dt.triangle(fac).squared_area()) << " " << endl;
//                edge_features << edge_count++ << " " << sqrt(Dt.triangle(fac).squared_area()) << " " << endl;
//            }
//        }
//        // print the inside/outside distance and count
//        labels << fci->info().inside_count << " ";
//        labels << fci->info().outside_count << " ";
//        labels << fci->info().inside_dist << " ";
//        labels << fci->info().outside_dist << " ";

//        // print the two labels
//        labels << fci->info().gc_label << endl;

//    }
//    labels.close();
//    adjacency_ij.close();
//    adjacency_ji.close();
//    edge_features.close();


//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    std::cout << "Exported " << nc << " tetrahedron labels and" << std::endl;
//    std::cout << "2 x " << edge_count << " edges to " << out_path << " in " << duration.count() << "s\n" << std::endl;
//}







//GraphLearning::GraphLearning(const GraphLearningOptions& options,
//                             const Delaunay& Dt)
//{

//}
