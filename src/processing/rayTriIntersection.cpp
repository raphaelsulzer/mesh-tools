#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/make_surface_mesh.h>

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_3<Kernel>              Delaunay;
typedef Delaunay::Cell_handle                               Cell_handle;
typedef Delaunay::Vertex_handle                             Vertex_handle;
typedef Delaunay::Locate_type                               Locate_type;
typedef Delaunay::Point_3                                     Point;
typedef Delaunay::Facet                                     Facet;

typedef Kernel::Vector_3 Vector;

typedef std::pair<Point, Vector> PointVectorPair;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


std::vector<PointVectorPair> estimateNormalsFun(const std::vector<Point> points)
{

    std::vector<PointVectorPair> pointVectorPairs(points.size());

    for(std::size_t i=0; i < points.size(); ++i)
    {
//        std::cout << "here" << std::endl;
        pointVectorPairs[i].first = points[i];
    }
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>());

    const int nb_neighbors = 3; // K-nearest neighbors = 3 rings
    CGAL::pca_estimate_normals<Concurrency_tag>
      (pointVectorPairs, nb_neighbors,
       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
       normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    for(std::size_t i=0; i < points.size(); ++i)
    {
        std::cout << pointVectorPairs[i].second << std::endl;
    }

    return pointVectorPairs;
};


int rayTriIntersectionFun(const char* ofn)
{

    std::ofstream fo;
    fo.open (ofn);


//    std::list<Point> L;
//    L.push_front(Point(0,2,2));
//    L.push_front(Point(0,7,2));
//    L.push_front(Point(0,7,4));
//    L.push_front(Point(0,4,4));
//    L.push_front(Point(0,4,7));
//    L.push_front(Point(0,2,7));

//    L.push_front(Point(3,2,2));
//    L.push_front(Point(3,7,2));
//    L.push_front(Point(3,7,4));
//    L.push_front(Point(3,4,4));
//    L.push_front(Point(3,4,7));
//    L.push_front(Point(3,2,7));

    std::vector<Point> L(12);
    L[0]=Point(0,2,2);
    L[1]=Point(0,7,2);
    L[2]=Point(0,7,4);
    L[3]=Point(0,4,4);
    L[4]=Point(0,4,7);
    L[5]=Point(0,2,7);

    L[6]=Point(3,2,2);
    L[7]=Point(3,7,2);
    L[8]=Point(3,7,4);
    L[9]=Point(3,4,4);
    L[10]=Point(3,4,7);
    L[11]=Point(3,2,7);

    std::vector<PointVectorPair> pVP = estimateNormalsFun(L);

    Delaunay Dt(L.begin(), L.end());

//    Delaunay Dt(pVP.begin(), pVP.end());

//    //    Triangulation T(L.begin(), L.end());
//    Delaunay::size_type nv = Dt.number_of_vertices();
//    Delaunay::size_type nf = Dt.number_of_finite_facets();

//    std::cout << "Triangulation done.." << std::endl;

//    // start writing file
//    fo << "ply" << std::endl;
//    fo << "format ascii 1.0" << std::endl;
//    fo << "comment VCGLIB generated" << std::endl;
//    fo << "element vertex " << nv << std::endl;
//    fo << "property float x" << std::endl;
//    fo << "property float y" << std::endl;
//    fo << "property float z" << std::endl;
//    fo << "property float nx" << std::endl;
//    fo << "property float ny" << std::endl;
//    fo << "property float nz" << std::endl;
//    fo << "element face " << nf << std::endl;
//    fo << "property list uchar int vertex_indices" << std::endl;
//    fo << "end_header" << std::endl;
//    fo << std::setprecision(12);



//    // give every vertex from the triangulation an index starting at 0
//    // and already print the point coordinates, color and normal of the vertex to the PLY file
//    std::map<Vertex_handle, int> Vertices;
//    int index = 0;
//    Delaunay::Finite_vertices_iterator vft;
//    for (vft = Dt.finite_vertices_begin(); vft != Dt.finite_vertices_end(); vft++){
//        // assign index to vertex handle. this is needed later for finding the index a certain cell is constructed with
//        // TODO: maybe this could already be included in the info vector of every point/vertex, and then instead of using the finite_vertices_iterator,
//        // one could iterate over the index vector, so that the point order of the file is determined by the info vector
//        // and you would do the find operation to find the point in a std::map<Point, idx>
//        Vertices[vft] = index;
//        // print data to file
//        fo << vft->point() << std::endl;
////        fo << vft->point() << " "                           // coordinates
////           << int(std::get<1>(vft->info())[0]) << " "       // red
////           << int(std::get<1>(vft->info())[1]) << " "       // green
////           << int(std::get<1>(vft->info())[2]) << " "       // blue
////           << std::get<0>(vft->info()) << std::endl;        // normal
//        index++;
//    }

//    // Save the facets to the PLY file
//    //std::cout << "iterate over finite triangles: " << std::endl;
//    Delaunay::Finite_facets_iterator fft;
//    int vidx;
//    // initialise cell and vertex handle
//    Cell_handle c;
//    Vertex_handle v;
//    for(fft = Dt.finite_facets_begin() ; fft != Dt.finite_facets_end() ; fft++){

//        // facet fft is represented by std::pair(cell c, int vidx). vidx is the vertex opposite to the cell.
//        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite.
//        // see: https://doc.cgal.org/latest/Triangulation_3/index.html


//        c = fft->first;         // cell
//        vidx = fft->second;     // vertex index


//        //std::cout << "is infinite: " << T.is_infinite(v) << std::endl;
//        //std::cout << "opposite vertex: " << vidx << std::endl;
//        fo << 3 << ' ';
//        // if opposite vertex vidx is 2, we start at j = vidx + 1 = 3, 3%4 = 3
//        // next iteration: j = 4, 4%4 = 0, next iteration: j = 5, 5%4 = 1;
//        // so we exactely skip 2 - the opposite vertex.
//        for(int j = vidx + 1 ; j <= vidx + 3 ; j++){


//            //std::cout << "modulo: " << j%4;

//            // in the first and second iteration I am calling vertex() with the same value for j%4,
//            // but get a different vertex idx. Reason is that the cell c is different (see std::cout of cell points).
//            // for some reason the first facet is seen from the infinite cell.
//            v = c->vertex(j%4);

//            // print the indicies of each cell to the file
//            fo << Vertices.find(v)->second << ' ';

//        }
//        // new cell
//        fo << std::endl;
//    }

    fo.close();

    return 0;


}




