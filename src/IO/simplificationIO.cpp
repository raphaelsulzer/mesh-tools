#include <base/cgal_typedefs.h>

#include <IO/fileIO.h>
#include <IO/simplificationIO.h>
#include <util/geometricOperations.h>
#include <util/helper.h>
#include <CGAL/Point_set_3/IO.h>
#include <rPLY/ply.h>

using namespace std;


void exportStructured(dirHolder& dir, PSS& pss){

    //////////////////////////////
    ////// Export STRUCTURE //////
    //////////////////////////////
    cout << "\t-" << pss.size () << " structured points generated." << endl;
    vector<Point> pcorner, pedge, pplane, pfree, pall;
    vector<vertex_info> icorner, iedge, iplane, ifree, iall;

    for(int i = 0; i < pss.size(); i++){
        vector<EPICK::Plane_3> pi;
        pss.adjacency(i,back_inserter(pi));
        pall.push_back(pss.point(i));
        vertex_info inf;
        if(pi.size() == 0){
            pfree.push_back(pss.point(i));
            inf.color = CGAL::white();
            ifree.push_back(inf);
        }
        else if(pi.size() == 1){
            pplane.push_back(pss.point(i));
            inf.color = CGAL::blue();
            iplane.push_back(inf);
        }
        else if(pi.size() == 2){
            pedge.push_back(pss.point(i));
            inf.color = CGAL::red();
            iedge.push_back(inf);
        }
        else if(pi.size() >= 3){
            pcorner.push_back(pss.point(i));
            inf.color = CGAL::yellow();
            icorner.push_back(inf);
        }
        iall.push_back(inf);
    }
    exportOptions eo;
    eo.color = true;

//    dir.suffix = "_structured_corners";
//    exportPLY(dir,pcorner,icorner,eo);

//    dir.suffix = "_structured_edges";
//    exportPLY(dir,pedge,iedge,eo);

//    dir.suffix = "_structured_planes";
//    exportPLY(dir,pplane,iplane,eo);

//    dir.suffix = "_structured_free";
//    exportPLY(dir,pfree,ifree,eo);

    dir.suffix = "_structured";
    exportPLY(dir,pall,iall,eo);

}

void exportAlphaShape(dirHolder& dir, vector<Point>& points){


    ofstream out(dir.path+dir.write_file+dir.suffix+".off");
    out << "OFF\n" << points.size() << " 1 0\n";

    stringstream face;
    face << points.size();

    for(int i = 0; i < points.size(); i++){
        out << points[i] << endl;
        face << " " << i;
    }
    out << face.str();

}

void exportPolygon(dirHolder& dir, Polygon_2& poly, EPICK::Plane_3& plane){


    vector<Segment> segs;
    for(int i = 0; i < poly.size()-1; i++){
        Segment seg(plane.to_3d(poly.vertex(i)),plane.to_3d(poly.vertex(i+1)));
        segs.push_back(seg);
    }
    exportEdges(dir,segs);
}

void exportEdges(dirHolder& dir, vector<Segment>& segs){

    bool debug = 0;

    cout << "\nExport edges..." << endl;
    cout << "\t-to " << dir.write_file+dir.suffix+".ply" << endl;

    stringstream pts, ind, header;

    int i = 0;
    for(const auto s : segs){
        // points
        pts << s.start() << " 255 0 0" << endl;
        if(!debug)
            pts << s.end() << " 255 0 0" << endl;
        // edge indices
        if(debug)
            ind << i << " " << i+1 << endl;
        else
            ind << 2*i << " " << 2*i+1 << endl;
        i++;
    }
//    pts << segs[i-1].end() << " 0 255 0" << endl;

    // header
    header << "ply" << endl;
    header << "format ascii 1.0" << endl;
    if(debug)
        header << "element vertex " << i+2 << endl;
    else
        header << "element vertex " << 2*i << endl;
    header << "property float x" << endl;
    header << "property float y" << endl;
    header << "property float z" << endl;
    header << "property uchar red" << endl;
    header << "property uchar green" << endl;
    header << "property uchar blue" << endl;
    header << "element edge " << segs.size() << endl;
    header << "property int vertex1" << endl;
    header << "property int vertex2" << endl;
    header << "end_header" << endl;


    ofstream out(dir.path+dir.write_file+dir.suffix+".ply");
    out << header.str();
    out << pts.str();
    out << ind.str();
}


void exportConvexHull2d(dirHolder& dir,DT2& as,CGAL::Color& col){

    // idea of exporting the points in the order of which they appear
    // on the facet_iterator is comming from here: https://stackoverflow.com/a/15917845
    // stupid thing about this approach is that lots of dublicate points (not vertices) are saved
    // so probably still better to loop, index and export points first, and then loop over facets.

    stringstream pts, ind;
    int i=0;
    auto fit = as.finite_faces_begin();
    while(fit != as.finite_faces_end()){

        /// if you care about orientation of the facets
//      int indices[3]={
//        (facets[i].second+1)%4,
//        (facets[i].second+2)%4,
//        (facets[i].second+3)%4,
//      };
      /// according to the encoding of vertex indices, this is needed to get
      /// a consistent orienation
//      if ( facets[i].second%2==0 ) std::swap(indices[0], indices[1]);

        pts <<
        fit->vertex(0)->info() << "\n" <<
        fit->vertex(1)->info() << "\n" <<
        fit->vertex(2)->info() << "\n";
        ind <<
        "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << " " <<
        to_string(col.r()) << " " << to_string(col.g()) << " " << to_string(col.b()) <<  "\n";

        i++, fit++;
    }

    ofstream out(dir.path+dir.write_file+dir.suffix+".off");
    out << "OFF "<< 3*i << " " << i << " 0\n";
    out << pts.str();
    out << ind.str();

}

