#include <base/cgal_typedefs.h>
#include <util/helper.h>
//#include <processing/tetIntersection.h>

#include <random>
#include <iterator>

using namespace std;

#ifdef RECONBENCH
#include <modeling/Vector3.h>
const Vector3RB point2Vector3RB(const Point p){

    return Vector3RB(p.x(),p.y(),p.z());

}
#endif

const Eigen::Vector3d cgal2Eigen(const Point p){

    return Eigen::Vector3d(p.x(), p.y(), p.z());
}

const Point eigen2Cgal(const Eigen::Vector3f p){

    return Point(p.x(), p.y(), p.z());
}

vector<Eigen::Vector3d> cgal2Eigen(vector<Point> in){

    vector<Eigen::Vector3d> out;
    for(const auto& point : in){
        Eigen::Vector3d p(point.x(), point.y(), point.z());
        out.push_back(p);
    }
    return out;
}

vector<Point> eigen2Cgal(const vector<Eigen::Vector3d> in){

    vector<Point> out;
    for(const auto& point : in){
        Point p(point.x(), point.y(), point.z());
        out.push_back(p);
    }
    return out;
}


#include <CGAL/Polyhedron_copy_3.h>
Polyhedron_Exact inexact2exactPolyhedron(Polyhedron in){


    CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> modifier(in);
    Polyhedron_Exact ex;
    ex.delegate(modifier);
    return ex;

}

Polyhedron exact2inexactPolyhedron(Polyhedron_Exact ex){


    CGAL::Polyhedron_copy_3<Polyhedron_Exact, Polyhedron::HalfedgeDS> modifier(ex);
    Polyhedron in;
    in.delegate(modifier);
    return in;

}


size_t splitString(const std::string &txt, std::vector<std::string> &strs, char ch)
{
    size_t pos = txt.find( ch );
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while( pos != std::string::npos ) {
        strs.push_back( txt.substr( initialPos, pos - initialPos ) );
        initialPos = pos + 1;

        pos = txt.find( ch, initialPos );
    }

    // Add the last one
    strs.push_back( txt.substr( initialPos, std::min( pos, txt.size() ) - initialPos + 1 ) );

    return strs.size();
}

string double2string(double regularization_weight){
    ostringstream out;
    if(regularization_weight == 0.0)
        out.precision(0);
    else if(regularization_weight >= 10)
        out.precision(0);
    else if(regularization_weight >= 1)
        out.precision(1);
    else if(regularization_weight >= 0.1)
        out.precision(1);
    else if(regularization_weight >= 0.01)
        out.precision(2);
    else
        out.precision(3);

    out << std::fixed << regularization_weight;
    string rw_string = "_"+out.str();
    return rw_string;
}

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}


//string getCurrentPath(){
//    char szTmp[32];
//    sprintf(szTmp, "/proc/%d/exe", getpid());
//    char pBuf[256];
//    size_t len = sizeof(pBuf);
//    int bytes = MIN(readlink(szTmp, pBuf, len), len - 1);
//    if(bytes >= 0)
//        pBuf[bytes] = '\0';
//    return pBuf;
//}


string whichOS()
{
    #ifdef _WIN32
    return "Windows 32-bit";
    #elif _WIN64
    return "Windows 64-bit";
    #elif __APPLE__ || __MACH__
    return "Mac OSX";
    #elif __linux__
    return "Linux";
    #elif __FreeBSD__
    return "FreeBSD";
    #elif __unix || __unix__
    return "Unix";
    #else
    return "Other";
    #endif
}

//void remove(std::vector<Point> &v)
//{
//    auto end = v.end();
//    for (auto it = v.begin(); it != end; ++it) {
//        end = std::remove(it + 1, end, *it);
//    }

//    v.erase(end, v.end());
//}
