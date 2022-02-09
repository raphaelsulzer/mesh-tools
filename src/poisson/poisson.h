#ifndef POISSON_H
#define POISSON_H

#include "base/cgal_typedefs.h"
#include <IO/fileIO.h>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// POISSON RECONSTRUCTION ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class PoissonReconstructor{

public:
    dirHolder dir;
    exportOptions eo;
    PoissonReconstructor(dirHolder a, exportOptions b);

void makePoissonMesh(dataHolder& data);

void run(runningOptions options);



}; // end of poisson class

#endif // POISSON_H
