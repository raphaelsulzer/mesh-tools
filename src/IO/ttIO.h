#ifndef TTIO_H
#define TTIO_H

#include <base/cgal_typedefs.h>
#include <exe/exeOptions.h>


//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////

int importTanksAndTemplesIs(dirHolder dir, dataHolder& data, Point sensor_position);
int importTanksAndTemples(dirHolder dir, dataHolder& data, runningOptions options);
int importTanksAndTemplesScannerLocations(dirHolder dir, dataHolder& data, runningOptions options);

#endif

