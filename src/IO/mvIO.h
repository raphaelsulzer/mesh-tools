#ifdef OpenMVS
#include <IO/fileIO.h>

// OpenMVS includes
//#include "MVS/Common.h"
//#include "MVS/Scene.h"

// fix with PI from here: https://github.com/cdcseacave/openMVS/issues/643
// otherwise nameclash with some CGAL headers
#pragma push_macro("PI")
#undef PI
#include "MVS.h"
#pragma pop_macro("PI")


//using namespace MVS;

int loadOMVSScene(dirHolder dir, dataHolder& data);

int omvsCleanMesh(dirHolder& dir, dataHolder& data);

#endif
