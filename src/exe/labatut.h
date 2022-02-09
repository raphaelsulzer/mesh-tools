#ifndef SURFACERECON_H
#define SURFACERECON_H

#include <IO/fileIO.h>
#include <exe/exeOptions.h>

//void surfaceMerging(string path, string filename1, string filename2,
//                    string scoreType, string regularization_term, double regularization_weight);
//void surfaceReconstruction(string path, string filename,
//                           int sample,
//                           int number_of_images,
//                           string scoreType, string regularization_term, double regularization_weight);
//void poissonReconstruction(string path, string filename1, string filename2);

int surfaceReconstruction(dirHolder& dir, dataHolder& data, runningOptions& options, exportOptions& exportO);

#endif // SURFACERECON_H
