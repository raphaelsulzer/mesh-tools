# mesh-tools
Tools for reconstructing, processing and evaluating a surface mesh from a point cloud.



## Common interface for all tools

All tools have a couple of common command line parameters. To get a more complete and detailed list of all possible parameters use
the help command.

| Command               | Description                                                                             |
|:----------------------|:----------------------------------------------------------------------------------------|
| `-h`                  | Display a help message, explaining all command line options of the specific tool.       |
| `-w pathToWorkingDir` | All input and output files are located inside the specified directory.                  |
| `-i inputFileName`    | The file name of the input for the tool. The file ending does not have to be specified. |
| `-o outputFileName`   | The file name of the output file(s). If left empty, the input file name is used.        |
| `-e exportOptions`    | Specify which files and which properties should be exported.                            |



## Available tools

### scan

Synthetically scan a polygonal mesh, e.g. for [Deep Surface Reconstruction from Point Clouds with Visibility Information](https://github.com/raphaelsulzer/dsrv-data).

### omvs2npz

Extract points, normals and sensors from an [OpenMVS](https://github.com/cdcseacave/openMVS) project file and save them as arrays in an `.npz` file.


### labatut

Reconstruct a mesh from a point cloud using the algorithm presented in [Labatut et al. 2009](https://diglib.eg.org/handle/10.2312/CGF.v28i8pp2275-2290).

### feat

Extract features from a visibility augmented 3D Delaunay triangulation for [DGNN](https://github.com/raphaelsulzer/dgnn).
If a ground truth file is specified, ground truth occupancies for each cell of the 3D Delaunay triangulation will also be exported.


### occ2mesh

Transform ground truth cell occupancies as computed by `feat` to a surface mesh.

### normal

Estimate normals for a mesh.

### sample

Sample a mesh.

### collapse

Apply the edge collapse algorithm.


## Installation

```bash
git clone git@github.com:raphaelsulzer/mesh-tools.git
cd mesh-tools
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make -j
```
