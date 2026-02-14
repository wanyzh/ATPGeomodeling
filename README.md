# ATPGeomodeling

VTK-based geomodeling tool that generates 3D prism meshes from layered surfaces with fracture constraints, and exports results to VTK formats.

## Project Layout

- `src/` Core library (`CATPGeoModel`)
- `example/` Demo program and example input data
- `CMakeLists.txt` CMake build entry

## Requirements

- CMake (3.15+)
- A C++17 compiler
- VTK (development package) discoverable by CMake (`find_package(VTK CONFIG REQUIRED)`)
- Eigen3 (optional; header-only)

## Build (CMake)

Configure and generate:

```bash
cmake -S . -B build
```

Build:

```bash
cmake --build build --target ATPExtrudeWithConstrain
```

## Run Demo

Run from the output directory so the demo data files can be found:

```bash
cd build/bin
./ATPExtrudeWithConstrain
```

Outputs are written under `build/bin/` (e.g. `geomodelGrid.vtk`, `surfOutput.vtm`, `geomodelFracs.vtm`, `surfRefined.vtm`).

## Notes

- If CMake cannot find VTK, set `VTK_DIR` to the directory containing `VTKConfig.cmake`, or add the VTK prefix to `CMAKE_PREFIX_PATH`.
