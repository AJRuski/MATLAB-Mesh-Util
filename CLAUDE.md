# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

MATLAB-Mesh-Util is a library of MATLAB functions for working with 3D models represented as faces and vertices (face/vertex meshes).

There is no build or lint step. Unit tests are in `tests/MeshUtilTests.m` (a `matlab.unittest` class covering every function). Run them from the repo root in MATLAB:

```matlab
results = runtests("tests");                      % all tests
results = runtests("tests/MeshUtilTests.m", "Name", "MeshUtilTests/SplitBoxAlongAxis");   % one test
```

The tests were written without being able to run MATLAB, so a failure may be a bug in the test as well as in the library. New functions should get tests in the same file, using the shape builders at the bottom of the class (`MakeBox`, `MakeCylinder`, `MakeSphere`, `MakeLPrism`) and `VerifyClosedPiece`.

## Environment

- Target MATLAB version: R2024b or newer.
- Deployment has all toolboxes available, so toolbox functions may be used freely.
- MATLAB is **not installed** in the development environment, so code cannot be run or tested here. Write carefully, and do not claim code has been tested or verified by execution.

## Code conventions

- All code is human-reviewed: favor readable, easy-to-follow code over clever or compact code.
- Include comments explaining what the code does and why.
- Every function must begin with an `arguments` block that validates its inputs.
- Function names and variable names use PascalCase (e.g. `ComputeFaceNormals`, `VertexCount`).

## Mesh conventions

- Meshes are passed as `Faces` (Nx3 triangle vertex indices, 1-based) and `Verts` (Mx3 positions). Meshes are never mixed tri/quad.
- Only `QuadsToTris` accepts Nx4 quads; every other function takes triangles, so quads must be converted first.
- Meshes passed to the functions are assumed watertight. `FixFaceNormals` errors if not; use `IsMeshWatertight` to find the problem edges.
- Correct normals means counter-clockwise winding seen from outside (outward). `ComputeMeshVolume` relies on this, and a negative volume means inverted normals.
- Functions that conceptually transform a mesh return `Verts` as an output so the pair can be passed along together. It is unchanged unless the function edits vertices (`FixDegenerateTriangles` moves and renumbers them).
