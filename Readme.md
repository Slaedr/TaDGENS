TADGENS
=======

A Taylor-basis discontinuous Galerkin solver for compressible flow problems. See ./doc/theory for a description of the method used.

Building
--------
A CMake build system is used. GCC g++ is the compiler regularly built against, but clang++ is theoretically supported as well.

For compiling, the variable "EIGEN_DIR" has to be set to the root directory of Eigen 3 library (version 3.3.3 is used in testing). Either set the environment variable, or pass it to make as "EIGEN_DIR=/path/to/eigen-3.3.3" and to cmake as "-DEIGEN\_DIR=/path/to/eigen-3.3.3". 

Running
-------
The executables should be called with the path to a control file as input. Note that in control files, the locations of mesh files and output files should be relative to the directory from which the executable is called.

Roadmap (in order of priority)
------------------------------
- Explicit TVD RK time-stepping scheme for both steady and unsteady inviscid flow problems with smooth solutions
- Reconstruction DG scheme for cost-effective higher-order accuracy
- P-multigrid scheme with explicit smoothers for steady problems
- P-multigrid scheme with implicit smoother at P0 level for steady problems
- Viscous fluxes for compressible viscous flow
- High-order implicit time-stepping in physical time coupled with p-multigrid in pseudo time


Copyright 2017 Aditya Kashi. See LICENSE.md for terms of redistribution with/without modification and those of linking.
