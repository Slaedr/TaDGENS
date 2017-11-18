TADGENS
=======

A discontinuous Galerkin solver aimed towards fluid flow problems. One relatively rare feature is the option of using Taylor basis finite elements. See ./doc/theory and Doxygen comments in the code for a description of the methods used.

Building
--------
A CMake build system is used. GCC g++ is the compiler regularly built against, but clang++ is theoretically supported as well.

For compiling, the variable "EIGEN\_DIR" has to be set to the root directory of Eigen 3 library (version 3.3.3 is used in testing). Either set the environment variable, or pass it to make as "EIGEN\_DIR=/path/to/eigen-3.3.3" and to cmake as "-DEIGEN\_DIR=/path/to/eigen-3.3.3". 

Running
-------
The executables should be called with the path to a control file as input. Note that in control files, the locations of mesh files and output files should be relative to the directory from which the executable is called.


Copyright (C) 2017 Aditya Kashi. See LICENSE.md for terms of redistribution with/without modification and those of linking.

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.
