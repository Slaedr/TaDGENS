#ifndef __ACONSTANTS_H

#define __ACONSTANTS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

#define  PI   3.14159265358979323846
#define SQRT2 1.4142135623730950
#define SQRT3 1.73205080756887729353

/// tolerance to check if something is zero, ie, machine epsilon
#define ZERO_TOL 2.2e-16

/// A small number likely smaller than most convergence tolerances
#define SMALL_NUMBER 1e-12

#define NDIM 2
#define NVARS 4

#ifndef MESHDATA_DOUBLE_PRECISION
#define MESHDATA_DOUBLE_PRECISION 20
#endif

namespace acfd
{
	typedef double a_real;
	typedef int a_int;

	/// A type encoding the shape
	enum Shape {LINE, TRIANGLE, QUADRANGLE};
}

#ifndef EIGEN_CORE_H
#include <Eigen/Core>

/// Matrix type to be used for dense linear algebra
typedef Eigen::Matrix<acfd::a_real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

/// Fixed-size ndim x ndim matrix
typedef Eigen::Matrix<acfd::a_real, NDIM, NDIM, Eigen::RowMajor> MatrixDim;

/// Vector type to be used for dense linear algebra
typedef Eigen::Matrix<acfd::a_real, Eigen::Dynamic, 1> Vector;

#endif


#endif
