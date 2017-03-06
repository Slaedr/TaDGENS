/** @file aquadrature.hpp
 * @brief Provides Gauss quadrature rules
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#ifndef __AQUADRATURE_H
#define __AQUADRATURE_H

#ifndef __ACONSTANTS_H
#include "aconstants.hpp"
#endif

#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

namespace acfd {

/// Base class for defining quadrature rules
class QuadratureRule
{
protected:
	int ngauss;
	amat::Array2d<acfd_real> gweights;
	amat::Array2d<acfd_real> gpoints;
public:
	QuadratureRule(const int n_gauss) : nguass(n_guass) { }
	
	const amat::Array2d<acfd_real>& weight() const {
		return gweights;
	}

	const amat::Array2d<acfd_real>& point() const {
		return gpoints;
	}
};

/// 1D Gauss-Legendre quadrature
class Quadrature1D : public QuadratureRule
{
public:
	Quadrature1D(const int n_guass);
};

class Quadrature2D : public QuadratureRule
{

};

/// Integration over the reference square
/** Note that currently, this is restricted to having the same number of quadrature points in the x- and y-directions.
 */
class Quadrature2DSquare : public Quadrature2D
{
public:
	Quadrature2DSquare(const int n_guass);
};

/// Integration over the reference triangle [(0,0), (1,0), (0,1)]
class Quadrature2DTriangle : public Quadrature2D
{
public:
	Quadrature2DTriangle(const int n_guass);
};

} // end namespace acfd
#endif
