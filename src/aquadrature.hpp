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
	int ngauss;										///< Number of quadrature points
	int nPoly;										///< Degree of polynomial to integrate exactly
	amat::Array2d<acfd_real> gweights;
	amat::Array2d<acfd_real> gpoints;
public:
	virtual void initialize(const int n_poly) = 0;
	
	const amat::Array2d<acfd_real>& weights() const {
		return gweights;
	}

	const amat::Array2d<acfd_real>& points() const {
		return gpoints;
	}

	int numGauss() const {
		return ngauss;
	}
};

/// 1D Gauss-Legendre quadrature
class Quadrature1D : public QuadratureRule
{
public:
	void initialize(const int n_poly);
};

class Quadrature2D : public QuadratureRule
{
public:
	virtual void initialize(const int n_poly) = 0;
};

/// Integration over the reference square
/** Note that currently, this is restricted to having the same number of quadrature points in the x- and y-directions.
 */
class Quadrature2DSquare : public Quadrature2D
{
public:
	void initialize(const int n_poly);
};

/// Integration over the reference triangle [(0,0), (1,0), (0,1)]
class Quadrature2DTriangle : public Quadrature2D
{
public:
	void initialize(const int n_poly);
};

} // end namespace acfd
#endif
