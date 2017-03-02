/** @file elements.hpp
 * @brief Classes for geometric mappings and finite elements.
 * @author Aditya Kashi
 * @date 2017 March 1
 */

#ifndef __AELEMENTS_H
#define __AELEMENTS_H

#ifndef __ACONSTANTS_H
#include "aconstants.hpp"
#endif

#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

#ifndef EIGEN_CORE_H
#include <Eigen/Core>
typedef Eigen::Matrix<acfd::acfd_real, Eigen::dynamic, Eigen::dynamic> Matrix;
#endif

namespace acfd {

/// Abstract geometric mapping between a physical element and a reference element
class GeomMapping
{
protected:
	int order;
	Array2d<acfd_real> phyNodes;
	std::vector<Matrix> jaco;
	std::vector<Matrix> jacoinv;
	std::vector<acfd_real> jacodet;

public:
	/// Set polynomial order of mapping
	void setOrder(const int ord) {
		order = ord;
	}

	/// Return the order
	int getOrder() const {
		return order;
	}

	/// Sets the coordinates of physical nodes of the element
	void setPhysicalNodes(const Array2d& physicalnodes) {
		phyNodes = physicalnodes;
	}

	/// Sets quadrature points of the reference element
	/*void setQuadraturePoints(const Array2d& quadraturePoints, const Array2d& quadratureWeights) {
		quadpoints = quadraturePoints;
		quadWeights = quadratureWeights;
	}*/

	/// Sets the jacobian, jacobian inverse and jacobian determinant corresponding to some reference coordinates
	virtual void computeJacobians(const acfd_real xi[NDIM], Matrix& jac, Matrix& jacinv, acfd_real& jdet) = 0;

	/// Read-only access to jacobians
	virtual const Matrix& jac(int ipoint) const = 0;
	
	/// Read-only access to inverse of jacobians
	virtual const Matrix& jacInv(const int ipoint) const = 0;

	/// Jacobian determinant
	virtual acfd_real jacDet(const int ipoint) const = 0;
};

class LagrangeMapping2DTriangle : public GeomMapping
{
public:
};

class LagrangeMapping2DQuadrangle : public GeomMapping
{
public:
};

class LagrangeMapping1D : public GeomMapping
{
public:
};

/// Abstract finite element
class Element
{
protected:
	int order;
	std::vector<Array2d<acfd_real>> basis;
	std::vector<Array2d<acfd_real>> basisGrad;
	GeomMapping* map;

public:
};

/// Element described by Taylor basis functions and Lagrange geometric mapping
class TaylorElement : public Element
{
public:
};

}
#endif
