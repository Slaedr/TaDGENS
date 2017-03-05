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

namespace acfd {

/// Abstract geometric mapping between a 2D physical element and a reference element
class GeomMapping2D
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

	/// Returns the value of the mapping function at given reference coordinates
	virtual void evaluateMapping(const acfd_real xi[NDIM], acfd_real xg[NDIM]) = 0;
	
	/// Sets the jacobian, jacobian inverse and jacobian determinant corresponding to some reference coordinates
	virtual void computeJacobians(const acfd_real xi[NDIM], Matrix& jac, Matrix& jacinv, acfd_real& jdet) = 0;

	/// Read-only access to jacobians
	const Matrix& jac(int ipoint) const {
		return jaco[ipoint];
	}
	
	/// Read-only access to inverse of jacobians
	const Matrix& jacInv(const int ipoint) const {
		return jacoinv[ipoint];
	}

	/// Jacobian determinant
	acfd_real jacDet(const int ipoint) const {
		return jacodet[ipoint];
	}
};

class LagrangeMapping2DTriangle : public GeomMapping2D
{
public:
};

class LagrangeMapping2DQuadrangle : public GeomMapping2D
{
public:
};

/// Abstract class for a vector function of a single variable for parameterizing boundary curves
/** The reference nodes are at equally-spaced points between (and including) -1 and 1.
 */
class GeomMapping1D
{
protected:
	int order;
	amat::Array2d<acfd_real> phyNodes;		///< Physical locations of the nodes
	std::vector<acfd_real> speeds;			///< Magnitude of derivative of the curve w.r.t. reference coordinate at several points
	std::vector<Vector> normals;			///< Normals at several points
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
	
	/// Returns the value of the mapping function at a given reference coordinate
	virtual void evaluateMapping(const acfd_real xi, acfd_real xg[NDIM]) = 0;

	const Vector& normal(const int ipoin) const {
		return normals[ipoin];
	}

	acfd_real speed(const int ipoin) const {
		return speeds[ipoin];
	}
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
