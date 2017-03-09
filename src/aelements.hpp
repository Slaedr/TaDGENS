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

/// Abstract class for a vector function of a single variable for parameterizing boundary curves
/** The reference nodes are at equally-spaced points between (and including) -1 and 1.
 */
class GeomMapping1D
{
protected:
	int order;
	amat::Array2d<acfd_real> phyNodes;		///< Physical locations of the nodes
	std::vector<acfd_real> speeds;			///< Magnitude of derivative of the curve w.r.t. reference coordinate at quadrature points
	std::vector<Vector> normals;			///< Normals at quadrature points
	std::vector<Vector> mapping;			///< Physical coordinates of the quadrature points, ie the mapping evaluated at the quadrature points
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

	const Vector& normal(const int ipoin) const {
		return normals[ipoin];
	}

	acfd_real speed(const int ipoin) const {
		return speeds[ipoin];
	}

	/// Read-only access to basis function values
	const Vector& map(int ipoin) const {
		return mapping[ipoin];
	}

	/// Computes the curve tangents and normals at a list of points in the reference space
	virtual void computeAll(const amat::Array2d<acfd_real>& points) = 0;
};

/// 1D Lagrange basis parameterization of a boundary
/** The normal to the curve is defined as
 * \f[
 * \hat{n} = \hat{t} \times \hat{k}
 * \f]
 * where \f$ \hat{t} \f$ is the unit tangent vector and \f$ \hat{k} \f$ is the "upward" unit vector normal to the 2D domain.
 * If the nodes of the boundary (in intfac) are oriented properly (which would be the order they would appear in phyNodes as well),
 * the unit normal vector should point from the left element to the right element.
 */
class LagrangeMapping1D : public GeomMapping1D
{
public:
	void computeAll(const amat::Array2d<acfd_real>& points);
};

/// Abstract geometric mapping between a 2D physical element and a reference element
/** For some kinds of elements we may not need Jacobian matrices but only the determinant. computeMapping is for that case.
 */
class GeomMapping2D
{
protected:
	int order;
	amat::Array2d<acfd_real> phyNodes;			///< Physical coordinates of the nodes
	std::vector<Matrix> jaco;					///< Jacobian matrix of the mapping
	std::vector<Matrix> jacoinv;				///< Inverse of the Jacobian matrix
	std::vector<acfd_real> jacodet;				///< Determinant of the Jacobian matrix
	std::vector<Vector> mapping;				///< Physical coords of the quadrature points

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

	/// Sets the basis function values, jacobians, jacobian inverses and jacobian determinants corresponding to a list of reference points
	virtual void computeAll(/*const acfd_real xi[NDIM], Matrix& jac, Matrix& jacinv, acfd_real& jdet*/ const amat::Array2d<acfd_real>& points) = 0;

	/// Computes basis function values at quadrature points
	virtual void computeMappingAndJacobianDet(const amat::Array2d<acfd_real>& points) = 0;

	/// Read-only access to the mapping evaluated at quadrature points
	const Vector& map(int ipoin) const {
		return mapping[ipoin];
	}

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

/// Lagrange geometric mapping on the reference triangle
/** The ref triangle is the one having vertices (0,0), (1,0), (0,1) in that order.
 */
class LagrangeMapping2DTriangle : public GeomMapping2D
{
public:
	LagrangeMapping2DTriangle();
	void computeJacobians(const amat::Array2d<acfd_real>& points);
};

/// Lagrange geometric mapping on the reference square
/** The reference square's vertices are (-1,-1), (1,-1), (1,1), (-1,1) in that order.
 */
class LagrangeMapping2DQuadrangle : public GeomMapping2D
{
public:
	LagrangeMapping2DQuadrangle();
	void computeJacobians(const amat::Array2d<acfd_real>& points);
};

/// Abstract finite element
class Element
{
protected:
	int order;
	std::vector<Array2d<acfd_real>> basis;
	std::vector<Matrix> basisGrad;
	GeomMapping* map;

public:
};

/// Element described by Taylor basis functions and Lagrange geometric mapping
/** The Taylor basis functions are defined on the \emph{physical} element.
 */
class TaylorElement : public Element
{
public:
};

}
#endif
