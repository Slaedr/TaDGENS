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
	int degee;								///< Polynomial degree of the mapping
	amat::Array2d<acfd_real> phyNodes;		///< Physical locations of the nodes
	std::vector<Vector> normals;			///< Normals at quadrature points. Note that these are NOT unit vectors; they're the 1D analogue of "area" vectors
	std::vector<Vector> mapping;			///< Physical coordinates of the quadrature points, ie the mapping evaluated at the quadrature points

public:
	/// Return the order
	int getDegree() const {
		return degree;
	}

	/// Sets the polynomial degree and coordinates of physical nodes of the element
	void setAll(const int deg, const Array2d& physicalnodes) {
		degree = deg;
		phyNodes = physicalnodes;
	}

	/// Read-only access to physical node locations
	const amat::Array2d<acfd_real>& getPhyNodes() const {
		return phyNodes;
	}

	/// Read-only access to the "area" vectors normal to the face at quadrature points
	const Vector& normal(const int ipoin) const {
		return normals[ipoin];
	}

	/// Read-only access to physical coords of the quadrature points
	const Vector& map(int ipoin) const {
		return mapping[ipoin];
	}

	/// Computes the curve normals at a list of points in the reference space
	/** \param[in] points is an npoin x ndim array, ie, whose each row contains the coordinates of one of the quadrature points
	 */
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
	int degree;									///< Polynomial degree of the map
	amat::Array2d<acfd_real> phyNodes;			///< Physical coordinates of the nodes
	std::vector<Matrix> jaco;					///< Jacobian matrix of the mapping
	std::vector<Matrix> jacoinv;				///< Inverse of the Jacobian matrix
	std::vector<acfd_real> jacodet;				///< Determinant of the Jacobian matrix
	std::vector<Vector> mapping;				///< Physical coords of the quadrature points

public:
	/// Return the order
	int getDegree() const {
		return degree;
	}

	/// Sets the polynomial degree and coordinates of physical nodes of the element
	void setAll(const int deg, const Array2d<acfd_real>& physicalnodes) {
		degree = deg;
		phyNodes = physicalnodes;
	}

	/// Sets the basis function values, jacobians, jacobian inverses and jacobian determinants corresponding to a list of reference points
	/** This function also allocates storage for all member data.
	 * \param[in] points is an npoin x ndim array, ie, whose each row contains the coordinates of one of the quadrature points
	 */
	virtual void computeAll(/*const acfd_real xi[NDIM], Matrix& jac, Matrix& jacinv, acfd_real& jdet*/ const amat::Array2d<acfd_real>& points) = 0;

	/// Computes basis function values at quadrature points
	/** Note that storage is allocated only for mapping and jacodet.
	 * \param[in] points is an npoin x ndim array, ie, whose each row contains the coordinates of one of the quadrature points
	 */
	virtual void computeMappingAndJacobianDet(const amat::Array2d<acfd_real>& points) = 0;

	/// Read-only access to physical node locations
	const amat::Array2d<acfd_real>& getPhyNodes() const {
		return phyNodes;
	}

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
	void computeAll(const amat::Array2d<acfd_real>& points);
	void computeMappingAndJacobianDet(const amat::Array2d<acfd_real>& points);
};

/// Lagrange geometric mapping on the reference square
/** The reference square's vertices are (-1,-1), (1,-1), (1,1), (-1,1) in that order.
 */
class LagrangeMapping2DQuadrangle : public GeomMapping2D
{
public:
	void computeAll(const amat::Array2d<acfd_real>& points);
	void computeMappingAndJacobianDet(const amat::Array2d<acfd_real>& points);
};

/// Abstract finite element
class Element
{
protected:
	int degree;									///< Polynomial degree
	int ngauss;									///< Number of quadrature points
	std::vector<Array2d<acfd_real>> basis;		///< Values of basis functions at quadrature points
	std::vector<Matrix> basisGrad;				///< Values of derivatives of the basis functions at the quadrature points
	const GeomMapping2D* gmap;					///< The 2D geometric map which maps this element to the reference element
	const Quadrature2D* quad;					///< Numerical integration context

public:
	/// Set the data, compute geom map, and compute basis and basis grad
	/** \param[in] geommap The geometric mapping should be initialized and all values computed externally; we'll not do that here
	 */
	virtual void initialize(int degr, int nquadpoin, const Quadrature2D* q, const GeomMapping2D* geommap) = 0;

	/// Read-only access to basis at a given quadrature point
	const Array2d<acfd_real>& basisFunctions(const int ipoin) {
		return basis[ipoin];
	}

	/// Read-only access to basis gradients at a given quadrature point
	const Matrix& basisGradients(const int ipoin) {
		return basisGrad[ipoin];
	}
};

/// Element described by Taylor basis functions and Lagrange geometric mapping
/** The Taylor basis functions are defined on the \emph{physical} element.
 */
class TaylorElement : public Element
{
public:
	void initialize(int degr, int nquadpoin, const Quadrature2D* q, const Array2d<acfd_real>& phynodes);
};

/// An interface "element" between 2 adjacent finite elements
/** 
 * In future, perhaps intfac data could be stored in this class.
 */
class FaceElement
{
	const Element* leftel;							///< "Left" element
	const Element* rightel;							///< "Right" element
	std::vector<Array2d<acfd_real>> leftbasis;		///< Values of the left element's basis functions at the face quadrature points
	std::vector<Array2d<acfd_real>> rightbasis;		///< Values of the left element's basis functions at the face quadrature points
	const GeomMapping1D* gmap;						///< 1D geometric mapping (parameterization) of the face
	const Quadrature1D* quad;						///< Numerical integration context

public:
	/// Sets data; computes basis function values of left and right element at each quadrature point
	/** NOTE: Call only after element data has been precomputed, ie, by calling the compute function on the elements, first!
	 * \param[in] geommap The geometric mapping must be initialized externally; we don't do it here
	 */
	void initialize(int degr, int nquadpoin, const Element* lelem, const Element* relem, const Quadrature1D* q, const GeomMapping1D* gmap);

	/// Read-only access to basis function values from left element
	const Array2d<acfd_real>& leftBasis(const int ipoin) {
		return leftbasis[ipoin];
	}

	/// Read-only access to basis function values from right element
	const Array2d<acfd_real>& rightBasis(const int ipoin) {
		return rightbasis[ipoin];
	}
};

}
#endif
