/** @file elements.hpp
 * @brief Classes for geometric mappings and finite elements.
 *
 * Note that the geometric mapping classes do not let one compute the Jacobian etc. at any point;
 * in other words, they do not provide functions into which reference/physical coordinates can be plugged in.
 * Rather, they compute and store values of these functions at the quadrature points.
 *
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

#ifndef __AQUADRATURE_H
#include "aquadrature.hpp"
#endif

namespace acfd {

/// A type encoding the shape of 2D elements
enum Shape2d {TRIANGLE, QUADRANGLE};

/// Abstract class for a vector function of a single variable for parameterizing boundary curves
/** The reference nodes are at equally-spaced points between (and including) -1 and 1.
 */
class GeomMapping1D
{
protected:
	int degree;								///< Polynomial degree of the mapping
	amat::Array2d<acfd_real> phyNodes;		///< Physical locations of the nodes
	std::vector<Vector> normals;			///< Normals at quadrature points. Note that these are NOT unit vectors; they're the 1D analogue of "area" vectors
	amat::Array2d<acfd_real> mapping;		///< Physical coordinates of the quadrature points, ie the mapping evaluated at the quadrature points
	const Quadrature1D* quadrature;			///< Gauss points and weights for integrating quantities

public:
	/// Return the order
	int getDegree() const {
		return degree;
	}

	/// Sets the polynomial degree and coordinates of physical nodes of the element and the quadrature context
	void setAll(const int deg, const amat::Array2d<acfd_real>& physicalnodes, const Quadrature1D* q) {
		degree = deg;
		phyNodes = physicalnodes;
		quadrature = q;
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
	const amat::Array2d<acfd_real>& map() const {
		return mapping;
	}

	/// Computes the curve normals at a list of points in the reference space
	/** \param[in] points is an npoin x ndim array, ie, whose each row contains the coordinates of one of the quadrature points
	 */
	virtual void computeAll(const amat::Array2d<acfd_real>& points) = 0;

	/// Access to quadrature context
	const Quadrature1D* getQuadrature() const {
		return quadrature;
	}
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
	Shape2d shape;								///< Shape of the element
	int degree;									///< Polynomial degree of the map
	amat::Array2d<acfd_real> phyNodes;			///< Physical coordinates of the nodes
	std::vector<Matrix> jaco;					///< Jacobian matrix of the mapping
	std::vector<Matrix> jacoinv;				///< Inverse of the Jacobian matrix
	std::vector<acfd_real> jacodet;				///< Determinant of the Jacobian matrix
	amat::Array2d<acfd_real> mapping;			///< Physical coords of the quadrature points
	const Quadrature2D* quadrature;				///< Gauss points and weights for integrating quantities

public:
	/// Return the order
	int getDegree() const {
		return degree;
	}

	Shape2d getShape() const {
		return shape;
	}

	/// Sets the polynomial degree, coordinates of physical nodes of the element and the integration context
	void setAll(const int deg, const amat::Array2d<acfd_real>& physicalnodes, const Quadrature2D* const quad) {
		degree = deg;
		phyNodes = physicalnodes;
		quadrature = quad;
	}

	/// Sets the basis function values, jacobians, jacobian inverses and jacobian determinants corresponding to the quadrature points
	/** This function also allocates storage for all member data except phyNodes, which must be set beforehand.
	 */
	virtual void computeAll() = 0;

	/// Computes basis function values at quadrature points
	/** Note that storage is allocated only for mapping and jacodet.
	 */
	virtual void computeMappingAndJacobianDet() = 0;

	/// Read-only access to physical node locations
	const amat::Array2d<acfd_real>& getPhyNodes() const {
		return phyNodes;
	}

	/// Read-only access to the mapping evaluated at quadrature points
	const amat::Array2d<acfd_real>& map() const {
		return mapping;
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

	/// Access to quadrature context
	const Quadrature2D* getQuadrature() const {
		return quadrature;
	}
};

/// Lagrange geometric mapping on the reference triangle
/** The ref triangle is the one having vertices (0,0), (1,0), (0,1) in that order.
 */
class LagrangeMapping2DTriangle : public GeomMapping2D
{
public:
	void computeAll();
	void computeMappingAndJacobianDet();
};

/// Lagrange geometric mapping on the reference square
/** The reference square's vertices are (-1,-1), (1,-1), (1,1), (-1,1) in that order.
 */
class LagrangeMapping2DQuadrangle : public GeomMapping2D
{
public:
	void computeAll();
	void computeMappingAndJacobianDet();
};

/// Abstract finite element
class Element
{
protected:
	int degree;										///< Polynomial degree
	int ndof;										///< Number of local DOFs
	std::vector<amat::Array2d<acfd_real>> basis;	///< Values of basis functions at quadrature points
	std::vector<Matrix> basisGrad;					///< Values of derivatives of the basis functions at the quadrature points
	const GeomMapping2D* gmap;						///< The 2D geometric map which maps this element to the reference element

public:
	/// Set the data, compute geom map, and compute basis and basis grad
	/** \param[in] geommap The geometric mapping should be initialized and all values computed beforehand; we'll not do that here
	 */
	virtual void initialize(int degr, const GeomMapping2D* geommap) = 0;

	/// Read-only access to basis at a given quadrature point
	const amat::Array2d<acfd_real>& bFunc(const int ipoin) {
		return basis[ipoin];
	}

	/// Read-only access to basis gradients at a given quadrature point
	const Matrix& bGrad(const int ipoin) {
		return basisGrad[ipoin];
	}

	int getDegree() const {
		return degree;
	}

	int getNumDOFs() const {
		return ndof;
	}
};

/// Abstract element defined on the physical element
/** Taylor basis element is an implementation of this, while Lagrange element is generally not.
 */
class Element_PhysicalSpace : public Element
{
public:
	virtual void initialize(int degr, const GeomMapping2D* geommap) = 0;

	/// Computes values of basis functions at a given point in physical space
	virtual void computeBasis(const Vector& point, std::vector<acfd_real>& basisvalues) = 0;
};

/// Element described by Taylor basis functions
/** The Taylor basis functions are defined on the \emph{physical} element.
 *
 * Note that the first (p0) DOF is not the value at the element center, but the average value over the element.
 * We thus need to compute offsets from Taylor polynomials for terms associated with P2 and higher.
 */
class TaylorElement : public Element_PhysicalSpace
{
	acfd_real area;										///< Area of the element
	acfd_real center[NDIM];								///< Physical location of element's geometric center
	acfd_real delta[NDIM];								///< Maximum extent of the element in the coordinate directions
	std::vector<std::vector<acfd_real>> basisOffset;	///< The quantities by which the basis functions are offset from actual Taylor polynomial basis
public:
	void initialize(int degr, const GeomMapping2D* geommap);
	void computeBasis(const acfd_real* point, acfd_real* basisvalues);
};

/// An interface "element" between 2 adjacent finite elements with basis defined on physical elements
/** 
 * In future, perhaps intfac data could be stored in this class.
 */
class FaceElement_PhysicalSpace
{
	const Element_PhysicalSpace* leftel;				///< "Left" element
	const Element_PhysicalSpace* rightel;				///< "Right" element
	amat::Array2d<acfd_real> leftbasis;					///< Values of the left element's basis functions at the face quadrature points
	amat::Array2d<acfd_real> rightbasis;				///< Values of the left element's basis functions at the face quadrature points
	const GeomMapping1D* gmap;							///< 1D geometric mapping (parameterization) of the face

public:
	/// Sets data; computes basis function values of left and right element at each quadrature point
	/** NOTE: Call only after element data has been precomputed, ie, by calling the compute function on the elements, first!
	 * \param[in] geommap The geometric mapping must be initialized externally; we don't do it here
	 */
	void initialize(int degr, const Element_PhysicalSpace* lelem, const Element_PhysicalSpace* relem, const GeomMapping1D* geommap);

	/// Read-only access to basis function values from left element
	const amat::Array2d<acfd_real>& leftBasis() {
		return leftbasis;
	}

	/// Read-only access to basis function values from right element
	const amat::Array2d<acfd_real>& rightBasis() {
		return rightbasis;
	}
};

}
#endif
