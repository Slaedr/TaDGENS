/** @file elements.hpp
 * @brief Classes for geometric mappings and finite elements.
 *
 * Note that the geometric mapping classes do not let one compute the Jacobian etc. at any point;
 * in other words, they do not provide functions into which reference/physical coordinates can be plugged in.
 * Rather, they compute and store values of these functions at the quadrature points.
 *
 * TODO: Verify P2 geometric map for triangles, P1 and P2 map for quads
 * TODO: Verify P1 Taylor basis on P2 triangles, P2 Taylor basis for triangles, P1 and P2 for quads
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

/// Abstract class for a vector function of a single variable for parameterizing boundary curves
/** The reference nodes are at equally-spaced points between (and including) -1 and 1.
 */
class GeomMapping1D
{
protected:
	int degree;								///< Polynomial degree of the mapping
	Matrix phyNodes;						///< Physical locations of the nodes
	std::vector<Vector> normals;			///< Unit normals at quadrature points
	std::vector<a_real> speeds;				///< Magnitude of tangent vectors to the curve at quadrature points
	Matrix mapping;							///< Physical coordinates of the quadrature points, ie the mapping evaluated at the quadrature points
	const Quadrature1D* quadrature;			///< Gauss points and weights for integrating quantities

public:
	/// Return the order
	int getDegree() const {
		return degree;
	}

	/// Sets the polynomial degree and coordinates of physical nodes of the element and the quadrature context
	void setAll(const int deg, const Matrix& physicalnodes, const Quadrature1D* q) {
		degree = deg;
		phyNodes = physicalnodes;
		quadrature = q;
	}

	/// Read-only access to physical node locations
	const Matrix& getPhyNodes() const {
		return phyNodes;
	}

	/// Read-only access to the curve "speed"
	const std::vector<a_real>& speed() const {
		return speeds;
	}

	/// Read-only access to unit vectors normal to the face at quadrature points
	const std::vector<Vector>& normal() const {
		return normals;
	}

	/// Read-only access to physical coords of the quadrature points
	const Matrix& map() const {
		return mapping;
	}

	/// Computes the curve normals at quadrature points in the reference space
	virtual void computeAll() = 0;

	/// Access to quadrature context
	const Quadrature1D* getQuadrature() const {
		return quadrature;
	}

	virtual ~GeomMapping1D() { }
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
	void computeAll();
};

/// Abstract geometric mapping between a 2D physical element and a reference element
class GeomMapping2D
{
protected:
	Shape shape;								///< Shape of the element
	int degree;									///< Polynomial degree of the map
	Matrix phyNodes;							///< Physical coordinates of the nodes
	std::vector<MatrixDim> jaco;				///< Jacobian matrix of the mapping
	std::vector<MatrixDim> jacoinv;				///< Inverse of the Jacobian matrix
	std::vector<a_real> jacodet;				///< Determinant of the Jacobian matrix
	Matrix mapping;								///< Physical coords of the quadrature points
	const Quadrature2D* quadrature;				///< Gauss points and weights for integrating quantities

public:
	/// Return the order
	int getDegree() const {
		return degree;
	}

	Shape getShape() const {
		return shape;
	}

	/// Sets the polynomial degree, coordinates of physical nodes of the element and the integration context
	void setAll(const int deg, const Matrix& physicalnodes, const Quadrature2D* const quad) {
		degree = deg;
		phyNodes = physicalnodes;
		quadrature = quad;
	}

	/// Allocates and sets the Jacobian inverses and Jacobian determinants corresponding to the quadrature points
	/** Call this function only after [setting up](@ref setAll).
	 */
	virtual void computeForReferenceElement() = 0;

	/// Computes basis function values and Jacobian determinants at quadrature points
	/** Note that storage is allocated only for mapping and jacodet, not for Jacobian matrix or its inverse.
	 */
	virtual void computeForPhysicalElement() = 0;

	/// Computes the Jacobian inverse and its determinant at a set of reference coordinates
	virtual void calculateJacobianDetAndInverse(const Matrix& refpoints, std::vector<MatrixDim>& jacoi, std::vector<a_real>& jacod) const = 0;

	/// Calculate physical coords of domain quadrature points separately, in case needed by referential elements
	virtual void computePhysicalCoordsOfDomainQuadraturePoints() = 0;
	
	/// Computes physical locations of points given their reference coordinates
	virtual void calculateMap(const Matrix& points, Matrix& maps) const = 0;

	/// Read-only access to physical node locations
	const Matrix& getPhyNodes() const {
		return phyNodes;
	}

	/// Read-only access to the mapping evaluated at quadrature points \sa mapping
	const Matrix& map() const {
		return mapping;
	}

	/// Read-only access to jacobians at domain quadrature points
	const std::vector<MatrixDim>& jac() const {
		return jaco;
	}

	/// Read-only access to inverse of jacobians at domain quadrature points
	const std::vector<MatrixDim>& jacInv() const {
		return jacoinv;
	}

	/// Jacobian determinant at domain quadrature points
	const std::vector<a_real>& jacDet() const {
		return jacodet;
	}

	/// Access to quadrature context
	const Quadrature2D* getQuadrature() const {
		return quadrature;
	}

	virtual ~GeomMapping2D() { }
};

/// Lagrange geometric mapping on the reference triangle
/** The ref triangle is the one having vertices (0,0), (1,0), (0,1) in that order.
 * The reference square's vertices are (-1,-1), (1,-1), (1,1), (-1,1) in that order.
 *
 * \note The Jacobian matrix is not actually stored, instead the inverse is stored.
 */
class LagrangeMapping2D: public GeomMapping2D
{
public:
	void computeForReferenceElement();

	void computeForPhysicalElement();
	
	void calculateJacobianDetAndInverse(const Matrix& po, std::vector<MatrixDim>& jacoi, std::vector<a_real>& jacod) const;
	
	void computePhysicalCoordsOfDomainQuadraturePoints();
	
	void calculateMap(const Matrix& __restrict__ points, Matrix& __restrict__ maps) const;
};

/** \brief A type defining whether basis functions are defined in reference space or physical space
 * or whether it's a dummy.
 *
 * For example, Lagrange basis is generally defined as a function of reference coordinates
 * while Taylor basis is defined as a function of physical coordinates.
 */
enum BasisType {REFERENTIAL, PHYSICAL, NONEXISTENT};

/// Abstract finite element
class Element
{
protected:
	BasisType type;									///< Where are the basis functions defined? \sa BasisType
	int degree;										///< Polynomial degree
	int ndof;										///< Number of local DOFs
	Matrix basis;									///< Values of basis functions at quadrature points
	std::vector<Matrix> basisGrad;					///< Values of derivatives of the basis functions at the quadrature points
	const GeomMapping2D* gmap;						///< The 2D geometric map which maps this element to the reference element

public:
	/// Set the data, compute geom map, and compute basis and basis grad
	/** \param[in] geommap The geometric mapping should be initialized beforehand;
	 * however, the computation of required geometric quantities such as the Jacobian is done here.
	 */
	virtual void initialize(int degr, GeomMapping2D* geommap) = 0;

	/// Computes values of basis functions at given points in either reference space or physical space
	virtual void computeBasis(const Matrix& points, Matrix& basisvalues) const = 0;

	/// Computes basis functions' gradients at given points in either reference space or physical space
	virtual void computeBasisGrads(const Matrix& points, std::vector<Matrix>& basisgrads) const = 0;

	virtual ~Element() { }

	/// Computes interpolated values at the quadrature point with index ig from given DOF values
	a_real interpolate(const int ig, const Vector& dofs) const
	{
		a_real val = 0;
		for(int i = 0; i < ndof; i++)
			val += dofs[i]*basis(ig,i);
		return val;
	}

	/// Computes values of a function at domain quadrature points using DOFs supplied
	/** \param[in] dofs The rowmajor matrix of local DOFs (1 row per physical variable)
	 * \param[in|out] values Values of the function at each quadrature point. Each row contains all DOFs of a physical variable
	 * CANNOT be the same as dofs, the first argument.
	 * 
	 * \todo Notice that when we have a system of PDEs, we might want to use a BLAS 3 call to compute all values.
	 */
	void interpolateAll(const Matrix& __restrict__ dofs, Matrix& __restrict__ values) const
	{
		values.noalias() = basis * dofs.transpose();
	}
	
	/// Computes values of the specified component at domain quadrature points using DOFs supplied
	/** \param[in] comp specifies the row to use in the matrix of DOFs
	 */
	void interpolateComponent(const int comp, const Matrix& __restrict__ dofs, Vector& __restrict__ values) const
	{
		values.noalias() = basis * dofs.row(comp).transpose();
	}

	/// Read-only access to basis at a given quadrature point
	const Matrix& bFunc() const {
		return basis;
	}

	/// Read-only access to basis gradients at the element's domain quadrature point
	const std::vector<Matrix>& bGrad() const {
		return basisGrad;
	}

	int getDegree() const {
		return degree;
	}

	int getNumDOFs() const {
		return ndof;
	}

	BasisType getType() const {
		return type;
	}

	const GeomMapping2D* getGeometricMapping() const {
		return gmap;
	}
};

/// Element described by Taylor basis functions
/** The Taylor basis functions are defined on the \emph{physical} element.
 *
 * Note that the first (p0) DOF is not the value at the element center, but the average value over the element.
 * We thus need to compute offsets from Taylor polynomials for terms associated with P2 and higher.
 *
 * Computation of basis gradients does not need the Jacobian of the geometric mapping. They are simply computed as
 * \f[ \nabla B(x) = \nabla B(F(\xi)) \f]
 */
class TaylorElement : public Element
{
	a_real area;										///< Area of the element
	a_real center[NDIM];								///< Physical location of element's geometric center
	a_real delta[NDIM];									///< Maximum extent of the element in the coordinate directions
	std::vector<std::vector<a_real>> basisOffset;		///< The quantities by which the basis functions are offset from actual Taylor polynomial basis
public:
	TaylorElement() {
		type = PHYSICAL;
	}

	/// Sets data, computes geometric map data and computes basis functions and their gradients
	void initialize(int degr, GeomMapping2D* geommap);
	
	/// Computes values of basis functions at a given point in physical space
	void computeBasis(const Matrix& points, Matrix& basisvalues) const;
	
	/// Computes basis functions' gradients at given points in physical space
	void computeBasisGrads(const Matrix& points, std::vector<Matrix>& basisgrads) const;

	void printDetails() const {
		std::printf("  (%f,%f), %f, %f, %f\n", center[0], center[1], delta[0], delta[1], area);
	}

	a_real getArea() const {
		return area;
	}

	const a_real* getCenter() const {
		return center;
	}

	const a_real* getDelta() const {
		return delta;
	}

	const std::vector<std::vector<a_real>>& getBasisOffsets() const {
		return basisOffset;
	}
};

/// Lagrange finite element with equi-spaced nodes
/** Computation of basis function gradients requires geometric Jacobian.
 * \f[ 
 * \nabla B(x) = \nabla \hat{B}(F^{-1}(x)) = J^{-T} \nabla_\xi \hat{B}(F^{-1}(F(\xi)))
 * = \nabla_\xi \hat{B}(\xi)
 * \f]
 */
class LagrangeElement : public Element
{
public:
	LagrangeElement() {
		type = REFERENTIAL;
	}

	/// Sets data and computes basis functions and their gradients
	void initialize(int degr, GeomMapping2D* geommap);
	
	/// Computes values of basis functions at a given point in reference space
	void computeBasis(const Matrix& points, Matrix& basisvalues) const;
	
	/// Computes basis functions' gradients at given points in reference space
	void computeBasisGrads(const Matrix& points, std::vector<Matrix>& basisgrads) const;
	
	/// Returns the locations of nodes in reference space
	Matrix getReferenceNodes() const;
};

/// Just that - a dummy element
/** Used for `ghost' elements on boundary faces.
 */
class DummyElement : public Element
{
public:
	void initialize(int degr, GeomMapping2D* geommap) { type = NONEXISTENT; }
	void computeBasis(const Matrix& points, Matrix& basisvalues) const { };
	void computeBasisGrads(const Matrix& points, std::vector<Matrix>& basisgrads) const { };
};

/// An interface "element" between 2 adjacent finite elements
/** In future, perhaps intfac data could be stored in this class.
 * \todo TODO: Make interpolation functions more efficient for nodal basis functions.
 */
class FaceElement
{
	const Element* leftel;								///< "Left" element
	const Element* rightel;								///< "Right" element
	int llfn, rlfn;										///< Local face number of this face w.r.t the left and right elements
	Matrix leftbasis;									///< Values of the left element's basis functions at the face quadrature points
	Matrix rightbasis;									///< Values of the left element's basis functions at the face quadrature points
	std::vector<Matrix> leftbgrad;						///< left element's basis gradients at face quadrature points
	std::vector<Matrix> rightbgrad;						///< right element's basis gradients at face quadrature points
	const GeomMapping1D* gmap;							///< 1D geometric mapping (parameterization) of the face

	/// Computes 2D reference coordinates on the face of an element that shares this face corresponding to face reference points
	/** \param[in] facepoints 1D coordinate on the face
	 * \param[in] elem Pointer to element
	 * \param[in] lfn Local face number of this face in element elem
	 * \param[in] isright A flag that's 1 if elem is the left element and -1 if it's the right.
	 * \praram[in|out] lpoints Contains 2D reference coordinates of face quadrature points in element elem on output.
	 * \return Returns the magnitude of the speed of the face of the  reference element w.r.t. the reference face.
	 */
	a_real getElementRefCoords(const Matrix& facepoints, const Element *const elem,
		const int lfn, const int isright, Matrix& lpoints);
public:
	/// Sets data; computes basis function values of left and right element at each quadrature point
	/** \note Call only after element data has been precomputed, ie, by calling the compute function on the elements, first!
	 * \param[in] geommap The geometric mapping must be [initialized](@ref GeomMapping1D::setAll) and map and normals [computed externally](@ref GeomMapping1D::computeAll)
	 * \param[in] l_localface The local face number of this face as seen from the left element
	 * \param[in] r_localface The local face number of this face as seen from the right element
	 */
	void initialize(const Element *const lelem, const Element *const relem, const GeomMapping1D *const geommap, const int l_localface, const int r_localface);
	
	/// Computes gradients of the left- and right-elements' basis functions at face quadrature points
	/** To be called only after [initializing](@ref initialize) the face element.
	 * TODO: Should probably merge this with initialize.
	 */
	void computeBasisGrads();

	/// Read-only access to basis function values from left element
	const Matrix& leftBasis() {
		return leftbasis;
	}

	/// Read-only access to basis function values from right element
	const Matrix& rightBasis() {
		return rightbasis;
	}

	/// Read access to left basis gradients
	const std::vector<Matrix>& leftBasisGrad() {
		return leftbgrad;
	}

	/// Read access to right basis gradients
	const std::vector<Matrix>& rightBasisGrad() {
		return rightbgrad;
	}

	/// Interpolates values from left element at the face quadrature points
	a_real interpolate_left(const int ig, const Vector& dofs)
	{
		a_real val = 0;
		for(int i = 0; i < leftel->getNumDOFs(); i++)
			val += dofs[i]*leftbasis(ig,i);
		return val;
	}

	/// Interpolates values from right element at the face quadrature points
	a_real interpolate_right(const int ig, const Vector& dofs)
	{
		a_real val = 0;
		for(int i = 0; i < rightel->getNumDOFs(); i++)
			val += dofs[i]*rightbasis(ig,i);
		return val;
	}

	void interpolateAll_left(const Matrix& __restrict__ dofs, Matrix& __restrict__ values) {
		values.noalias() = leftbasis*dofs.transpose();
	}

	void interpolateAll_right(const Matrix& __restrict__ dofs, Matrix& __restrict__ values) {
		values.noalias() = rightbasis*dofs.transpose();
	}
};

}
#endif
