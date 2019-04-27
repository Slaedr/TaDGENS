/** @file aspatial.hpp
 * @brief Residual and Jacobian computations
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#ifndef ASPATIAL_H
#define ASPATIAL_H 1

#include "aconstants.hpp"
#include "utilities/aarray2d.hpp"
#include "mesh/amesh2dh.hpp"
#include "fem/aelements.hpp"
#include "areconstruction.hpp"

#include <Eigen/LU>

namespace acfd {

/// Base class for spatial discretization and integration of weak forms of PDEs
/**
 * Provides residual computation, and potentially residual Jacobian evaluation, interface for all solvers.
 * \note Make sure compute_topological() has been called on the mesh object prior to initialzing 
 * an object of any subclass.
 */
class SpatialBase
{
protected:
	/// Mesh context; requires compute_topological() and compute_boundary_maps() to have been called
	const UMesh2dh* m;

	std::vector<Matrix> minv;             ///< Inverse of mass matrix for each variable of each element
	int p_degree;                         ///< Polynomial degree of trial/test functions
	a_int ntotaldofs;                     ///< Total number of DOFs in the discretization per physical variable)
	char basis_type;                      ///< Type of basis to use - Lagrange ('l') or Taylor ('t')
	bool reconstruct;                     ///< Use reconstruction or not

	Quadrature2DTriangle* dtquad;				///< Domain quadrature context
	Quadrature2DSquare* dsquad;					///< Domain quadrature context
	Quadrature1D* bquad;						///< Boundary quadrature context
	LagrangeMapping2D* map2d;					///< Array containing geometric mapping data for each element
	LagrangeMapping1D* map1d;					///< Array containing geometric mapping data for each face
	Element** elems;							///< List of finite elements
	Element* dummyelem;							///< Empty element used for ghost elements
	FaceElement* faces;							///< List of face elements

	amat::Array2d<a_real> scalars;				///< Holds scalar variables for each mesh point
	amat::Array2d<a_real> velocities;			///< Holds velocity components for each mesh point

	/// Sets up geometric maps, elements and mass matrices 
	void computeFEData();
	
	/// Computes the L2 error in a FE function on an element
	/** \param[in] comp The index of the row of ug whose error is to be computed
	 */
	a_real computeElemL2Error2(const int ielem, const int comp, const Matrix& ug, 
			const double time) const;

	/// Computes the L2 norm of a FE function on an element
	a_real computeElemL2Norm2(const int ielem, const Vector& ug) const;

	/// Intended to provide a test source term for a verification case
	/** This does not need to be virtual as it will only be used in derived classes.
	 * The implemenatation in this base class does nothing.
	 */
	virtual a_real source_term(const a_real position[NDIM], const a_real time) const = 0;

public:
	/// Constructor
	/** \param[in] mesh is the mesh context
	 * \param _p_degree is the polynomial degree for FE basis functions
	 */
	SpatialBase(const UMesh2dh* mesh, const int _p_degree, char basistype);

	virtual ~SpatialBase();

	/// Allocate arrays and precompute some finite-element data
	/** Compute all finite element data, including mass matrix, 
	 * and allocates solution, residual and time-step arrays
	 */
	virtual void spatialSetup(std::vector<Matrix>& u, std::vector<Matrix>& res,
	                          std::vector<a_real>& mets);

	/// Computes L2 norm of the the specified component of some vector quantity w
	a_real computeL2Norm(const std::vector<Matrix> w, const int comp) const;

	/// Inverse of mass matrix
	const std::vector<Matrix>& massInv() const {
		return minv;
	}

	a_int numTotalDOFs() const { return ntotaldofs; }

	/// Calls functions to add contribution to the RHS, and also compute max time steps
	virtual void update_residual(const std::vector<Matrix>& u, 
	                             std::vector<Matrix>& res, 
	                             std::vector<a_real>& mets) = 0;

	/// Compute quantities to export
	virtual void postprocess(const std::vector<Matrix>& u) = 0;

	/// Read-only access to output quantities
	virtual const amat::Array2d<a_real>& getOutput() const = 0;

	/// Intended to provide the exact solution for a verification case
	virtual a_real exact_solution(const a_real position[NDIM], const a_real time) const = 0;

	/// Computes the norm of the difference between a FE solution and an analytically defined function
	/** \param exact A function that takes (x,y,t) and returns the solution at that point
	 */
	a_real computeL2Error(const double time, const std::vector<Matrix>& u) const;

	/// Sets initial conditions using a function describing a variable
	void setInitialConditionNodal( const int comp, double (**const init)(a_real, a_real), 
	                               std::vector<Matrix>& u);

	/// Sets initial conditions using functions for a variable and its space derivatives
	void setInitialConditionModal( const int comp, double (**const init)(a_real, a_real), 
	                               std::vector<Matrix>& u);
};

}	// end namespace
#endif
