/** @file aspatial.hpp
 * @brief Residual and Jacobian computations
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#ifndef ASPATIAL_H
#define ASPATIAL_H 1

#include "aconstants.hpp"
#include "aarray2d.hpp"
#include "amesh2dh.hpp"
#include "aelements.hpp"
#include "areconstruction.hpp"

#include <Eigen/LU>

namespace acfd {

/// Base class for spatial discretization and integration of weak forms of PDEs
/**
 * Provides residual computation, and potentially residual Jacobian evaluation, interface for all solvers.
 * The template parameter nvars is the number of variables in the PDE system.
 * \note Make sure compute_topological() has been called on the mesh object prior to initialzing 
 * an object of any subclass.
 */
template <short nvars>
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

	/* Reconstruction-related stuff - currently not implemented
	//Reconstruction* rec;						///< Reconstruction context
	//FaceDataComputation* lim;					///< Limiter context

	// Ghost cell centers
	amat::Array2d<a_real> ghc;

	/// Ghost elements' flow quantities
	std::vector<Vector> ug;

	/// computes ghost cell centers assuming symmetry about the midpoint of the boundary face
	void compute_ghost_cell_coords_about_midpoint();

	/// computes ghost cell centers assuming symmetry about the face
	void compute_ghost_cell_coords_about_face();
	*/

	/// Sets up geometric maps, elements and mass matrices 
	void computeFEData();
	
	/// Computes the L2 error in a FE function on an element
	/** \param[in] comp The index of the row of ug whose error is to be computed
	 */
	a_real computeElemL2Error2(const int ielem, const int comp, const Matrix& ug, 
			a_real (* const exact)(a_real, a_real, a_real), const double time) const;

	/// Computes the L2 norm of a FE function on an element
	a_real computeElemL2Norm2(const int ielem, const Vector& ug) const;

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
	void spatialSetup(std::vector<Matrix>& u, std::vector<Matrix>& res, std::vector<a_real>& mets);

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

	/// Adds source term contribution to residual
	/** As implemented in this class, does nothing.
	 */
	virtual void add_source( a_real (*const rhs)(a_real, a_real, a_real), a_real t, std::vector<Matrix>& res);

	/// Compute quantities to export
	virtual void postprocess(const std::vector<Matrix>& u) = 0;

	/// Read-only access to output quantities
	virtual const amat::Array2d<a_real>& getOutput() const = 0;

	/// Computes the norm of the difference between a FE solution and an analytically defined function
	a_real computeL2Error(double (*const exact)(double,double,double), const double time, 
			const std::vector<Matrix>& u) const;

	/// Sets initial conditions using a function describing a variable
	void setInitialConditionNodal( const int comp, double (**const init)(a_real, a_real), 
			std::vector<Matrix>& u);

	/// Sets initial conditions using functions for a variable and its space derivatives
	void setInitialConditionModal( const int comp, double (**const init)(a_real, a_real), 
			std::vector<Matrix>& u);
};

}	// end namespace
#endif
