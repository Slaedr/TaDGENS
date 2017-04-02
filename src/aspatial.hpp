/** @file aspatial.hpp
 * @brief Residual and Jacobian computations
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#ifndef __ASPATIAL_H
#define __ASPATIAL_H 1


#ifndef __ACONSTANTS_H
#include "aconstants.hpp"
#endif

#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

#ifndef __AMESH2DH_H
#include "amesh2dh.hpp"
#endif

#ifndef __ANUMERICALFLUX_H
#include "anumericalflux.hpp"
#endif

#ifndef __AELEMENTS_H
#include "aelements.hpp"
#endif

namespace acfd {

/// Base class for spatial discretization and integration of weak forms of PDEs
/**
 * Provides residual computation, and potentially residual Jacobian evaluation, interface for all solvers.
 * The template parameter nvars is the number of variables in the PDE system.
 * \note Make sure compute_topological() has been called on the mesh object prior to initialzing an object of any subclass.
 */
template <int nvars>
class SpatialBase
{
protected:
	const UMesh2dh* m;									///< Mesh context
	std::vector<std::array<Matrix, nvars>> m_inv;		///< Inverse of mass matrix for each variable of each element
	std::vector<Vector> residual;						///< Right hand side for boundary integrals and source terms
	int p_degree;										///< Polynomial degree of trial/test functions

	/// stores (for each cell i) \f$ \sum_{j \in \partial\Omega_I} \int_j( |v_n| + c) d \Gamma \f$, where v_n and c are average values of the cell faces
	std::vector<acfd_real> integ;

	Quadrature2DTriangle* dtquad;				///< Domain quadrature context
	Quadrature2DSquare* dsquad;					///< Domain quadrature context
	Quadrature1D* bquad;						///< Boundary quadrature context
	LagrangeMapping2D* map2d;					///< Array containing geometric mapping data for each element
	LagrangeMapping1D* map1d;					///< Array containing geometric mapping data for each face
	Element_PhysicalElement* elems;				///< List of finite elements
	FaceElement_PhysicalSpace* faces;			///< List of face elements

	/// Integral of fluxes across each face for all dofs
	/** The entries corresponding to different DOFs of a given flow variable are stored contiguously.
	 */
	std::vector<Array2d<acfd_real>> faceintegral;

	/// vector of unknowns
	/** Each Eigen3 (E3) Vector contains the DOF values for an element.
	 */
	std::vector<Vector> u;

	amat::Array2d<acfd_real> scalars;			///< Holds density, Mach number and pressure for each mesh point
	amat::Array2d<acfd_real> velocities;		///< Holds velocity components for each mesh point

	/* Reconstruction-related stuff - currently not implemented
	//Reconstruction* rec;						///< Reconstruction context
	//FaceDataComputation* lim;					///< Limiter context

	// Ghost cell centers
	amat::Array2d<acfd_real> ghc;

	/// Ghost elements' flow quantities
	std::vector<Vector> ug;

	/// computes ghost cell centers assuming symmetry about the midpoint of the boundary face
	void compute_ghost_cell_coords_about_midpoint();

	/// computes ghost cell centers assuming symmetry about the face
	void compute_ghost_cell_coords_about_face();
	*/

public:
	/// Constructor
	/** \param[in] mesh is the mesh context
	 * \param _p_degree is the polynomial degree for FE basis functions
	 */
	SpatialBase(const UMesh2dh* mesh, const int _p_degree);

	virtual ~SpatialBase();

	/// Compute all finite element data, including mass matrix, needed for the spatial discretization
	void computeFEData();

	/// Computes ghost states of the PDE variables at boundary quadrature points
	virtual void compute_boundary_states(const std::vector<Vector>& instates, std::vector<Vector>& bounstates) = 0;

	/// Calls functions to add contribution to the [right hand side](@ref residual)
	virtual void update_residual() = 0;

	/// Access to vector of unknowns
	std::vector<Vector>& unk() {
		return u;
	}

	/// Access to residual vector
	std::vector<Vector>& res() {
		return residual;
	}

	/// Mass matrix
	const std::vector<std::array<Matrix, nvars>>& mass() {
		return m_inv;
	}

	/// Compute quantities to export
	virtual void postprocess() = 0;

	/// Read-only access to output quantities
	virtual const amat::Array2d<acfd_real>& getoutput() const = 0;
};

/// Spatial discretization for 2D Euler equations
/** Most functions are virtual so that other solvers can be subclassed
 * and the functionality here can be selectively used.
 */
class InviscidFlow : public SpatialBase<4>
{
protected:
	Vector uinf;										///< Free-stream/reference condition
	int slip_wall_id;									///< Boundary marker corresponding to solid wall
	int inflow_outflow_id;								///< Boundary marker corresponding to inflow/outflow
	int periodic_id;									///< Boundary marker for periodic boundary

	/// Add contribution of inviscid numerical flux to flux storage
	void inviscidFluxContribution();

	/// Computes flow variables at boundaries (either Gauss points or ghost cell centers) using the interior state provided
	/** \param[in] instates provides the left (interior state) for each boundary face
	 * \param[out] bounstates will contain the right state of boundary faces
	 *
	 * Currently does not use characteristic BCs.
	 * \todo Implement and test characteristic BCs
	 */
	virtual void compute_boundary_states(const std::vector<Vector>& instates, std::vector<Vector>& bounstates);

public:
	/// Calls functions to assemble the [right hand side](@ref residual)
	virtual void compute_residual();

	/// Compute norm of entropy production
	/// Call after computing pressure etc \sa postprocess_cell
	acfd_real compute_entropy();

	/// Compute nodal quantities to export, based on area-weighted averaging (which takes into account ghost cells as well)
	virtual void postprocess();

	/// Read-only access to output quantities
	virtual const amat::Array2d<acfd_real>& getoutput() const;
};

}	// end namespace
#endif
