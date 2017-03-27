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

/// Base class for spatial discretization and integration
/** 
 * Provides residual computation, and potentially residual Jacobian evaluation, interface for all solvers.
 * \note Make sure compute_topological(), compute_face_data() and compute_jacobians() have been called on the mesh object prior to initialzing an object of any subclass.
 */
class SpatialBase
{
protected:
	const UMesh2dh* m;							///< Mesh context
	std::vector<Matrix> m_inv;					///< Inverse of mass matrix
	std::vector<Vector> residual;				///< Right hand side for boundary integrals and source terms
	int nvars;									///< number of conserved variables ** deprecated, use the preprocessor constant NVARS instead **
	int p_degree;								///< Polynomial degree of trial/test functions
	Vector uinf;								///< Free-stream/reference condition
	acfd_real g;								///< adiabatic index

	/// stores (for each cell i) \f$ \sum_{j \in \partial\Omega_I} \int_j( |v_n| + c) d \Gamma \f$, where v_n and c are average values of the cell faces
	std::vector<acfd_real> integ;
	
	InviscidFlux* inviflux;						///< Flux (boundary integral) calculation context
	//Reconstruction* rec;						///< Reconstruction context
	//FaceDataComputation* lim;					///< Limiter context
	Quadrature2D* dquad;						///< Domain quadrature context
	Quadrature1D* bquad;						///< Boundary quadrature context
	Element_PhysicalElement* elems;				///< List of finite elements
	FaceElement_PhysicalSpace* faces;			///< List of face elements

	// Cell centers
	//amat::Array2d<acfd_real> rc;

	/// Ghost element centers
	amat::Array2d<acfd_real> rcg;
	/// Ghost elements' flow quantities
	std::vector<Vector> ug;

	/// Integral of fluxes across each face for all dofs
	/** The entries corresponding to different DOFs of a given flow variable are stored contiguously.
	 */
	std::vector<Array2d<acfd_real>> faceintegral;
	/// Left state at each face
	amat::Array2d<acfd_real> uleft;
	/// Rigt state at each face
	amat::Array2d<acfd_real> uright;

	/// vector of unknowns
	/** Each Eigen3 (E3) Vector contains the DOF values for an element.
	 */
	std::vector<Vector> u;


	amat::Array2d<acfd_real> scalars;			///< Holds density, Mach number and pressure for each mesh point
	amat::Array2d<acfd_real> velocities;		///< Holds velocity components for each mesh point

	int slip_wall_id;							///< Boundary marker corresponding to solid wall
	int inflow_outflow_id;						///< Boundary marker corresponding to inflow/outflow
	int periodic_id;							///< Boundary marker for periodic boundary

	/// computes ghost cell centers assuming symmetry about the midpoint of the boundary face
	void compute_ghost_cell_coords_about_midpoint();

	/// computes ghost cell centers assuming symmetry about the face
	void compute_ghost_cell_coords_about_face();

	/// Computes the left and right states at each face
	/** If applicable, the [reconstruction](@ref rec) and [limiter](@ref limiter) objects are used too.
	 */
	void compute_face_states();

	/// Add contribution of inviscid numerical flux to flux storage
	void inviscidFluxContribution();

public:
	SpatialBase(const UMesh2dh* mesh, const int _p_degree, const InviscidFlux* invflux/*, const Reconstruction* reconst*/);
	~SpatialBase();

	/// Computes flow variables at boundaries (either Gauss points or ghost cell centers) using the interior state provided
	/** \param[in] instates provides the left (interior state) for each boundary face
	 * \param[out] bounstates will contain the right state of boundary faces
	 *
	 * Currently does not use characteristic BCs.
	 * \todo Implement and test characteristic BCs
	 */
	void compute_boundary_states(const std::vector<Vector>& instates, std::vector<Vector>& bounstates);

	/// Calls functions to assemble the [right hand side](@ref residual)
	void compute_residual();

	/// Computes the L2 norm of a cell-centered quantity
	acfd_real l2norm(const amat::Array2d<acfd_real>* const v);

	/// Access to vector of unknowns
	std::vector<Vector>& unk();

	/// Access to residual vector
	std::vector<Vector>& res();

	/// Mass matrix
	std::vector<Matrix>& mass();

	/// Compute cell-centred quantities to export
	void postprocess_cell();

	/// Compute nodal quantities to export, based on area-weighted averaging (which takes into account ghost cells as well)
	void postprocess_point();

	/// Compute norm of entropy production
	/// Call after computing pressure etc \sa postprocess_cell
	acfd_real compute_entropy();

	amat::Array2d<acfd_real> getscalars() const;
	amat::Array2d<acfd_real> getvelocities() const;
};

class InviscidFlow : public SpatialBase
{

};

}	// end namespace
#endif
