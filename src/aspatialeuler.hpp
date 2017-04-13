/** @file aeulerspatial.hpp
 * @brief Spatial discretization ofr Euler equations
 * @auhor Aditya Kashi
 * @date 2017 April 14
 */


#ifndef __ANUMERICALFLUX_H
#include "anumericalflux.hpp"
#endif

/// Spatial discretization for 2D Euler equations
/** Most functions are virtual so that other solvers can be subclassed
 * and the functionality here can be selectively used.
 */
class EulerFlow : public SpatialBase
{
protected:
	a_real g;										///< Adiabatic index
	Vector uinf;									///< Free-stream/reference condition
	Vector uin;										///< Inflow condition
	Vector uout;									///< Outflow condition
	int slipwall_id;								///< Boundary marker corresponding to solid wall
	int inflow_id;									///< Boundary marker corresponding to inflow
	int outflow_id;									///< Boundary marker corresponding to outflow
	int farfield_id;								///< Boundary marker corresponding to far field
	int periodic_id;								///< Boundary marker for "periodic" boundary
	int symmetry_id;								///< Boundary marker for "symmetry" boundary

	/// Add contribution of inviscid numerical flux to flux storage
	void inviscidFluxContribution();

	/// Computes flow variables at boundaries (either Gauss points or ghost cell centers) using the interior state provided
	/** \param[in] ins provides the left (interior state) for each boundary face
	 * \param[in] n the unit normal vector to the face
	 * \param[in] iface The face index in the [face structure](@ref UMesh2dh::intfac)
	 * \param[out] bs will contain the right state of boundary faces
	 *
	 * Currently does not use characteristic BCs.
	 * \todo Implement and test characteristic BCs
	 */
	virtual void compute_boundary_states(const a_real ins[NVARS], const Vector& n, int iface, a_real bs[NVARS]);

public:
	/// Constructor
	/** \param[in] mesh is the mesh context
	 * \param _p_degree is the polynomial degree for FE basis functions
	 * \param gamma Adiabatic index
	 * \param u_inf Farfield conditions - density, velocity magnitude, angle of attack and Mach number
	 * \param u_in Inlet conditions - total pressure, total temperature, angle of attack and Mach number
	 * \param u_out Outflow conditions - static pressure
	 * \param boun_ids Boundary flags for the following types of boundary conditions in that order:
	 * 0: slip wall, 1: inflow, 2: outflow, 3: far field, 4: periodic boundary, 5: symmetry boundary.
	 * Some simulations may not require one or more of these, in which case dummy unallocated Vectors should be passed.
	 */
	EulerFlow(const UMesh2dh* mesh, const int _p_degree, const a_real gamma, Vector& u_inf, Vector& u_in, Vector& u_out, int boun_ids[6]);

	/// Calls functions to update the [right hand side](@ref residual)
	virtual void update_residual();

	/// Compute norm of entropy production
	a_real compute_entropy();

	/// Compute nodal quantities to export, based on area-weighted averaging (which takes into account ghost cells as well)
	virtual void postprocess();

	/// Read-only access to output quantities
	virtual const amat::Array2d<a_real>& getoutput() const;
};
