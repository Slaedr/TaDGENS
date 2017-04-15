/** @file aspatialadvection.hpp
 * @brief Spatial discretization for linear advection
 * @author Aditya Kashi
 * @date 2017 April 14
 */

#ifndef __ASPATIALADVECTION_H
#define __ASPATIALADVECTION_H

#include "aspatial.hpp"

namespace acfd {

/// Residual computation for linear advection
/** \note Make sure to call both compute_topological and compute_boundary_maps on the mesh object before using an object of this class!
 */
class LinearAdvection : public SpatialBase
{
protected:
	Vector a;								///< Advection velocity
	a_real bval;							///< Value to be imposed at inflow
	int inoutflow_flag;						///< Boundary flag at faces where inflow or outflow is required
	int extrapolation_flag;					///< Boundary flag for extrapolation condition

	void computeFaceTerms();				///< Computes face integrals
	void computeDomainTerms();				///< Computes domain integrals

	/// Computes boundary (ghost) states depending on face marker
	void computeBoundaryState(const int iface, const Vector& instate, Vector& bstate);

public:
	LinearAdvection(const UMesh2dh* mesh, const int _p_degree, const char basis, const Vector vel, const a_real b_val, const int inoutflag, const int extrapflag);
	
	/// Calls functions to add contribution to the [right hand side](@ref residual)
	void update_residual();

	/// Compute quantities to export
	void postprocess();

	/// Read-only access to output quantities
	const amat::Array2d<a_real>& getOutput() const;
};

}
#endif
