/** @file atimesteady.hpp
 * @brief Steady-state pesudo-time stepping schemes
 * @author Aditya Kashi
 * @date 2017 April 18
 */

#ifndef __ATIMESTEADY_H
#define __ATIMESTEADY_H

#ifndef __ASPATIAL_H
#include "aspatial.hpp"
#endif

#ifndef __AODECOEFFS_H
#include "aodecoeffs.hpp"
#endif

namespace acfd {
	
/// TVD RK explicit time stepping
/** The initial condition must be specified in the spatial discretization context [elsewhere](@ref SpatialBase::setInitialConditionModal).
 */
class SteadyBase
{
protected:
	const UMesh2dh *const m;						///< Mesh context
	SpatialBase* spatial;							///< Spatial discretization context
	int order;										///< Desird temporal order of accuracy
	double cfl;										///< CFL number
	double tol;										///< Tolerance for residual
	int maxiter;									///< Max number of iterations
	bool source;									///< Whether or not to use source term
	a_real (*const rhs)(a_real, a_real, a_real);	///< Function describing source term

public:
	SteadyBase(const UMesh2dh*const mesh, SpatialBase* s, a_real cflnumber, double toler, int max_iter, bool use_source);
	
	void set_source( a_real (*const source)(a_real, a_real, a_real)) {
		rhs = source;
	}

	/// Carries out the time stepping process
	virtual void integrate() = 0;
};

/// Explicit forward-Euler scheme with local time stepping
class SteadyExplicit : public SteadyBase
{
public:
	SteadyExplicit(const UMesh2dh*const mesh, SpatialBase* s, a_real cflnumber, double toler, int max_iter, bool use_source);

	/// Carries out the time stepping process
	void integrate();
};

}
#endif
