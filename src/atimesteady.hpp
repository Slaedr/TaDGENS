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
	a_real (*rhs)(a_real, a_real, a_real);			///< Function describing source term

public:
	SteadyBase(const UMesh2dh*const mesh, SpatialBase* s, a_real cflnumber, double toler, int max_iter, bool use_source);
	
	/// Sets the forcing function for the source term
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
	/** \param[in] mesh The mesh context
	 * \param[in] s The spatial discretization context
	 * \param[in] cflnumber
	 * \param[in] toler Tolerance for the relative residual
	 * \param[in] max_iter Maximum number of iterations
	 * \param[in] use_source True if source term integration is required
	 */
	SteadyExplicit(const UMesh2dh*const mesh, SpatialBase* s, a_real cflnumber, double toler, int max_iter, bool use_source);

	/// Carries out the time stepping process
	void integrate();
};

}
#endif
