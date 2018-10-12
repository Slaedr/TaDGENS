/** @file atimetvdrk.hpp
 * @brief Explicit total variation diminishing Runge-Kutta (TVDRK) time stepping scheme
 * @author Aditya Kashi
 * @date 2017 April 15
 */

#ifndef ATIMETVDRK_H
#define ATIMETVDRK_H

#include "spatial/aspatial.hpp"
#include "aodecoeffs.hpp"

namespace acfd {

/// TVD RK explicit time stepping
/** The initial condition must be specified in the spatial discretization context 
 * [elsewhere](@ref SpatialBase::setInitialConditionModal).
 */
class TVDRKStepping
{
protected:
	const UMesh2dh *const m;						///< Mesh context
	SpatialBase* spatial;							///< Spatial discretization context
	int order;										///< Desird temporal order of accuracy
	double cfl;										///< CFL number
	double ftime;									///< Physical time up to which simulation should proceed

	/// 'c' or 'a' for constant or automatic (non-constant) time steps respectively
	char tch;										
	double timestep;								///< Fixed time step, if tch was 'c'

public:
	TVDRKStepping(const UMesh2dh*const mesh, SpatialBase *const s, const int timeorder,
	              a_real final_time, a_real cflnumber,
	              const char tc, const double time_step);

	void setInitialCondition();

	/// Carries out the time stepping process and returns the final time
	double integrate();
	
	double integrate_ForwardEuler();
};

}
#endif
