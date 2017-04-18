/** @file atimetvdrk.hpp
 * @brief Explicit total variation diminishing Runge-Kutta (TVDRK) time stepping scheme
 * @author Aditya Kashi
 * @date 2017 April 15
 */

#ifndef __ATIMETVDRK_H
#define __ATIMETVDRK_H

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
class TVDRKStepping
{
protected:
	const UMesh2dh *const m;						///< Mesh context
	SpatialBase* spatial;							///< Spatial discretization context
	int order;										///< Desird temporal order of accuracy
	double cfl;										///< CFL number
	double ftime;									///< Physical time up to which simulation should proceed
	char tch;										///< 'c' or 'a' for constant or automatic (non-constant) time steps respectively
	double timestep;								///< Fixed time step, if tch was 'c'

public:
	TVDRKStepping(const UMesh2dh*const mesh, SpatialBase* s, const int timeorder, a_real final_time, a_real cflnumber,
			const char tc, const double time_step);
	void setInitialCondition();

	/// Carries out the time stepping process and returns the final time
	double integrate();
	
	double integrate_ForwardEuler();
};

TVDRKStepping::TVDRKStepping(const UMesh2dh *const mesh, SpatialBase* s, const int timeorder, a_real final_time, a_real cflnumber, 
		const char tc, const double time_step)
	: m(mesh), spatial(s), order(timeorder), cfl(cflnumber), ftime(final_time), tch(tc), timestep(time_step)
{
}

double TVDRKStepping::integrate()
{
	int step = 0; double time = 0; double tsg = timestep;
	std::vector<Matrix>& u = spatial->unk();
	
	std::vector<Matrix> ustage(m->gnelem());
	for(int iel = 0; iel < m->gnelem(); iel++)
		ustage[iel] = Matrix::Zero(u[iel].rows(), u[iel].cols());

	std::vector<Matrix>& R = spatial->residual();
	const std::vector<Matrix>& Mi = spatial->massInv();
	const std::vector<a_real>& tsl = spatial->maxExplicitTimeStep();
	std::printf(" TVDRKStepping: integrate: Time step = %f, option = %c, order = %d\n", tsg, tch, order);
	initializeOdeCoeffs();
	amat::Array2d<a_real> tvdrk;
	if(order == 1)
		tvdrk = tvdrk1;
	else if(order == 2)
		tvdrk = tvdrk2;
	else if(order == 3)
		tvdrk = tvdrk3;
	else {
		std::printf(" TVDRKStepping: integrate: Order not supported! Using 3.\n");
		tvdrk = tvdrk3;
		order = 3;
	}
	
	for(int iel = 0; iel < m->gnelem(); iel++)
		ustage[iel] = u[iel];
	
	while(time < ftime-SMALL_NUMBER)
	{
		for(int istage = 0; istage < order; istage++)
		{
			for(int iel = 0; iel < m->gnelem(); iel++) {
				R[iel] = Matrix::Zero(R[iel].rows(), R[iel].cols());
			}

			spatial->update_residual(ustage);
			
			if(istage == 0) {
				// get global time step
				if(tch == 'a') {
					for(int iel = 0; iel < m->gnelem(); iel++) {
						if(tsl[iel] < tsg)
							tsg = tsl[iel];
					}
					tsg = cfl*tsg;
				}
			}

			// step
			for(int iel = 0; iel < m->gnelem(); iel++)
			{
				ustage[iel] = tvdrk[istage][0]*u[iel] + tvdrk[istage][1]*ustage[iel] - tvdrk[istage][2] * tsg*R[iel]*Mi[iel];
			}
		}

		for(int iel = 0; iel < m->gnelem(); iel++)
			u[iel] = ustage[iel];

		time += tsg; step++;
		if(step % 20 == 0)
			std::printf("  TVDRKStepping: integrate: Step %d, time = %f\n", step, time);
	}
	return time;
}

double TVDRKStepping::integrate_ForwardEuler()
{
	int step = 0; double time = 0; double tsg = timestep;
	std::vector<Matrix>& u = spatial->unk();

	std::vector<Matrix>& R = spatial->residual();
	const std::vector<Matrix>& Mi = spatial->massInv();
	const std::vector<a_real>& tsl = spatial->maxExplicitTimeStep();
	std::printf(" TVDRKStepping: integrate: Time step = %f, option = %c, order = %d\n", tsg, tch, order);

	// compare P1 area inverses
	std::vector<a_real> mi(m->gnelem());
	for(int i = 0; i < m->gnelem(); i++)
	{
		mi[i] = 2.0/( m->gcoords(m->ginpoel(i,0),0)*(m->gcoords(m->ginpoel(i,1),1) - m->gcoords(m->ginpoel(i,2),1)) 
				- m->gcoords(m->ginpoel(i,0),1)*(m->gcoords(m->ginpoel(i,1),0)-m->gcoords(m->ginpoel(i,2),0)) 
				+ m->gcoords(m->ginpoel(i,1),0)*m->gcoords(m->ginpoel(i,2),1) - m->gcoords(m->ginpoel(i,2),0)*m->gcoords(m->ginpoel(i,1),1) );
	}

	while(time < ftime-SMALL_NUMBER)
	{
		for(int iel = 0; iel < m->gnelem(); iel++) {
			R[iel] = Matrix::Zero(R[iel].rows(), R[iel].cols());
		}

		spatial->update_residual(u);
	
		// get global time step
		if(tch == 'a') {
			for(int iel = 0; iel < m->gnelem(); iel++) {
				if(tsl[iel] < tsg)
					tsg = tsl[iel];
			}
			tsg = cfl*tsg;
		}

		// step
		for(int iel = 0; iel < m->gnelem(); iel++)
		{
			u[iel](0,0) = u[iel](0,0) - tsg*R[iel](0,0)*Mi[iel](0,0);
		}

		time += tsg; step++;
		if(step % 20 == 0)
			std::printf("  TVDRKStepping: integrate: Step %d, time = %f\n", step, time);
	}
	return time;
}

}
#endif
