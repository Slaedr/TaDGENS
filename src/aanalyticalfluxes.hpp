/** @file aanalyticalfluxes.hpp
 * @brief Classes for computation of analytical flux function vectors
 * @author Aditya Kashi
 * @date 2016, 2017-03-10
 */

#ifndef __AANALYTICALFLUXES_H
#define __AANALYTICALLUXES_H

#ifndef __ACONSTANTS_H
#include "aconstants.hpp"
#endif

#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

namespace acfd {

/// Computation of the single-phase ideal gas Euler flux corresponding to any given state and along any given face-normal
class NormalEulerFluxFunction
{
protected:
	const a_real gamma;
public:
	NormalEulerFluxFunction (a_real _gamma) : gamma(_gamma)
	{ }

	void evaluate_flux(const amat::Array2d<a_real>& state, const a_real* const n, amat::Array2d<a_real>& flux) const
	{
		a_real vn = (state.get(1)*n[0] + state.get(2)*n[1])/state.get(0);
		a_real p = (gamma-1.0)*(state.get(3) - 0.5*(state.get(1)*state.get(1) + state.get(2)*state.get(2))/state.get(0));
		flux(0) = state.get(0) * vn;
		flux(1) = vn*state.get(1) + p*n[0];
		flux(2) = vn*state.get(2) + p*n[1];
		flux(3) = vn*(state.get(3) + p);
	}
	
	void evaluate_flux_2(const amat::Array2d<a_real>& state, const int ielem, const a_real* const n, amat::Array2d<a_real>& flux, const int iside) const
	{
		a_real vn = (state.get(ielem,1)*n[0] + state.get(ielem,2)*n[1])/state.get(ielem,0);
		a_real p = (gamma-1.0)*(state.get(ielem,3) - 0.5*(state.get(ielem,1)*state.get(ielem,1) + state.get(ielem,2)*state.get(ielem,2))/state.get(ielem,0));
		flux(iside,0) = state.get(ielem,0) * vn;
		flux(iside,1) = vn*state.get(ielem,1) + p*n[0];
		flux(iside,2) = vn*state.get(ielem,2) + p*n[1];
		flux(iside,3) = vn*(state.get(ielem,3) + p);
	}
};

/// Abstract context for computation of analytical flux vectors
class FluxFunction
{
public:
	virtual void evaluate_flux(const a_real* const u, a_real* const fx, a_real* const fy) = 0;
};

/// Sort of a functor for evaluation of analytical inviscid flux vector
class EulerFlux : public FluxFunction
{
	const a_real g;				///< Adiabatic index
public:
	EulerFlux (a_real _gamma) : g(_gamma)
	{ }
	void evaluate_flux(const a_real *const __restrict__ u, a_real *const __restrict__ fx, a_real *const __restrict__ fy)
	{
		a_real p = (g-1)*(u[3] - 0.5*(u[1]*u[1]+u[2]*u[2])/u[0]);
		fx[0] = u[1];
		fx[1] = u[1]*u[1]/u[0] + p;
	}
};

}
#endif
