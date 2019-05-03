/** @file atimesteady.cpp
 * @brief Implementation of steady-state solvers
 * @author Aditya Kashi
 * @date 2017 April 15
 */

#include <iostream>
#include "atimesteady.hpp"

namespace tadgens {

SteadyBase::SteadyBase(const UMesh2dh *const mesh, SpatialBase *const s,
                              a_real cflnumber, double toler, int max_iter)
	: m(mesh), spatial(s), cfl(cflnumber), tol(toler), maxiter(max_iter)
{
	std::cout << " SteadyBase: CFL = " << cfl << std::endl;

	spatial->spatialSetup(u, R, tsl);

	for(a_int iel = 0; iel < m->gnelem(); iel++)
	{
		for(int i = 0; i < u[iel].rows(); i++)
			for(int j = 0; j < u[iel].cols(); j++) {
				u[iel](i,j) = 1.0;
				R[iel](i,j) = 0;
			}
	}
}

SteadyExplicit::SteadyExplicit(const UMesh2dh*const mesh, SpatialBase *const s,
                                      a_real cflnumber, double toler, int max_iter)
	: SteadyBase(mesh, s, cflnumber, toler, max_iter)
{
}

void SteadyExplicit::solve()
{
	int step = 0;
	double relresnorm = 1.0, resnorm0 = 1.0;
	const std::vector<Matrix>& Mi = spatial->massInv();

	while((relresnorm > tol && step < maxiter))
	{
		for(int iel = 0; iel < m->gnelem(); iel++) {
			for(int i = 0; i < R[iel].rows(); i++)
				for(int j = 0; j < R[iel].cols(); j++)
					R[iel](i,j) = 0.0;
		}

		spatial->update_residual(u, R, tsl);

		// step
		for(int iel = 0; iel < m->gnelem(); iel++)
		{
			u[iel] = u[iel] - cfl*tsl[iel]*R[iel]*Mi[iel];
		}

		//double resnorm = spatial->computeL2Norm(R, 0);
		double resnorm = 0;
		for(int iel = 0; iel < m->gnelem(); iel++) {
			for(int j = 0; j < R[iel].rows(); j++)
			{
				resnorm += R[iel](0,j)*R[iel](0,j);
			}
		}
		resnorm = std::sqrt(resnorm);
		if(step == 0) resnorm0 = resnorm;
		else relresnorm = resnorm/resnorm0;

		step++;
		if(step % 50 == 0)
			std::cout << "  SteadyExplicit: integrate: Step " << step
			          << ", res = " << resnorm << ", rel res = " << relresnorm << std::endl;
	}

	std::printf(" SteadyExplicit: integrate: Total steps %d, final rel res = %e\n", step, relresnorm);
}

SteadyImplicit::SteadyImplicit(const UMesh2dh*const mesh, SpatialBase *const s,
                                      a_real cflnumber, double toler, int max_iter)
	: SteadyBase(mesh, s, cflnumber, toler, max_iter)
{
}

void SteadyImplicit::solve()
{
	int step = 0;
	double relresnorm = 1.0, resnorm0 = 1.0;
	const std::vector<Matrix>& Mi = spatial->massInv();

	while((relresnorm > tol && step < maxiter))
	{
		for(int iel = 0; iel < m->gnelem(); iel++) {
			for(int i = 0; i < R[iel].rows(); i++)
				for(int j = 0; j < R[iel].cols(); j++)
					R[iel](i,j) = 0.0;
		}

		spatial->update_residual(u, R, tsl);

		// step
		for(int iel = 0; iel < m->gnelem(); iel++)
		{
			u[iel] = u[iel] - cfl*tsl[iel]*R[iel]*Mi[iel];
		}

		double resnorm = spatial->computeL2Norm(R, 0);
		if(step == 0) resnorm0 = resnorm;
		else relresnorm = resnorm/resnorm0;

		step++;
		if(step % 50 == 0)
			std::cout << "  SteadyExplicit: integrate: Step " << step
				<< ", rel res = " << relresnorm << std::endl;
	}

	std::printf(" SteadyExplicit: integrate: Total steps %d, final rel res = %e\n", step, relresnorm);
}

}
