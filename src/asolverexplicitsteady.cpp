/** @file asolverexplicitsteady.cpp
 * @brief Implementation of explicit TVDRK solver
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#include "asolverexplicitsteady.hpp"

void SolverExplicitSteady::solve_rk1_steady()
{
	int step = 0;
	acfd_real resi = 1.0;
	acfd_real initres = 1.0;
	amat::Array2d<acfd_real> res(NVARS,1);
	res.ones();
	amat::Array2d<acfd_real> dtm(m->gnelem(), 1);		// for local time-stepping
	amat::Array2d<acfd_real> uold(u.rows(), u.cols());

	while(resi/initres > tol && step < maxiter)
	{
		//std::cout << "EulerFV: solve_rk1_steady(): Entered loop. Step " << step << std::endl;

		//calculate fluxes
		compute_RHS();		// this invokes Flux calculating function after zeroing the residuals, also computes max wave speeds integ

		acfd_real err[NVARS];
		acfd_real errmass = 0;
		for(int i = 0; i < NVARS; i++)
			err[i] = 0;

		//calculate dt based on CFL
#pragma omp parallel default(shared)
		{
#pragma omp for simd
			for(int iel = 0; iel < m->gnelem(); iel++)
			{
				dtm(iel) = cfl*(0.5*m->gjacobians(iel)/integ(iel));
			}

#pragma omp for simd
			for(int iel = 0; iel < m->gnelem(); iel++)
			{
				for(int i = 0; i < NVARS; i++)
				{
					//uold(iel,i) = u(iel,i);
					u(iel,i) += dtm.get(iel)*m_inverse.get(iel)*residual.get(iel,i);
				}
			}

#pragma omp for simd reduction(+:errmass)
//#pragma omp for simd reduction(+:err[:NVARS])
			for(int iel = 0; iel < m->gnelem(); iel++)
			{
				/*for(int i = 0; i < NVARS; i++)
				{
					err[i] += residual(iel,i)*residual(iel,i)*m->garea(iel);
				}*/
				errmass += residual(iel,0)*residual(iel,0)*m->garea(iel);
			}
		} // end parallel region

		/*resi = 2e-15;
		for(int i = 0; i < NVARS; i++)
			if(err[i] > resi*resi)
				resi = err[i];
		resi = sqrt(resi);*/
		resi = sqrt(errmass);

		if(step == 0)
			initres = resi;

		if(step % 50 == 0)
			std::cout << "EulerFV: solve_rk1_steady(): Step " << step << ", rel residual " << resi/initres << std::endl;

		step++;
		/*acfd_real totalenergy = 0;
		for(int i = 0; i < m->gnelem(); i++)
			totalenergy += u(i,3)*m->jacobians(i);
			std::cout << "EulerFV: solve(): Total energy = " << totalenergy << std::endl;*/
		//if(step == 10000) break;
	}

	if(step == maxiter)
		std::cout << "ExplicitSolver: solve_rk1_steady(): Exceeded max iterations!" << std::endl;
}
