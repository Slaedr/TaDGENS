/** @file atimesteady.hpp
 * @brief Implementation of steady-state solvers
 * @author Aditya Kashi
 * @date 2017 April 15
 */

#include "atimesteady.hpp"

namespace acfd {

SteadyBase::SteadyBase(const UMesh2dh *const mesh, SpatialBase* s, a_real cflnumber, double toler, int max_iter, bool use_source)
	: m(mesh), spatial(s), cfl(cflnumber), tol(toler), maxiter(max_iter), source(use_source)
{
	std::cout << " SteadyBase: CFL = " << cfl << ", use source? " << source << std::endl;
}

SteadyExplicit::SteadyExplicit(const UMesh2dh*const mesh, SpatialBase* s, a_real cflnumber, double toler, int max_iter, bool use_source)
	: SteadyBase(mesh, s, cflnumber, toler, max_iter, use_source)
{
}

void SteadyExplicit::integrate()
{
	int step = 0;
	double relresnorm = 1.0, resnorm0 = 1.0;
	std::vector<Matrix>& u = spatial->unk();

	std::vector<Matrix>& R = spatial->residual();
	const std::vector<Matrix>& Mi = spatial->massInv();
	const std::vector<a_real>& tsl = spatial->maxExplicitTimeStep();

	while((relresnorm > tol && step < maxiter))
	{
		for(int iel = 0; iel < m->gnelem(); iel++) {
			R[iel] = Matrix::Zero(R[iel].rows(), R[iel].cols());
		}

		spatial->update_residual(u);
		if(source)
			spatial->add_source(rhs,0);

		// step
		for(int iel = 0; iel < m->gnelem(); iel++)
		{
			u[iel] = u[iel] - cfl*tsl[iel]*R[iel]*Mi[iel];
		}

		double resnorm = spatial->computeL2Norm(R, 0);
		if(step == 0) resnorm0 = resnorm;
		else relresnorm = resnorm/resnorm0;

		step++;
		if(step % 20 == 0)
			std::printf("  SteadyExplicit: integrate: Step %d, rel res = %f\n", step, relresnorm);
	}
	
	std::printf(" SteadyExplicit: integrate: Total steps %d, final rel res = %f\n", step, relresnorm);
}

}
