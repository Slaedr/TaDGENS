/** @file aspatialadvection.cpp
 * @brief Implementatio of spatial discretization for linear advection
 * @author Aditya Kashi
 * @date 2017 April 14
 */

#include "aspatialadvection.hpp"

namespace acfd {

LinearAdvection::LinearAdvection(const UMesh2dh* mesh, const int _p_degree, const char basis, 
                                 const int inoutflag, const int extrapflag)
	: SpatialBase(mesh, _p_degree, basis), inoutflow_flag(inoutflag), 
	  extrapolation_flag(extrapflag)
{
	a.resize(NDIM,1);
	a[0] = 1.0; a[1] = 0.0;
	
	std::cout << " LinearAdvection: Velocity is (" << a(0) << ", " << a(1) << ")\n";
	amag = std::sqrt(a[0]*a[0]+a[1]*a[1]);
}

void LinearAdvection::computeBoundaryState(const int iface, const Matrix& instate, 
	Matrix& bstate)
{
	if(m->gintfacbtags(iface, 0) == inoutflow_flag)
	{
		// compute normal velocity and decide whether to extrapolate 
		// or impose specified boundary value at each quadrature point
		const std::vector<Vector>& n = map1d[iface].normal();
		const Matrix& phypoints = map1d[iface].map();
		for(size_t ig = 0; ig < n.size(); ig++)
		{
			const a_real phycoords[] = {phypoints(ig,0), phypoints(ig,1)};
			a_real bval = bcfunc(phycoords);
			if(a.dot(n[ig]) >= 0)
				bstate.row(ig) = instate.row(ig);
			else
				bstate.row(ig)(0) = bval;
		}
	}
	else
	{
		bstate = instate;
	}
}

void LinearAdvection::computeNumericalFlux(
	const a_real* const uleft, const a_real* const uright, const a_real* const n, 
	a_real* const flux)
{
	a_real adotn = a[0]*n[0]+a[1]*n[1];
	if(adotn >= 0)
		flux[0] = adotn*uleft[0];
	else
		flux[0] = adotn*uright[0];
}

void LinearAdvection::update_residual(const std::vector<Matrix>& u, std::vector<Matrix>& res,
                                      std::vector<a_real>& mets)
{
#pragma omp parallel for default(shared)
	for(a_int iface = 0; iface < m->gnbface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);
		int ng = map1d[iface].getQuadrature()->numGauss();
		const std::vector<Vector>& n = map1d[iface].normal();
		const Matrix& lbasis = faces[iface].leftBasis();

		Matrix linterps(ng,NVARS), rinterps(ng,NVARS);
		Matrix fluxes(ng,NVARS);
		
		faces[iface].interpolateAll_left(u[lelem], linterps);
		computeBoundaryState(iface, linterps, rinterps);

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightandsp = map1d[iface].getQuadrature()->weights()(ig) * map1d[iface].speed()[ig];

			computeNumericalFlux(&linterps(ig,0), &rinterps(ig,0), &n[ig](0), &fluxes(ig,0));

			for(int ivar = 0; ivar < NVARS; ivar++) {
				for(int idof = 0; idof < elems[lelem]->getNumDOFs(); idof++)
					res[lelem](ivar,idof) += fluxes(ig,ivar) * lbasis(ig,idof) * weightandsp;
			}
		}
	}
	
#pragma omp parallel for default(shared)
	for(a_int iface = m->gnbface(); iface < m->gnaface(); iface++)
	{
		const a_int lelem = m->gintfac(iface,0);
		const a_int relem = m->gintfac(iface,1);
		const int ng = map1d[iface].getQuadrature()->numGauss();
		const std::vector<Vector>& n = map1d[iface].normal();
		const Matrix& lbasis = faces[iface].leftBasis();
		const Matrix& rbasis = faces[iface].rightBasis();

		Matrix linterps(ng,NVARS), rinterps(ng,NVARS);
		Matrix fluxes(ng,NVARS);
		
		faces[iface].interpolateAll_left(u[lelem], linterps);
		faces[iface].interpolateAll_right(u[relem], rinterps);

		for(int ig = 0; ig < ng; ig++)
		{
			const a_real weightandsp = map1d[iface].getQuadrature()->weights()(ig) * map1d[iface].speed()[ig];

			computeNumericalFlux(&linterps(ig,0), &rinterps(ig,0), &n[ig](0), &fluxes(ig,0));

			for(int ivar = 0; ivar < NVARS; ivar++) {
				for(int idof = 0; idof < elems[lelem]->getNumDOFs(); idof++)
#pragma omp atomic update
					res[lelem](ivar,idof) += fluxes(ig,ivar) * lbasis(ig,idof) * weightandsp;
				for(int idof = 0; idof < elems[relem]->getNumDOFs(); idof++)
#pragma omp atomic update
					res[relem](ivar,idof) -= fluxes(ig,ivar) * rbasis(ig,idof) * weightandsp;
			}
		}
	}

#pragma omp parallel for default(shared)
	for(a_int iel = 0; iel < m->gnelem(); iel++)
	{
		if(p_degree > 0)
		{
			const int ng = map2d[iel].getQuadrature()->numGauss();
			const int ndofs = elems[iel]->getNumDOFs();
			const std::vector<Matrix>& bgrads = elems[iel]->bGrad();
			const Matrix& bas = elems[iel]->bFunc();
			const Matrix& pts = elems[iel]->getGeometricMapping()->map();

			Matrix xflux(ng, NVARS), yflux(ng, NVARS);
			elems[iel]->interpolateAll(u[iel], xflux);
			yflux = a[1]*xflux;
			xflux *= a[0];
			Matrix term = Matrix::Zero(NVARS, ndofs);

			for(int ig = 0; ig < ng; ig++)
			{
				const a_real weightjacdet = map2d[iel].jacDet()[ig] * map2d[iel].getQuadrature()->weights()(ig);

				// add flux
				for(int ivar = 0; ivar < NVARS; ivar++)
					for(int idof = 0; idof < ndofs; idof++)
						term(ivar,idof) += (xflux(ig,ivar)*bgrads[ig](idof,0)
						                    + yflux(ig,ivar)*bgrads[ig](idof,1)) * weightjacdet;

				// add source term
				const a_real ptcoords[] = {pts(ig,0), pts(ig,1)};
				for(int idof = 0; idof < ndofs; idof++)
					term(0,idof) += source_term(ptcoords,0) * bas(ig,idof) * weightjacdet;
			}

			res[iel] -= term;
		}
		
		a_real hsize = 1.0;

		for(int ifa = 0; ifa < m->gnfael(iel); ifa++) {
			const a_int iface = m->gelemface(iel,ifa);
			if(hsize > m->gedgelengthsquared(iface)) hsize = m->gedgelengthsquared(iface);
		}

		mets[iel] = std::sqrt(hsize)/amag;
	}
}

// void LinearAdvection::add_source( a_real (*const rhs)(a_real, a_real, a_real),
//                                   std::vector<Matrix>& res)
// {
// 	//#pragma omp parallel for default(shared)
// 	for(a_int iel = 0; iel < m->gnelem(); iel++)
// 	{
// 		const int ng = map2d[iel].getQuadrature()->numGauss();
// 		const int ndofs = elems[iel]->getNumDOFs();
// 		const Matrix& bas = elems[iel]->bFunc();
// 		const Matrix& pts = elems[iel]->getGeometricMapping()->map();
// 		Matrix term = Matrix::Zero(NVARS, ndofs);

// 		for(int ig = 0; ig < ng; ig++)
// 		{
// 			const a_real weightjacdet = map2d[iel].jacDet()[ig] * map2d[iel].getQuadrature()->weights()(ig);
// 			const a_real ptcoords[] = {pts(ig,0), pts(ig,1)};
// 			for(int idof = 0; idof < ndofs; idof++)
// 				term(0,idof) += source_term(ptcoords,t) * bas(ig,idof) * weightjacdet;
// 		}

// 		res[iel] -= term;
// 	}
// }

// very crude
void LinearAdvection::postprocess(const std::vector<Matrix>& u)
{
	output.resize(m->gnpoin(),1);
	output.zeros();
	std::vector<int> surelems(m->gnpoin(),0);

	for(int iel = 0; iel < m->gnelem(); iel++)
	{
		if(basis_type == 'l')
		{
			//int ndofs = elems[iel]->getNumDOFs();
			for(int ino = 0; ino < m->gnfael(iel); ino++) {
				output(m->ginpoel(iel,ino)) += u[iel](0,ino);
				surelems[m->ginpoel(iel,ino)] += 1;
			}
			if(m->gnnode(iel) > m->gnfael(iel)) {
				for(int ino = m->gnfael(iel); ino < 2*m->gnfael(iel); ino++) {
					output(m->ginpoel(iel,ino)) += (u[iel](0,ino-m->gnfael(iel)) + u[iel](0, (ino-m->gnfael(iel)+1) % m->gnfael(iel)))/2.0;
					surelems[m->ginpoel(iel,ino)] += 1;
				}
				// for interior nodes, just use average of vertices
				for(int ino = 2*m->gnfael(iel); ino < m->gnnode(iel); ino++) {
					for(int jno = 0; jno < m->gnfael(iel); jno++)
						output(m->ginpoel(iel,ino)) += u[iel](0,jno);
					output(m->ginpoel(iel,ino)) /= m->gnfael(iel);
					surelems[m->ginpoel(iel,ino)] += 1;
				}
			}
		}
		
		// for Taylor, use only average values
		else
			for(int ino = 0; ino < m->gnnode(iel); ino++) {
				output(m->ginpoel(iel,ino)) += u[iel](0,0);
				surelems[m->ginpoel(iel,ino)] += 1;
			}
	}
	for(int ip = 0; ip < m->gnpoin(); ip++)
		output(ip) /= (a_real)surelems[ip];
}

a_real LinearAdvection::exact_solution(const a_real r[NDIM], const a_real t) const
{
	return sin(2*PI*r[1]);
}

}
