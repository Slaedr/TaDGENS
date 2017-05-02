/** @file aspatialadvection.cpp
 * @brief Implementatio of spatial discretization for linear advection
 * @author Aditya Kashi
 * @date 2017 April 14
 */

#include "aspatialadvection.hpp"

namespace acfd {

LinearAdvection::LinearAdvection(const UMesh2dh* mesh, const int _p_degree, const char basis, const Vector vel, const a_real b_val, const int inoutflag, const int extrapflag)
	: SpatialBase(mesh, _p_degree, basis), a(vel), bval(b_val), inoutflow_flag(inoutflag), extrapolation_flag(extrapflag)
{
	std::cout << " LinearAdvection: Velocity is (" << a(0) << ", " << a(1) << ")\n";
	if(a.rows() != NDIM)
		printf("! LinearAdvection: The advection velocity vector does not have dimension %d!\n", NDIM);
	amag = std::sqrt(a[0]*a[0]+a[1]*a[1]);

	computeFEData();
	u.resize(m->gnelem());
	res.resize(m->gnelem());
	mets.resize(m->gnelem());
	for(a_int iel = 0; iel < m->gnelem(); iel++) 
	{
		u[iel].resize(NVARS, elems[iel]->getNumDOFs());
		res[iel].resize(NVARS, elems[iel]->getNumDOFs());
		for(int i = 0; i < NVARS; i++)
			for(int j = 0; j < elems[iel]->getNumDOFs(); j++) {
				u[iel](i,j) = 1.0;
				res[iel](i,j) = 0;
			}
		
		a_real hsize = 1.0;

		for(int ifa = 0; ifa < m->gnfael(iel); ifa++) {
			a_int iface = m->gelemface(iel,ifa);
			if(hsize > m->gedgelengthsquared(iface)) hsize = m->gedgelengthsquared(iface);
		}

		mets[iel] = std::sqrt(hsize)/amag;
	}

	leftfaceterms.resize(m->gnaface());
	rightfaceterms.resize(m->gnaface());
	for(a_int iface = 0; iface < m->gnaface(); iface++)
		leftfaceterms[iface].resize(NVARS, elems[m->gintfac(iface,0)]->getNumDOFs());
	for(a_int iface = m->gnbface(); iface < m->gnaface(); iface++)
		rightfaceterms[iface].resize(NVARS, elems[m->gintfac(iface,1)]->getNumDOFs());
}

void LinearAdvection::computeBoundaryState(const int iface, const Matrix& instate, Matrix& bstate)
{
	if(m->gintfacbtags(iface, 0) == inoutflow_flag)
	{
		// compute normal velocity and decide whether to extrapolate or impose specified boundary value at each quadrature point
		const std::vector<Vector>& n = map1d[iface].normal();
		for(size_t ig = 0; ig < n.size(); ig++) {
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

void LinearAdvection::computeNumericalFlux(const a_real* const uleft, const a_real* const uright, const a_real* const n, a_real* const flux)
{
	a_real adotn = a[0]*n[0]+a[1]*n[1];
	if(adotn >= 0)
		flux[0] = adotn*uleft[0];
	else
		flux[0] = adotn*uright[0];
}

void LinearAdvection::computeFaceTerms(const std::vector<Matrix>& u)
{
	for(a_int iface = 0; iface < m->gnbface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);
		int ng = map1d[iface].getQuadrature()->numGauss();
		const std::vector<Vector>& n = map1d[iface].normal();
		const Matrix& lbasis = faces[iface].leftBasis();

		Matrix linterps(ng,NVARS), rinterps(ng,NVARS);
		Matrix fluxes(ng,NVARS);
		leftfaceterms[iface] = Matrix::Zero(NVARS, elems[lelem]->getNumDOFs());
		
		faces[iface].interpolateAll_left(u[lelem], linterps);
		computeBoundaryState(iface, linterps, rinterps);

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightandsp = map1d[iface].getQuadrature()->weights()(ig) * map1d[iface].speed()[ig];

			computeNumericalFlux(&linterps(ig,0), &rinterps(ig,0), &n[ig](0), &fluxes(ig,0));

			for(int ivar = 0; ivar < NVARS; ivar++) {
				for(int idof = 0; idof < elems[lelem]->getNumDOFs(); idof++)
					leftfaceterms[iface](ivar,idof) += fluxes(ig,ivar) * lbasis(ig,idof) * weightandsp;
			}
		}
	
		// ! P0 only !
		/*a_real length = std::pow(m->gcoords(m->gintfac(iface,2),0)-m->gcoords(m->gintfac(iface,3),0),2);
		length += std::pow(m->gcoords(m->gintfac(iface,2),1)-m->gcoords(m->gintfac(iface,3),1),2);
		length = std::sqrt(length);
		leftfaceterms[iface](0,0) = fluxes(0,0)*length;*/
	}
	
	for(a_int iface = m->gnbface(); iface < m->gnaface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);
		a_int relem = m->gintfac(iface,1);
		int ng = map1d[iface].getQuadrature()->numGauss();
		const std::vector<Vector>& n = map1d[iface].normal();
		const Matrix& lbasis = faces[iface].leftBasis();
		const Matrix& rbasis = faces[iface].rightBasis();

		Matrix linterps(ng,NVARS), rinterps(ng,NVARS);
		Matrix fluxes(ng,NVARS);
		leftfaceterms[iface] = Matrix::Zero(NVARS, elems[lelem]->getNumDOFs());
		rightfaceterms[iface] = Matrix::Zero(NVARS, elems[relem]->getNumDOFs());
		
		faces[iface].interpolateAll_left(u[lelem], linterps);
		faces[iface].interpolateAll_right(u[relem], rinterps);

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightandsp = map1d[iface].getQuadrature()->weights()(ig) * map1d[iface].speed()[ig];

			computeNumericalFlux(&linterps(ig,0), &rinterps(ig,0), &n[ig](0), &fluxes(ig,0));

			for(int ivar = 0; ivar < NVARS; ivar++) {
				for(int idof = 0; idof < elems[lelem]->getNumDOFs(); idof++)
					leftfaceterms[iface](ivar,idof) += fluxes(ig,ivar) * lbasis(ig,idof) * weightandsp;
				for(int idof = 0; idof < elems[relem]->getNumDOFs(); idof++)
					rightfaceterms[iface](ivar,idof) += fluxes(ig,ivar) * rbasis(ig,idof) * weightandsp;
			}
		}
	
		// ! P0 only !
		/*a_real length = std::pow(m->gcoords(m->gintfac(iface,2),0)-m->gcoords(m->gintfac(iface,3),0),2);
		length += std::pow(m->gcoords(m->gintfac(iface,2),1)-m->gcoords(m->gintfac(iface,3),1),2);
		length = std::sqrt(length);
		leftfaceterms[iface](0,0) = fluxes(0,0)*length;
		rightfaceterms[iface](0,0) = fluxes(0,0)*length;*/
	}
}

void LinearAdvection::update_residual(const std::vector<Matrix>& u)
{
	computeFaceTerms(u);

	for(a_int iel = 0; iel < m->gnelem(); iel++)
	{
		if(p_degree > 0) {	
			int ng = map2d[iel].getQuadrature()->numGauss();
			int ndofs = elems[iel]->getNumDOFs();
			const std::vector<Matrix>& bgrads = elems[iel]->bGrad();

			Matrix xflux(ng, NVARS), yflux(ng, NVARS);
			elems[iel]->interpolateAll(u[iel], xflux);
			yflux = a[1]*xflux;
			xflux *= a[0];
			Matrix term = Matrix::Zero(NVARS, ndofs);

			for(int ig = 0; ig < ng; ig++)
			{
				a_real weightjacdet = map2d[iel].jacDet()[ig] * map2d[iel].getQuadrature()->weights()(ig);
				for(int ivar = 0; ivar < NVARS; ivar++)
					for(int idof = 0; idof < ndofs; idof++)
						term(ivar,idof) += (xflux(ig,ivar)*bgrads[ig](idof,0) + yflux(ig,ivar)*bgrads[ig](idof,1)) * weightjacdet;
			}

			res[iel] -= term;
		}

		for(int ifa = 0; ifa < m->gnfael(iel); ifa++) {
			a_int iface = m->gelemface(iel,ifa);
			a_int nbdelem = m->gesuel(iel,ifa);
			if(iel < nbdelem) {
				res[iel] += leftfaceterms[iface];
			}
			else {
				res[iel] -= rightfaceterms[iface];
			}
		}
	}
}

void LinearAdvection::add_source( a_real (*const rhs)(a_real, a_real, a_real), a_real t)
{
	for(a_int iel = 0; iel < m->gnelem(); iel++)
	{
		int ng = map2d[iel].getQuadrature()->numGauss();
		int ndofs = elems[iel]->getNumDOFs();
		const Matrix& bas = elems[iel]->bFunc();
		const Matrix& pts = elems[iel]->getGeometricMapping()->map();
		Matrix term = Matrix::Zero(NVARS, ndofs);

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightjacdet = map2d[iel].jacDet()[ig] * map2d[iel].getQuadrature()->weights()(ig);
			for(int idof = 0; idof < ndofs; idof++)
				term(0,idof) += rhs(pts(ig,0),pts(ig,1),t) * bas(ig,idof) * weightjacdet;
		}

		res[iel] -= term;
	}
}

a_real LinearAdvection::computeL2Error(double (*const exact)(double,double,double), const double time) const
{
	double l2error = 0;
	for(int iel = 0; iel < m->gnelem(); iel++)
	{
		l2error += computeElemL2Error2(iel, 0, u[iel], exact, time);
	}
	return sqrt(l2error);
}

// very crude
void LinearAdvection::postprocess()
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

}
