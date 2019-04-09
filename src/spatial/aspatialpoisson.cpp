/** @file aspatialpoisson.cpp
 * @brief Implements spatial discretization for Poisson equation
 * @author Aditya Kashi
 * @date 2016-04-10
 */

#include "spatial/aspatialpoisson.hpp"

namespace acfd {

LaplaceSIP::LaplaceSIP(const UMesh2dh* mesh, const int _p_degree, const a_real stab)
	: SpatialBase(mesh, _p_degree, 'l')
{
	computeFEData();

	for(int iel = 0; iel < m->gnelem(); iel++)
		map2d[iel].computePhysicalCoordsOfDomainQuadraturePoints();

	for(int iface = 0; iface < m->gnaface(); iface++)
		faces[iface].computeBasisGrads();

	int ndofs = elems[0]->getNumDOFs();
	dirdofflags.resize(m->gnelem()*elems[0]->getNumDOFs(), 0);
	for(int iel = 0; iel < m->gnelem(); iel++)
	{
		for(int ino = 0; ino < m->gnnode(iel); ino++) {
			a_int pno = m->ginpoel(iel,ino);
			if(m->gflag_bpoin(pno) == 1)
				dirdofflags[iel*ndofs+ino] = 1;
		}
	}
	ndirdofs = 0;
	for(size_t i = 0; i < dirdofflags.size(); i++)
		ndirdofs += dirdofflags[i];
	printf(" LaplaceSIP: No. of Dirichlet DOFs = %d\n", ndirdofs);

	cbig = 1.0e30;
	nu=1.0;
}

void LaplaceSIP::assemble()
{
	printf(" LaplaceSIP: solve: Assembling LHS and RHS\n");

	// declare LHS in coordinate (triplet) form for assembly
	typedef Eigen::Triplet<a_real> COO;
	std::vector<COO> coo; 
	const int ndofs = elems[0]->getNumDOFs();
	bg = Vector::Zero(ntotaldofs);
	
	// domain integral and RHS
	for(int ielem = 0; ielem < m->gnelem(); ielem++)
	{
		const Matrix& basis = elems[ielem]->bFunc();
		const std::vector<Matrix>& bgrad = elems[ielem]->bGrad();
		const GeomMapping2D* gmap = elems[ielem]->getGeometricMapping();
		const int ng = gmap->getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = gmap->getQuadrature()->weights();
		const Matrix& quadp = map2d[ielem].map();

		Matrix A = Matrix::Zero(ndofs,ndofs);
		Vector bl = Vector::Zero(ndofs);

		for(int ig = 0; ig < ng; ig++)
		{
			const a_real weightAndJDet = wts(ig)*map2d[ielem].jacDet()[ig];
			for(int i = 0; i < ndofs; i++) 
			{
				const a_real coords[] = {quadp(ig,0), quadp(ig,1)};
				//bl(i) += rhs(quadp(ig,0),quadp(ig,1)) * basis(ig,i) * weightAndJDet;
				bl(i) += source_term(coords,0) * basis(ig,i) * weightAndJDet;
				for(int j = 0; j < ndofs; j++) {
					A(i,j) += nu * bgrad[ig].row(i).dot(bgrad[ig].row(j)) * weightAndJDet;
				}
			}
		}

		for(int i = 0; i < ndofs; i++)
		{
			bg(ielem*ndofs+i) = bl(i);
			for(int j = 0; j < ndofs; j++) {
				coo.push_back(COO(ielem*ndofs+i, ielem*ndofs+j, A(i,j)));
			}
		}
	}

	// face integrals
	const a_int nbf = m->gnbface();
	for(int iface = nbf; iface < m->gnaface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);
		a_int relem = m->gintfac(iface,1);

		// local matrices
		Matrix Bkk = Matrix::Zero(ndofs,ndofs),
			Bkkp = Matrix::Zero(ndofs,ndofs),
			Bkpk = Matrix::Zero(ndofs,ndofs),
			Bkpkp = Matrix::Zero(ndofs,ndofs);
		Matrix Skk = Matrix::Zero(ndofs,ndofs),
			Skkp = Matrix::Zero(ndofs,ndofs),
			Skpk = Matrix::Zero(ndofs,ndofs),
			Skpkp = Matrix::Zero(ndofs,ndofs);

		// inverse of (approx) measure of face
		const a_real hinv = 1.0/std::sqrt( std::pow(m->gcoords(m->gintfac(iface,2),0) -
		                                            m->gcoords(m->gintfac(iface,3),0),2) +
				std::pow(m->gcoords(m->gintfac(iface,2),1) - m->gcoords(m->gintfac(iface,3),1),2) );

		const int ng = map1d[iface].getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = bquad->weights();
		const std::vector<Vector>& n = map1d[iface].normal();
		const std::vector<Matrix>& lgrad = faces[iface].leftBasisGrad();
		const std::vector<Matrix>& rgrad = faces[iface].rightBasisGrad();
		const Matrix& lbas = faces[iface].leftBasis();
		const Matrix& rbas = faces[iface].rightBasis();

		for(int ig = 0; ig < ng; ig++)
		{
			const a_real weightandspeed = wts(ig) * map1d[iface].speed()[ig];
			for(int i = 0; i < ndofs; i++)
				for(int j = 0; j < ndofs; j++) {
					Bkk(i,j) +=   nu*0.5 * lgrad[ig].row(j).dot(n[ig]) * lbas(ig,i) * weightandspeed;
					Bkkp(i,j) +=  nu*0.5 * rgrad[ig].row(j).dot(n[ig]) * lbas(ig,i) * weightandspeed;
					Bkpk(i,j) +=  nu*0.5 * lgrad[ig].row(j).dot(n[ig]) * rbas(ig,i) * weightandspeed;
					Bkpkp(i,j) += nu*0.5 * rgrad[ig].row(j).dot(n[ig]) * rbas(ig,i) * weightandspeed;

					Skk(i,j) +=   eta*nu*hinv * lbas(ig,i)*lbas(ig,j) * weightandspeed;
					Skkp(i,j) +=  eta*nu*hinv * lbas(ig,i)*rbas(ig,j) * weightandspeed;
					Skpk(i,j) +=  eta*nu*hinv * rbas(ig,i)*lbas(ig,j) * weightandspeed;
					Skpkp(i,j) += eta*nu*hinv * rbas(ig,i)*rbas(ig,j) * weightandspeed;
				}
		}

		// add to global stiffness matrix
		for(int i = 0; i < ndofs; i++)
			for(int j = 0; j < ndofs; j++)
			{
				coo.push_back(COO( lelem*ndofs+i, lelem*ndofs+j,
				                   -Bkk(i,j)  +Bkpk(i,j) -Bkk(j,i)  -Bkkp(j,i) +Skk(i,j)  -Skpk(i,j) ));
				coo.push_back(COO( lelem*ndofs+i, relem*ndofs+j,
				                   -Bkkp(i,j) +Bkpkp(i,j)+Bkpk(j,i) +Bkpkp(j,i)+Skpkp(i,j)-Skkp(i,j) ));
				coo.push_back(COO( relem*ndofs+i, lelem*ndofs+j,
				                   Bkpk(i,j) -Bkk(i,j)  -Bkkp(j,i) -Bkk(j,i)  +Skk(i,j)  -Skpk(i,j) ));
				coo.push_back(COO( relem*ndofs+i, relem*ndofs+j,
				                   Bkpkp(i,j)-Bkkp(i,j) +Bkpkp(j,i)+Bkpk(j,i) +Skpkp(i,j)-Skkp(i,j) ));
			}
	}
	
	for(int iface = 0; iface < m->gnbface(); iface++)
	{
		const a_int lelem = m->gintfac(iface,0);

		Matrix Bkk = Matrix::Zero(ndofs,ndofs), Skk = Matrix::Zero(ndofs,ndofs);

		// inverse of (approx) measure of face
		const a_real hinv = 1.0/std::sqrt( std::pow(m->gcoords(m->gintfac(iface,2),0) 
					                                -m->gcoords(m->gintfac(iface,3),0),2) + 
				std::pow(m->gcoords(m->gintfac(iface,2),1) - m->gcoords(m->gintfac(iface,3),1),2) );

		const int ng = map1d[iface].getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = bquad->weights();
		const std::vector<Vector>& n = map1d[iface].normal();
		const std::vector<Matrix>& lgrad = faces[iface].leftBasisGrad();
		const Matrix& lbas = faces[iface].leftBasis();

		for(int ig = 0; ig < ng; ig++)
		{
			const a_real weightandspeed = wts(ig) * map1d[iface].speed()[ig];
			for(int i = 0; i < ndofs; i++)
				for(int j = 0; j < ndofs; j++) {
					Bkk(i,j) +=   nu*0.5 * lgrad[ig].row(j).dot(n[ig]) * lbas(ig,i) * weightandspeed;

					Skk(i,j) +=   eta*nu*hinv * lbas(ig,i)*lbas(ig,j) * weightandspeed;
				}
		}

		// add to global stiffness matrix
		for(int i = 0; i < ndofs; i++)
			for(int j = 0; j < ndofs; j++)
				coo.push_back(COO( lelem*ndofs+i, lelem*ndofs+j, -Bkk(i,j)-Bkk(j,i) +Skk(i,j) ));
	}

	Ag.resize(ntotaldofs, ntotaldofs);

	// assemble
	Ag.setFromTriplets(coo.begin(), coo.end());

	// apply Dirichlet penalties
	/*for(int i = 0; i < ntotaldofs; i++)
	{
		if(dirdofflags[i]) {
			Ag.coeffRef(i,i) *= cbig;
			bg(i) = 0;
		}
	}
	printf("\n");*/
	/*Matrix agmat = Matrix(Ag);
	std::cout << agmat;*/
}

void LaplaceSIP::solve()
{
	//std::cout << Ag;
	ug = Vector::Zero(ntotaldofs);
	
	// printf(" LaplaceSIP: solve: Removing Dirichlet rows and cols\n");
	// Eigen::SparseMatrix<a_real> Af; Af.resize(ntotaldofs-ndirdofs, ntotaldofs-ndirdofs);
	// Vector bf = Vector::Zero(ntotaldofs-ndirdofs);
	// Vector uf = Vector::Zero(ntotaldofs-ndirdofs);
	// a_int I=0, J=0;
	// for(int i = 0; i < ntotaldofs; i++)
	// {
	// 	if(!dirdofflags[i]) {
	// 		for(int j = 0; j < ntotaldofs; j++)
	// 		{
	// 			if(!dirdofflags[j]) {
	// 				Af.coeffRef(I,J) = Ag.coeff(i,j);
	// 				J++;
	// 			}
	// 		}
	// 		bf(I) = bg(i);
	// 		I++;
	// 		J = 0;
	// 	}
	// }
	
	printf(" LaplaceSIP: solve: Factoring LHS...\n");
	Eigen::SparseLU<Eigen::SparseMatrix<a_real>,Eigen::COLAMDOrdering<int>> solver;
	solver.compute(Ag);
	printf(" LaplaceSIP: solve: Solving\n");
	ug = solver.solve(bg);
	printf(" LaplaceSIP: solve: Done.\n");
	
	// I=0;
	// for(int i = 0; i < ntotaldofs; i++)
	// {
	// 	if(!dirdofflags[i]) {
	// 		// assign computed value
	// 		ug(i) = uf(I);
	// 		I++;
	// 	}
	// }
}

void LaplaceSIP::postprocess(const std::vector<Matrix>& u)
{
	output.resize(m->gnpoin(),1);
	output.zeros();
	std::vector<int> surelems(m->gnpoin(),0);
	int ndofs = elems[0]->getNumDOFs();

	for(int iel = 0; iel < m->gnelem(); iel++)
	{
		// iterate over vertices of element
		for(int ino = 0; ino < m->gnfael(iel); ino++) {
			output(m->ginpoel(iel,ino)) += ug(iel*ndofs+ino);
			surelems[m->ginpoel(iel,ino)] += 1;
		}
	}
	for(int ip = 0; ip < m->gnpoin(); ip++)
		output(ip) /= (a_real)surelems[ip];
}

void LaplaceSIP::computeErrors(a_real& __restrict__ l2error, a_real& __restrict__ siperror) const
{
	std::printf(" LaplaceSIP: computeErrors: Computing the L2 and SIP norm of the error\n");
	int ndofs = elems[0]->getNumDOFs();
	l2error = 0; siperror = 0;
	
	// domain integral
	for(int ielem = 0; ielem < m->gnelem(); ielem++)
	{
		const std::vector<Matrix>& bgrad = elems[ielem]->bGrad();
		const Matrix& bfunc = elems[ielem]->bFunc();
		const GeomMapping2D* gmap = elems[ielem]->getGeometricMapping();
		int ng = gmap->getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = gmap->getQuadrature()->weights();
		const Matrix& qp = gmap->map();

		for(int ig = 0; ig < ng; ig++)
		{
			a_real lu = 0, lux = 0, luy = 0;
			for(int j = 0; j < ndofs; j++) {
				lu += ug(ielem*ndofs+j)*bfunc(ig,j);
				lux += ug(ielem*ndofs+j)*bgrad[ig](j,0);
				luy += ug(ielem*ndofs+j)*bgrad[ig](j,1);
			}
			const a_real crds[] = {qp(ig,0), qp(ig,1)};
			l2error += std::pow(lu-exact_solution(crds,0),2) * wts(ig) * gmap->jacDet()[ig];
			siperror += ( std::pow(lux-exactgradx(crds),2)
			              + std::pow(luy-exactgrady(crds),2) ) * wts(ig) * gmap->jacDet()[ig];
		}
	}

	// face integrals
	a_int nbf = m->gnbface();
	for(int iface = nbf; iface < m->gnaface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);
		a_int relem = m->gintfac(iface,1);
		
		// inverse of (approx) measure of face
		a_real hinv = 1.0/std::sqrt( std::pow(m->gcoords(m->gintfac(iface,2),0) - m->gcoords(m->gintfac(iface,3),0),2) + 
				std::pow(m->gcoords(m->gintfac(iface,2),1) - m->gcoords(m->gintfac(iface,3),1),2) );

		int ng = map1d[iface].getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = bquad->weights();
		const Matrix& lbas = faces[iface].leftBasis();
		const Matrix& rbas = faces[iface].rightBasis();

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightandspeed = wts(ig) * map1d[iface].speed()[ig];
			a_real lu = 0;
			for(int j = 0; j < ndofs; j++) {
				lu += ug(lelem*ndofs+j)*lbas(ig,j) - ug(relem*ndofs+j)*rbas(ig,j);
			}
			siperror += hinv * lu*lu * weightandspeed;
		}
	}
	
	for(int iface = 0; iface < m->gnbface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);

		// inverse of (approx) measure of face
		a_real hinv = 1.0/std::sqrt( std::pow(m->gcoords(m->gintfac(iface,2),0)
		                                      - m->gcoords(m->gintfac(iface,3),0),2) +
				std::pow(m->gcoords(m->gintfac(iface,2),1) - m->gcoords(m->gintfac(iface,3),1),2) );

		int ng = map1d[iface].getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = bquad->weights();
		const Matrix& qp = map1d[iface].map();
		const Matrix& lbas = faces[iface].leftBasis();

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightandspeed = wts(ig) * map1d[iface].speed()[ig];
			a_real lu = 0;
			for(int j = 0; j < ndofs; j++) {
				lu += ug(lelem*ndofs+j)*lbas(ig,j);
			}

			const a_real coords[] = {qp(ig,0),qp(ig,1)};
			siperror += hinv * pow(lu-exact_solution(coords,0),2) * weightandspeed;
		}
	}

	l2error = std::sqrt(l2error); siperror = std::sqrt(siperror);
	std::printf(" LaplaceSIP: computeErrors: Done.\n");
}

a_real LaplaceSIP::source_term(const a_real r[NDIM], const a_real t) const
{
	return 2.0*PI*PI*sin(PI*r[0])*sin(PI*r[1]);
}

a_real LaplaceSIP::exact_solution(const a_real r[NDIM], const a_real t) const
{
	return sin(PI*r[0])*sin(PI*r[1]);
}

a_real LaplaceSIP::exactgradx(const a_real r[NDIM]) const {
	return PI*cos(PI*r[0])*sin(PI*r[1]);
}

a_real LaplaceSIP::exactgrady(const a_real r[NDIM]) const {
	return sin(PI*r[0])*PI*cos(PI*r[1]);
}

}	// end namespace
