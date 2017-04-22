/** @file aspatialpoisson.cpp
 * @brief Implements spatial discretization for Poisson equation
 * @author Aditya Kashi
 * @date 2016-04-10
 */

#include "aspatialpoisson.hpp"

namespace acfd {

LaplaceSIP::LaplaceSIP(const UMesh2dh* mesh, const int _p_degree, const a_real stab,
			a_real(*const f)(a_real,a_real), a_real(*const exact_sol)(a_real,a_real,a_real), 
			a_real(*const exact_gradx)(a_real,a_real), a_real(*const exact_grady)(a_real,a_real))
	: SpatialBase(mesh, _p_degree, 'l'), eta(stab), rhs(f), exact(exact_sol), exactgradx(exact_gradx), exactgrady(exact_grady)
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

	Ag.resize(ntotaldofs, ntotaldofs);
	bg = Vector::Zero(ntotaldofs);
	ug = Vector::Zero(ntotaldofs);

	cbig = 1.0e30;
	nu=1.0;
}

void LaplaceSIP::assemble()
{
	// declare LHS in coordinate (triplet) form for assembly
	typedef Eigen::Triplet<a_real> COO;
	std::vector<COO> coo; 
	int ndofs = elems[0]->getNumDOFs();
	
	// domain integral and RHS
	for(int ielem = 0; ielem < m->gnelem(); ielem++)
	{
		const Matrix& basis = elems[ielem]->bFunc();
		const std::vector<Matrix>& bgrad = elems[ielem]->bGrad();
		const GeomMapping2D* gmap = elems[ielem]->getGeometricMapping();
		int ng = gmap->getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = gmap->getQuadrature()->weights();
		const Matrix& quadp = map2d[ielem].map();

		Matrix A = Matrix::Zero(ndofs,ndofs);
		Vector bl = Vector::Zero(ndofs);

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightAndJDet = wts(ig)*map2d[ielem].jacDet()[ig];
			for(int i = 0; i < ndofs; i++) 
			{
				bl(i) += rhs(quadp(ig,0),quadp(ig,1)) * basis(ig,i) * weightAndJDet;
				for(int j = 0; j < ndofs; j++) {
					A(i,j) += nu * bgrad[ig].row(i).dot(bgrad[ig].row(j)) * weightAndJDet;
				}
			}
		}
		//std::cout << A << std::endl << std::endl;

		for(int i = 0; i < ndofs; i++)
		{
			bg(ielem*ndofs+i) = bl(i);
			for(int j = 0; j < ndofs; j++) {
					coo.push_back(COO(ielem*ndofs+i, ielem*ndofs+j, A(i,j)));
			}
		}
	}

	// face integrals
	a_int nbf = m->gnbface();
	for(int iface = nbf; iface < m->gnaface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);
		a_int relem = m->gintfac(iface,1);

		// local matrices
		Matrix Bkk = Matrix::Zero(ndofs,ndofs), Bkkp = Matrix::Zero(ndofs,ndofs), Bkpk = Matrix::Zero(ndofs,ndofs), Bkpkp = Matrix::Zero(ndofs,ndofs);
		Matrix Skk = Matrix::Zero(ndofs,ndofs), Skkp = Matrix::Zero(ndofs,ndofs), Skpk = Matrix::Zero(ndofs,ndofs), Skpkp = Matrix::Zero(ndofs,ndofs);

		// inverse of (approx) measure of face
		a_real hinv = 1.0/std::sqrt( std::pow(m->gcoords(m->gintfac(iface,2),0) - m->gcoords(m->gintfac(iface,3),0),2) + 
				std::pow(m->gcoords(m->gintfac(iface,2),1) - m->gcoords(m->gintfac(iface,3),1),2) );

		int ng = map1d[iface].getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = bquad->weights();
		const std::vector<Vector>& n = map1d[iface].normal();
		const std::vector<Matrix>& lgrad = faces[iface].leftBasisGrad();
		const std::vector<Matrix>& rgrad = faces[iface].rightBasisGrad();
		const Matrix& lbas = faces[iface].leftBasis();
		const Matrix& rbas = faces[iface].rightBasis();

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightandspeed = wts(ig) * map1d[iface].speed()[ig];
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
				coo.push_back(COO( lelem*ndofs+i, lelem*ndofs+j, -Bkk(i,j)  +Bkpk(i,j) -Bkk(j,i)  -Bkkp(j,i) +Skk(i,j)  -Skpk(i,j) ));
				coo.push_back(COO( lelem*ndofs+i, relem*ndofs+j, -Bkkp(i,j) +Bkpkp(i,j)+Bkpk(j,i) +Bkpkp(j,i)+Skpkp(i,j)-Skkp(i,j) ));
				coo.push_back(COO( relem*ndofs+i, lelem*ndofs+j,  Bkpk(i,j) -Bkk(i,j)  -Bkkp(j,i) -Bkk(j,i)  +Skk(i,j)  -Skpk(i,j) ));
				coo.push_back(COO( relem*ndofs+i, relem*ndofs+j,  Bkpkp(i,j)-Bkkp(i,j) +Bkpkp(j,i)+Bkpk(j,i) +Skpkp(i,j)-Skkp(i,j) ));
			}
	}
	
	for(int iface = 0; iface < m->gnbface(); iface++)
	{
		a_int lelem = m->gintfac(iface,0);

		Matrix Bkk = Matrix::Zero(ndofs,ndofs), Skk = Matrix::Zero(ndofs,ndofs);

		// inverse of (approx) measure of face
		a_real hinv = 1.0/std::sqrt( std::pow(m->gcoords(m->gintfac(iface,2),0) - m->gcoords(m->gintfac(iface,3),0),2) + 
				std::pow(m->gcoords(m->gintfac(iface,2),1) - m->gcoords(m->gintfac(iface,3),1),2) );

		int ng = map1d[iface].getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = bquad->weights();
		const std::vector<Vector>& n = map1d[iface].normal();
		const std::vector<Matrix>& lgrad = faces[iface].leftBasisGrad();
		const Matrix& lbas = faces[iface].leftBasis();

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightandspeed = wts(ig) * map1d[iface].speed()[ig];
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
	printf(" LaplaceSIP: solve: Assembling LHS and RHS\n");
	assemble();
	//std::cout << Ag;
	
	printf(" LaplaceSIP: solve: Removing Dirichlet rows and cols\n");
	Eigen::SparseMatrix<a_real> Af; Af.resize(ntotaldofs-ndirdofs, ntotaldofs-ndirdofs);
	Vector bf = Vector::Zero(ntotaldofs-ndirdofs);
	Vector uf = Vector::Zero(ntotaldofs-ndirdofs);
	a_int I=0, J=0;
	for(int i = 0; i < ntotaldofs; i++)
	{
		if(!dirdofflags[i]) {
			for(int j = 0; j < ntotaldofs; j++)
			{
				if(!dirdofflags[j]) {
					Af.coeffRef(I,J) = Ag.coeff(i,j);
					J++;
				}
			}
			bf(I) = bg(i);
			I++;
			J = 0;
		}
	}
	
	printf(" LaplaceSIP: solve: Analyzing and factoring LHS...\n");
	Eigen::SparseLU<Eigen::SparseMatrix<a_real>> solver;
	solver.compute(Af);
	printf(" LaplaceSIP: solve: Solving\n");
	uf = solver.solve(bf);
	printf(" LaplaceSIP: solve: Done.\n");
	
	I=0;
	//int ndofs = elems[0].getNumDOFs();
	for(int i = 0; i < ntotaldofs; i++)
	{
		if(!dirdofflags[i]) {
			// assign computed value
			ug(i) = uf(I);
			I++;
		}
	}
}

void LaplaceSIP::postprocess()
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
			l2error += std::pow(lu-exact(qp(ig,0),qp(ig,1),0),2) * wts(ig) * gmap->jacDet()[ig];
			siperror += ( std::pow(lux-exactgradx(qp(ig,0),qp(ig,1)),2) + std::pow(luy-exactgrady(qp(ig,0),qp(ig,1)),2) ) * wts(ig) * gmap->jacDet()[ig];
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
		a_real hinv = 1.0/std::sqrt( std::pow(m->gcoords(m->gintfac(iface,2),0) - m->gcoords(m->gintfac(iface,3),0),2) + 
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
			siperror += hinv * pow(lu-exact(qp(ig,0),qp(ig,1),0),2) * weightandspeed;
		}
	}

	l2error = std::sqrt(l2error); siperror = std::sqrt(siperror);
	std::printf(" LaplaceSIP: computeErrors: Done.\n");
}

LaplaceC::LaplaceC(const UMesh2dh* mesh, const int _p_degree, const a_real stab,
			a_real(*const f)(a_real,a_real), a_real(*const exact_sol)(a_real,a_real,a_real), 
			a_real(*const exact_gradx)(a_real,a_real), a_real(*const exact_grady)(a_real,a_real))
	: SpatialBase(mesh, _p_degree, 'l'), eta(stab), rhs(f), exact(exact_sol), exactgradx(exact_gradx), exactgrady(exact_grady)
{
	computeFEData();

	for(int iel = 0; iel < m->gnelem(); iel++)
		map2d[iel].computePhysicalCoordsOfDomainQuadraturePoints();

	for(int iface = 0; iface < m->gnaface(); iface++)
		faces[iface].computeBasisGrads();

	ntotaldofs = m->gnpoin() + (p_degree-1)*m->gnaface();
	std::cout << " LaplaceC: Total DOFs = " << ntotaldofs << std::endl;
	Ag.resize(ntotaldofs, ntotaldofs);
	bg = Vector::Zero(ntotaldofs);
	ug = Vector::Zero(ntotaldofs);

	cbig = 1.0e30;
	nu=1.0;

	int nlocdofs = elems[0]->getNumDOFs();
	std::cout << " LaplaceC: Local DOFs = " << nlocdofs << std::endl;
	dofmap.resize(m->gnelem(), nlocdofs);
	bflag = Vector::Zero(ntotaldofs);
	for(int i = 0; i < m->gnpoin(); i++)
		bflag(i) = m->gflag_bpoin(i);

	for(int i = 0; i < m->gnelem(); i++)
	{
		for(int j = 0; j < nlocdofs; j++) 
		{
			if(j < m->gnnode(i)) {
				dofmap(i,j) = m->ginpoel(i,j);
			}
			else if(j < m->gnnode(i)*2){
				a_int face = m->gelemface(i, j-m->gnnode(i));
				dofmap(i,j) = m->gnpoin()+face;
				if(face < m->gnbface()) bflag(dofmap(i,j)) = 1;
			}
			else {
				std::cout << "! LaplaceC: Dofmap not implemented for this kind of element!\n";
			}
			//std::cout << "  Dofmap " << i << "," << j << " = " << dofmap(i,j) << "\n";
		}
	}
}

void LaplaceC::assemble()
{
	// declare LHS in coordinate (triplet) form for assembly
	typedef Eigen::Triplet<a_real> COO;
	std::vector<COO> coo; 
	int ndofs = elems[0]->getNumDOFs();
	
	// domain integral
	for(int ielem = 0; ielem < m->gnelem(); ielem++)
	{
		const Matrix& basis = elems[ielem]->bFunc();
		const std::vector<Matrix>& bgrad = elems[ielem]->bGrad();
		const GeomMapping2D* gmap = elems[ielem]->getGeometricMapping();
		int ng = gmap->getQuadrature()->numGauss();
		const amat::Array2d<a_real>& wts = gmap->getQuadrature()->weights();
		const Matrix& quadp = map2d[ielem].map();

		Matrix A = Matrix::Zero(ndofs,ndofs);
		Vector bl = Vector::Zero(ndofs);

		for(int ig = 0; ig < ng; ig++)
		{
			a_real weightAndJDet = wts(ig)*map2d[ielem].jacDet()[ig];
			for(int i = 0; i < ndofs; i++) 
			{
				bl(i) += rhs(quadp(ig,0),quadp(ig,1)) * basis(ig,i) * weightAndJDet;
				for(int j = 0; j < ndofs; j++) {
					A(i,j) += nu * bgrad[ig].row(i).dot(bgrad[ig].row(j)) * weightAndJDet;
				}
			}
		}
		//std::cout << A << std::endl << std::endl;
		
		for(int i = 0; i < ndofs; i++)
		{
			bg(dofmap(ielem,i)) += bl(i);
			for(int j = 0; j < ndofs; j++) {
					coo.push_back(COO(dofmap(ielem,i), dofmap(ielem,j), A(i,j)));
			}
		}
	}
	
	// assemble
	Ag.setFromTriplets(coo.begin(), coo.end());

	// apply Dirichlet penalties
	for(int i = 0; i < ntotaldofs; i++)
	{
		if(bflag(i)) {
			Ag.coeffRef(i,i) *= cbig;
			bg(i) = 0;
		}
	}
}

void LaplaceC::solve()
{
	printf(" LaplaceC: solve: Assembling LHS and RHS\n");
	assemble();
	
	printf(" LaplaceC: solve: Analyzing and factoring LHS...\n");
	Eigen::SparseLU<Eigen::SparseMatrix<a_real>> solver;
	solver.compute(Ag);
	printf(" LaplaceC: solve: Solving\n");
	ug = solver.solve(bg);
	printf(" LaplaceC: solve: Done.\n");
}

void LaplaceC::postprocess()
{
	output.resize(m->gnpoin());
	for(int i = 0; i < m->gnpoin(); i++)
		output(i) = ug(i);
}

void LaplaceC::computeErrors(a_real& __restrict__ l2error, a_real& __restrict__ siperror) const
{
	std::printf(" LaplaceC: computeErrors: Computing the L2 and H1 norm of the error\n");
	int ndofs = elems[0]->getNumDOFs();
	l2error = 0; siperror = 0;
	
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
				lu += ug(dofmap(ielem,j))*bfunc(ig,j);
				lux += ug(dofmap(ielem,j))*bgrad[ig](j,0);
				luy += ug(dofmap(ielem,j))*bgrad[ig](j,1);
			}
			l2error += std::pow(lu-exact(qp(ig,0),qp(ig,1),0),2) * wts(ig) * gmap->jacDet()[ig];
			siperror += ( std::pow(lux-exactgradx(qp(ig,0),qp(ig,1)),2) + std::pow(luy-exactgrady(qp(ig,0),qp(ig,1)),2) ) * wts(ig) * gmap->jacDet()[ig];
		}
	}

	siperror += l2error;
	l2error = std::sqrt(l2error); siperror = std::sqrt(siperror);
	std::printf(" LaplaceC: computeErrors: Done.\n");
}

}	// end namespace
