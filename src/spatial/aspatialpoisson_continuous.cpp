/** @file aspatialpoisson_continuous.cpp
 * @brief Implements spatial discretization for Poisson equation
 * @author Aditya Kashi
 * @date 2016-04-10
 */

#include "spatial/aspatialpoisson_continuous.hpp"

namespace acfd {

LaplaceC::LaplaceC(const UMesh2dh* mesh, const int _p_degree, const a_real stab,
			a_real(*const f)(a_real,a_real), a_real(*const exact_sol)(a_real,a_real,a_real), 
			a_real(*const exact_gradx)(a_real,a_real), a_real(*const exact_grady)(a_real,a_real))
	: SpatialBase(mesh, _p_degree, 'l'), eta(stab), rhs(f), exact(exact_sol),
	  exactgradx(exact_gradx), exactgrady(exact_grady)
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

void LaplaceC::postprocess(const std::vector<Matrix>& u)
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

}
