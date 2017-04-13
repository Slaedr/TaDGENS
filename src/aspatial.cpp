/** @file aspatial.cpp
 * @brief Implements spatial discretizations.
 * @author Aditya Kashi
 * @date 2016-02-24
 */

#include "aspatial.hpp"

namespace acfd {

SpatialBase::SpatialBase(const UMesh2dh* mesh, const int _p_degree, char basistype) : m(mesh), p_degree(_p_degree)
{
	std::cout << " SpatialBase: Setting up spaital integrator for FE polynomial degree " << p_degree << std::endl;
	
	// set quadrature strength
	int dom_quaddegree = 4*p_degree;
	if(m->degree() == 2) dom_quaddegree += 1;
	int boun_quaddegree = 4*p_degree;
	if(m->degree() == 2) boun_quaddegree += 1;
	std::cout << " SpatialBase: Quadrature strengths for domain and boundary integrals set at " << dom_quaddegree << ", " << boun_quaddegree << std::endl;

	dtquad = new Quadrature2DTriangle();
	dtquad->initialize(dom_quaddegree);
	dsquad = new Quadrature2DSquare();
	dsquad->initialize(dom_quaddegree);
	bquad = new Quadrature1D();
	bquad->initialize(boun_quaddegree);

	map2d = new LagrangeMapping2D[m->gnelem()];
	if(basistype == 't')
		elems = new TaylorElement[m->gnelem()];
	else
		elems = new LagrangeElement[m->gnelem()];
	dummyelem = new DummyElement();

	map1d = new LagrangeMapping1D[m->gnaface()];
	faces = new FaceElement[m->gnaface()];
}

SpatialBase::~SpatialBase()
{
	delete dtquad;
	delete dsquad;
	delete bquad;
	delete [] map2d;
	delete [] map1d;
	delete [] faces;
	delete [] elems;
	delete dummyelem;
}

void SpatialBase::computeFEData()
{
	std::cout << " SpatialBase: computeFEData(): Computing basis functions, basis gradients and mass matrices for each element" << std::endl;
	minv.resize(m->gnelem());
	ntotaldofs = 0;

	// loop over elements to setup maps and elements and compute mass matrices
	for(int iel = 0; iel < m->gnelem(); iel++)
	{
		amat::Array2d<a_real> phynodes(m->gnnode(iel),NDIM);
		for(int i = 0; i < m->gnnode(iel); i++)
			for(int j = 0; j < NDIM; j++)
				phynodes(i,j) = m->gcoords(m->ginpoel(iel,i),j);

		if(m->gnnode(iel) == 4 || m->gnnode(iel) == 9)
			map2d[iel].setAll(m->degree(), phynodes, dsquad);
		else
			map2d[iel].setAll(m->degree(), phynodes, dtquad);

		elems[iel].initialize(p_degree, &map2d[iel]);
		ntotaldofs += elems[iel].getNumDOFs();

		// allocate mass matrix
		minv[iel] = Matrix::Zero(elems[iel].getNumDOFs(), elems[iel].getNumDOFs());

		// compute mass matrix
		const Quadrature2D* lquad = map2d[iel].getQuadrature();
		for(int ig = 0; ig < lquad->numGauss(); ig++)
		{
			for(int idof = 0; idof < elems[iel].getNumDOFs(); idof++)
				for(int jdof = 0; jdof < elems[iel].getNumDOFs(); jdof++)
					minv[iel](idof,jdof) += elems[iel].bFunc()(ig,idof)*elems[iel].bFunc()(ig,jdof)*map2d[iel].jacDet()[ig];
		}
	}
	std::printf(" SpatialBase: computeFEData: Total number of DOFs = %d\n", ntotaldofs);

	dummyelem->initialize(p_degree, &map2d[0]);

	// loop over faces
	for(int iface = 0; iface < m->gnbface(); iface++)
	{
		int lelem = m->gintfac(iface,0);
		amat::Array2d<a_real> phynodes(m->gnnofa(iface),NDIM);
		for(int i = 0; i < m->gnnofa(iface); i++)
			for(int j = 0; j < NDIM; j++)
				phynodes(i,j) = m->gcoords(m->gintfac(iface,2+i),j);

		map1d[iface].setAll(m->degree(), phynodes, bquad);
		map1d[iface].computeAll();

		faces[iface].initialize(&elems[lelem], dummyelem, &map1d[iface], m->gfacelocalnum(iface,0), m->gfacelocalnum(iface,1));
	}

	for(int iface = m->gnbface(); iface < m->gnaface(); iface++)
	{
		int lelem = m->gintfac(iface,0);
		int relem = m->gintfac(iface,1);
		amat::Array2d<a_real> phynodes(m->gnnofa(iface),NDIM);
		for(int i = 0; i < m->gnnofa(iface); i++)
			for(int j = 0; j < NDIM; j++)
				phynodes(i,j) = m->gcoords(m->gintfac(iface,2+i),j);

		map1d[iface].setAll(m->degree(), phynodes, bquad);
		map1d[iface].computeAll();

		faces[iface].initialize(&elems[lelem], &elems[relem], &map1d[iface], m->gfacelocalnum(iface,0), m->gfacelocalnum(iface,1));
	}
	
	std::cout << "               Mesh degree = " << m->degree() << ", geom map degee = " << map2d[0].getDegree() 
		 << ", element degree = " << elems[0].getDegree() << std::endl;
	std::cout << " SpatialBase: computeFEData(): Done." << std::endl;
}

}	// end namespace
