/** \file testelementtri.cpp
 * \brief Unit test for geometric mapping and Taylor basis functions on P2 triangles
 *
 * Currently, it only actually tests the locations of LG quadrature points.
 * 
 * \author Aditya Kashi
 * \date 2017 March 20
 */

#undef NDEBUG

#include <iostream>
#include "fem/aquadrature.hpp"
#include "mesh/amesh2dh.hpp"
#include "fem/aelements.hpp"

using namespace std;
using namespace amat;
using namespace tadgens;

int main()
{
	UMesh2dh m;
	m.readGmsh2("trimesh-skew_p2.msh",2);
	m.compute_topological();

	int p = 1;
	int geomdeg = 2;
	int nintp = 3;

	Quadrature2DTriangle integ;
	integ.initialize(nintp);

	const int elem = 13-1-m.gnface();
	//Array2d<a_real> elpoints(m.gnnode(elem),NDIM);
	Matrix elpoints(m.gnnode(elem),NDIM);
	for(int i = 0; i < m.gnnode(elem); i++)
		for(int j = 0; j < NDIM; j++)
			elpoints(i,j) = m.gcoords(m.ginpoel(elem,i),j);

	LagrangeMapping2D map;
	map.setAll(geomdeg, elpoints, &integ);
	map.computeForPhysicalElement();

	// check jacobian det and quadrature point coords
	/* The coords should be
	 * (2.0/3,1.4), (0.4,1.256), (0.8, 1.672), (0.8, 1.24)
	 */
	amat::Array2d<a_real> qc(4,2);
	qc(0,0) = 2/3.0; qc(0,1) = 1.4; qc(1,0) = 0.4; qc(1,1) = 1.256; qc(2,0) = 0.8; qc(2,1) = 1.672;
	qc(3,0) = 0.8; qc(3,1) = 1.24;
	//cout << std::setprecision(15);
	printf("Phy coords of quad points and Jacobian determinant of element %d:\n  ", elem+1+m.gnface());
	for(int ig = 0; ig < map.getQuadrature()->numGauss(); ig++) {
		//cout << "(" << map.map()(ig,0) << ", " << map.map()(ig,1) << "), ";
		assert(fabs(map.map()(ig,0)-qc(ig,0)) <= 10*SMALL_NUMBER
		       || fabs(map.map()(ig,1)-qc(ig,1)) <= 10*SMALL_NUMBER);
		printf("Test passed at phy coords of domain quadrature points.\n");
		cout << map.jacDet()[ig] << ".  ";
	}
	cout << endl;

	Element** tel;
	tel = new Element*[2];
	for(int i = 0; i < 2; i++)
		tel[i] = new TaylorElement();
	printf("Elem type = %d\n", tel[1]->getType());
	printf("Sizes %lu, %lu\n", sizeof(tel[0]), sizeof(tel[1]));
	tel[1]->initialize(p, &map);

	// check basis function values
	printf("Element center, delta x, delta y and area:\n");
	(reinterpret_cast<TaylorElement*>(tel[1]))->printDetails();
	int ng = tel[1]->getGeometricMapping()->getQuadrature()->numGauss(); printf("%d\n", ng);
	printf("Basis function and gradient values of element %d:\n", elem);
	for(int ig = 0; ig < ng; ig++) {
		printf("  Point %d with coords %f,%f:\n", ig, tel[1]->getGeometricMapping()->map()(ig,0),
		       tel[1]->getGeometricMapping()->map()(ig,1));
		for(int idof = 0; idof < tel[1]->getNumDOFs(); idof++) {
			printf("  %f", tel[1]->bFunc()(ig,idof));
			printf("  (%f,%f)", tel[1]->bGrad()[ig](idof,0), tel[1]->bGrad()[ig](idof,1));
		}
		printf("\n");
	}
	printf("\n");
	for(int i = 0; i < 2; i++)
		delete tel[i];
	delete [] tel;

	// 1D
	Quadrature1D i1d;
	nintp = 2;
	i1d.initialize(nintp);
	LagrangeMapping1D map1;
	int face = m.gelemface(elem, 2);
	Matrix fpoints(3,NDIM);
	printf("Face %d (point 0: %d. point 1: %d, left element: %d, right element: %d):\n",
	       face, m.gintfac(face,2)+1, m.gintfac(face,3)+1, m.gintfac(face,0)+1+m.gnface(),
	       m.gintfac(face,1)+1+m.gnface());
	for(int i = 0; i < m.gnnofa(face); i++) {
		printf(" (");
		for(int j = 0; j < NDIM; j++){
			fpoints(i,j) = m.gcoords(m.gintfac(face,2+i),j);
			printf("%f ", fpoints(i,j));
		}
		printf("), ");
	}
	printf("\n");

	map1.setAll(geomdeg, fpoints, &i1d);
	map1.computeAll();
	printf("  Phy coords of quadrature points and normals at those points:\n");
	for(int ig = 0; ig < map1.getQuadrature()->numGauss(); ig++) {
		printf("  Coord (%f,%f), normal (%f,%f)\n", map1.map()(ig,0), map1.map()(ig,1),
		       map1.normal()[ig][0], map1.normal()[ig][1]);
	}

	return 0;
}
