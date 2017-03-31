/** \file testelementtri.cpp
 * \brief Unit test for geometric mapping and Taylor basis functions on P2 triangles
 * \author Aditya Kashi
 * \date 2017 March 20
 */

#include "../aquadrature.hpp"
#include "../amesh2dh.hpp"
#include "../aelements.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	UMesh2dh m;
	m.readGmsh2("../../testcases/unittests/trimesh-skew_p2.msh",2);
	m.compute_topological();
	
	int p = 1;
	int geomdeg = 2;
	int nintp = 3;
	
	Quadrature2DTriangle integ;
	integ.initialize(nintp);

	int elem = 13-1-m.gnface();
	Array2d<acfd_real> elpoints(m.gnnode(elem),NDIM);
	for(int i = 0; i < m.gnnode(elem); i++)
		for(int j = 0; j < NDIM; j++)
			elpoints(i,j) = m.gcoords(m.ginpoel(elem,i),j);
	
	LagrangeMapping2DTriangle map;
	map.setAll(geomdeg, elpoints, &integ);
	map.computeMappingAndJacobianDet();
		
	// check jacobian det and quad point coords
	/* The coords should be
	 * (2.0/3,1.4), (0.4,1.256), (0.8, 1.672), (0.8, 1.24)
	 */
	printf("Phy coords of quad points and Jacobian determinant of element %d:\n  ", elem);
	for(int ig = 0; ig < map.getQuadrature()->numGauss(); ig++) {
		cout << "(" << map.map()(ig,0) << ", " << map.map()(ig,1) << "), ";
		cout << map.jacDet(ig) << ".  ";
	}
	cout << endl;

	TaylorElement tel;
	tel.initialize(p, &map);

	// check basis function values
	printf("Element center, delta x, delta y and area:\n");
	tel.printDetails();
	int ng = tel.getGeometricMapping()->getQuadrature()->numGauss();
	printf("Basis function and gradient values of element %d:\n", elem);
	for(int ig = 0; ig < ng; ig++) {
		printf("  Point %d with coords %f,%f:\n", ig, tel.getGeometricMapping()->map()(ig,0), tel.getGeometricMapping()->map()(ig,1));
		for(int idof = 0; idof < tel.getNumDOFs(); idof++) {
			printf("  %f", tel.bFunc(ig)(idof));
			printf("  (%f,%f)", tel.bGrad(ig)(idof,0), tel.bGrad(ig)(idof,1));
		}
		printf("\n");
	}
	printf("\n");
	
	// 1D
	Quadrature1D i1d;
	i1d.initialize(nintp);
	LagrangeMapping1D map1;
	int face = m.gelemface(elem, 2);
	
	return 0;
}
