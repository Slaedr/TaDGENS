#include "../aquadrature.hpp"
#include "../amesh2dh.hpp"
#include "../aelements.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	// global setup
	UMesh2dh m;
	m.readGmsh2("../../testcases/unittests/trimesh.msh",2);
	
	int p = 1;
	int geomdeg = 1;
	int nintp = 3;
	
	Quadrature2DTriangle integ;
	integ.initialize(nintp);

	// per-element setup
	
	int elem = 1;
	Array2d<acfd_real> elpoints(m.gnnode(elem),NDIM);
	for(int i = 0; i < m.gnnode(elem); i++)
		for(int j = 0; j < NDIM; j++)
			elpoints(i,j) = m.gcoords(m.ginpoel(elem,i),j);
	
	LagrangeMapping2DTriangle map;
	map.setAll(geomdeg, elpoints, &integ);
	map.computeMappingAndJacobianDet();

	TaylorElement tel;
	tel.initialize(p, &map);
		
	// check jacobian det
	printf("Jacobian determinant of element %d:\n  ", elem);
	for(int ig = 0; ig < map.getQuadrature()->numGauss(); ig++)
		cout << map.jacDet(ig) << " ";
	cout << endl;

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
	
	return 0;
}
