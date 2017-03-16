#include "../aquadrature.hpp"
#include "../amesh2dh.hpp"
#include "../aelement.hpp"

using namespace std, amat, acfd;

int main()
{
	UMesh2dh m;
	m.readGmsh2("../../testcases/unittests/trimesh.msh",2);
	
	int p = 1;
	int nintp = 2;
	
	Quadrature2DTriangle integ;
	integ.initialize(nintp);
	
	LagrangeMapping2DTriangle map;
	
	return 0;
}
