#include "../aquadrature.hpp"
#include "../amesh2dh.hpp"
#include "../aelements.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	UMesh2dh m;
	m.readGmsh2("../../testcases/unittests/trimesh.msh",2);
	
	int p = 1;
	int nintp = 3;
	
	Quadrature2DTriangle integ;
	integ.initialize(nintp);
	
	LagrangeMapping2DTriangle map;
	
	return 0;
}
