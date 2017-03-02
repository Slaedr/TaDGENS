#include "../mesh2dh.hpp"

using namespace amat;
using namespace acfd;
using namespace std;

int loworder_hybridmesh()
{
	UMesh2dh m;
	m.readGmsh2("../../testcases/unittests/squarequad2.msh",2);
	m.compute_topological();
	
	for(int i = 0; i < m.gnaface(); i++)
	{
	    cout << "nnofa = " << m.gnnofa(i) << "; ";
	    for(int j = 0; j < 2+m.gnnofa(i); j++)
	        cout << m.gintfac(i,j) << " ";
    	cout << endl;
    }
	
	int nbface = 8;
	
	cout << "Nbface = " << m.gnbface() << endl;
	int elem = 1;
	for(int i = 0; i < m.gnfael(elem); i++)
	    cout << m.gesuel(elem,i) << " ";
	cout << endl;
	
	/*int face = 10;
	cout << "L and R elems: " << m.gintfac(face,0) << " " << m.gintfac(face,1) << endl << "Nodes: ";
	for(int i = 0; i < m.gnnofa(face); i++)
	    cout << m.gintfac(face,2+i) << " ";
	cout << endl;*/
	
	cout << endl;
	return 0;
}

int main()
{
    loworder_hybridmesh();
    return 0;
}
