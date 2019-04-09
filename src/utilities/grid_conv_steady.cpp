/** @file ladvection.cpp
 * @brief Main function for DG linear advection solver
 * 
 * Note that convergence is plotted w.r.t. 1/sqrt(num DOFs).
 * 
 * @author Aditya Kashi
 * @date 2017 April 18
 */

#undef NDEBUG

#include "spatial/aspatialadvection.hpp"
#include "solvers/atimesteady.hpp"
#include "spatial/aoutput.hpp"

using namespace amat;
using namespace std;
using namespace acfd;

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printf("Please give a control file name.\n");
		return -1;
	}

	// Read control file
	ifstream control(argv[1]);

	string dum, meshprefix, outf, outerr;
	double cfl, tol;
	int sdegree, maxits, nmesh, extrapflag, inoutflag;
	char basistype;

	control >> dum; control >> nmesh;
	control >> dum; control >> meshprefix;
	control >> dum; control >> outf;
	control >> dum; control >> basistype;
	control >> dum; control >> sdegree;
	control >> dum; control >> cfl;
	control >> dum; control >> tol;
	control >> dum; control >> maxits;
	control >> dum; control >> inoutflag;
	control >> dum; control >> extrapflag;
	control.close();

	vector<string> mfiles(nmesh), sfiles(nmesh), exfiles(nmesh);
	vector<double> h(nmesh,0), l2err(nmesh,0);
	string names[] = {"passive-scalar"};

	for(int i = 0; i < nmesh; i++) {
		mfiles[i] = meshprefix + to_string(i) + ".msh";
		sfiles[i] = outf + to_string(i) + "-p"+to_string(sdegree)+".vtu";
		exfiles[i] = outf + to_string(i) + "-exact.vtu";
	}
	outerr = outf + "-p" + to_string(sdegree)+".txt";

	std::vector<double> l2slopes(nmesh-1);

	for(int imesh = 0; imesh < nmesh; imesh++)
	{
		UMesh2dh m;
		m.readGmsh2(mfiles[imesh], NDIM);
		m.compute_topological();
		m.compute_boundary_maps();
		
		const double hhactual = m.meshSizeParameter();
		printf("Mesh %d: h = %f\n", imesh, hhactual);

		LinearAdvection sd(&m, sdegree, basistype, inoutflag, extrapflag);
		const double hh = 1.0/sqrt(sd.numTotalDOFs());
		
		SteadyExplicit<1> td(&m, &sd, cfl, tol, maxits);
		
		td.integrate();

		sd.postprocess(td.solution());
		l2err[imesh] = sd.computeL2Error(0, td.solution());
		
		l2err[imesh] = log10(l2err[imesh]);
		h[imesh] = log10(hh);
		printf("Mesh %d: Log mesh size = %f, log L2 error = %f\n", imesh, h[imesh], l2err[imesh]);

		const Array2d<a_real>& u = sd.getOutput();
		Array2d<a_real> vecs;
		writeScalarsVectorToVtu_PointData(sfiles[imesh], m, u, names, vecs, "none");

		/*Array2d<a_real> exactout(m.gnpoin(),1);
		for(a_int i = 0; i < m.gnpoin(); i++)
			exactout(i) = exactsol(m.gcoords(i,0),m.gcoords(i,1),0);
		writeScalarsVectorToVtu_PointData(exfiles[imesh], m, exactout, names, vecs, "none");*/

		if(imesh > 0) {
			const double l2slope = (l2err[imesh]-l2err[imesh-1])/(h[imesh]-h[imesh-1]);
			printf("L2 error slope (%d) = %f\n", imesh, l2slope);
			l2slopes[imesh-1] = l2slope;
		}
		cout << endl;
	}

	/*ofstream convf(outerr);
	for(int i = 0; i < nmesh; i++)
		convf << h[i] << " " << l2err[i] << "\n";
	convf.close();*/

	/// \todo FIXME: Use problem-dependent tolerance
	assert(std::abs(l2slopes[nmesh-2] - (sdegree+1.0)) < 0.15);

	printf("---\n\n");
	return 0;
}
