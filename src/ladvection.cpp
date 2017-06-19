/** @file ladvection.cpp
 * @brief Main function for DG linear advection solver
 * 
 * Note that convergence is plotted w.r.t. 1/sqrt(num DOFs).
 * 
 * @author Aditya Kashi
 * @date 2017 April 18
 */

#include "aspatialadvection.hpp"
#include "atimesteady.hpp"
#include "aoutput.hpp"

using namespace amat;
using namespace std;
using namespace acfd;

double a0 = 1.0, a1 = 0.0;

// exact solution - sin in x
/*double exactsol(double x, double y, double t) {
	return 1.0+sin( 2*PI/3.0*(x+1.5));
}
// and corresponding RHS
double rhs(double x, double y, double t) {
	return 2*PI/3.0 * cos( 2*PI/3.0*(x+1.5));
}*/

double exactsol(double x, double y, double t) {
	return sin(2*PI*y);
}
// corresponding boundary distribution
a_real bcfunc(const a_real x, const a_real y)
{
	return std::sin(2*PI*y);
}

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

	for(int imesh = 0; imesh < nmesh; imesh++)
	{
		UMesh2dh m; m.readGmsh2(mfiles[imesh], NDIM); m.compute_topological(); m.compute_boundary_maps();
		
		double hh = m.meshSizeParameter();
		printf("Mesh %d: h = %f\n", imesh, hh);

		Vector a(2); a[0] = a0; a[1] = a1;
		LinearAdvection sd(&m, sdegree, basistype, a, inoutflag, extrapflag, bcfunc);
		hh = 1.0/sqrt(sd.numTotalDOFs());
		
		SteadyExplicit<1> td(&m, &sd, cfl, tol, maxits, false);
		//td.set_source(rhs);
		
		td.integrate();

		sd.postprocess(td.solution());
		l2err[imesh] = sd.computeL2Error(exactsol, 0, td.solution());
		
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
			double l2slope = (l2err[imesh]-l2err[imesh-1])/(h[imesh]-h[imesh-1]);
			printf("L2 error slope (%d) = %f\n", imesh, l2slope);
		}
		cout << endl;
	}

	ofstream convf(outerr);
	for(int i = 0; i < nmesh; i++)
		convf << h[i] << " " << l2err[i] << "\n";
	convf.close();

	printf("---\n\n");
	return 0;
}
