/** @file poisson.cpp
 * @brief Main function for continuous Galerkin Poisson solver
 * @author Aditya Kashi
 * @date 2017 April 11
 */

#undef NDEBUG
#include "aspatialpoisson_continuous.hpp"
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

	string dum, meshprefix, outf;
	int degree, nmesh, dirichlet_id, neumann_id;

	control >> dum; control >> nmesh;
	control >> dum; control >> meshprefix;
	control >> dum; control >> outf;
	control >> dum; control >> degree;
	control >> dum; control >> dirichlet_id;
	control >> dum; control >> neumann_id;
	control.close();

	vector<string> mfiles(nmesh), sfiles(nmesh);
	vector<double> h(nmesh,0), l2err(nmesh,0), siperr(nmesh,0);
	string names[] = {"poisson"};
	std::vector<Matrix> udum;

	for(int i = 0; i < nmesh; i++) {
		mfiles[i] = meshprefix + to_string(i) + ".msh";
		sfiles[i] = meshprefix + to_string(i) + "-p" + to_string(degree) + "-C.vtu";
	}

	std::vector<double> l2slopes(nmesh-1);

	for(int imesh = 0; imesh < nmesh; imesh++)
	{
		const UMesh2dh m = prepare_mesh(mfiles[imesh]);
		LaplaceC sd(&m, degree, dirichlet_id, neumann_id);
		sd.solve();
		sd.postprocess(udum);
		sd.computeErrors(l2err[imesh], siperr[imesh]);
		
		l2err[imesh] = log10(l2err[imesh]); siperr[imesh] = log10(siperr[imesh]);
		h[imesh] = log10(sqrt(1.0/m.gnelem()));
		printf("Mesh %d: Log mesh size = %f, log L2 error = %f, log H1 error = %f\n",
		       imesh, h[imesh], l2err[imesh], siperr[imesh]);

		const Array2d<a_real>& u = sd.getOutput();
		Array2d<a_real> vecs;
		writeScalarsVectorToVtu_PointData(sfiles[imesh], m, u, names, vecs, "none");

		if(imesh > 0) {
			double l2slope = (l2err[imesh]-l2err[imesh-1])/(h[imesh]-h[imesh-1]);
			double sipslope = (siperr[imesh]-siperr[imesh-1])/(h[imesh]-h[imesh-1]);
			printf("Mesh %d: L2 slope = %f, H1 slope = %f\n", imesh, l2slope, sipslope);
			l2slopes[imesh-1] = l2slope;
		}

		printf("\n"); fflush(stdout);
	}

	assert(std::abs(l2slopes[nmesh-2] - (degree+1.0)) < 0.15);

	ofstream convf(outf);
	for(int i = 0; i < nmesh; i++)
		convf << h[i] << " " << l2err[i] << " " << siperr[i] << "\n";
	convf.close();

	printf("---\n\n");
	return 0;
}
