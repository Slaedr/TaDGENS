/** @file ladvection.cpp
 * @brief Main function for DG linear advection solver
 * @author Aditya Kashi
 * @date 2017 April 15
 */

#include "aspatialadvection.hpp"
#include "atimetvdrk.hpp"
#include "aoutput.hpp"

using namespace amat;
using namespace std;
using namespace acfd;

// exact solution - Gaussian bump
double beta = 200.0, xc = -0.2, yc = 0;
double a0 = 1.0, a1 = 0.0;

double init(double x, double y) {
	return std::exp(beta*(-(x-xc)*(x-xc)-(y-yc)*(y-yc)));
}
double initgradx(double x, double y) {
	return -init(x,y)*beta*2*(x-xc);
}
double initgrady(double x, double y) {
	return -init(x,y)*beta*2*(y-yc);
}
double initgradxx(double x, double y) {
	return 4*beta*beta*init(x,y)*(x-xc)*(x-xc) - 2*beta*init(x,y);
}
double initgradyy(double x, double y) {
	return 4*beta*beta*init(x,y)*(y-yc)*(y-yc) - 2*beta*init(x,y);
}
double initgradxy(double x, double y) {
	return 4*beta*beta*init(x,y)*(x-xc)*(y-yc);
}

double exactsol(double x, double y, double t) {
	return init(x-a0*t, y-a1*t);
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

	string dum, meshprefix, outf;
	double cfl, tstep, ftime;
	int sdegree, tdegree, nmesh, extrapflag, inoutflag;
	char btype;

	control >> dum; control >> nmesh;
	control >> dum; control >> meshprefix;
	control >> dum; control >> outf;
	control >> dum; control >> btype;
	control >> dum; control >> sdegree;
	control >> dum; control >> tdegree;
	control >> dum; control >> ftime;
	control >> dum; control >> cfl;
	control >> dum; control >> inoutflag;
	control >> dum; control >> extrapflag;
	control.close();

	vector<string> mfiles(nmesh), sfiles(nmesh), exfiles(nmesh);
	vector<double> h(nmesh,0), l2err(nmesh,0), dt(nmesh,0);
	string names[] = {"passive-scalar"};

	for(int i = 0; i < nmesh; i++) {
		mfiles[i] = meshprefix + to_string(i) + ".msh";
		sfiles[i] = outf + to_string(i) + "-h"+to_string(sdegree) + "-t"+to_string(tdegree)+".vtu";
		exfiles[i] = meshprefix + to_string(i) + "-exact.vtu";
	}

	for(int imesh = 0; imesh < nmesh; imesh++)
	{
		UMesh2dh m; m.readGmsh2(mfiles[imesh], NDIM); m.compute_topological(); m.compute_boundary_maps();
		
		// fixed time step is a constant times mesh size
		//double hh = sqrt( 1.0/m.gnelem() );
		double hh = m.meshSizeParameter();
		tstep = cfl*hh;
		printf("Mesh %d: h = %f, time step = %f\n", imesh, hh, tstep);

		Vector a(2); a[0] = a0; a[1] = a1;
		LinearAdvection sd(&m, sdegree, btype, a, 0.0, inoutflag, extrapflag);
		
		double (* inits[6])(double,double);
		inits[0] = &init; inits[1] = &initgradx; inits[2] = &initgrady; inits[3] = &initgradxx; inits[4] = initgradyy; inits[5] = initgradxy;
		if(btype == 't')
			sd.setInitialConditionModal(0, inits);
		else
			sd.setInitialConditionNodal(0, inits);

		TVDRKStepping td(&m, &sd, tdegree, ftime, cfl, 'c', tstep);

		double actual_ftime = td.integrate();
		sd.postprocess();
		l2err[imesh] = sd.computeL2Error(exactsol, actual_ftime);
		
		l2err[imesh] = log10(l2err[imesh]);
		h[imesh] = log10(hh);
		dt[imesh] = log10(tstep);
		printf("Mesh %d: Log mesh size = %f, log time step = %f, log L2 error = %f\n", imesh, h[imesh], dt[imesh], l2err[imesh]);

		const Array2d<a_real>& u = sd.getOutput();
		Array2d<a_real> vecs;
		writeScalarsVectorToVtu_PointData(sfiles[imesh], m, u, names, vecs, "none");

		/*Array2d<a_real> exactout(m.gnpoin(),1);
		for(a_int i = 0; i < m.gnpoin(); i++)
			exactout(i) = exactsol(m.gcoords(i,0),m.gcoords(i,1),actual_ftime);
		writeScalarsVectorToVtu_PointData(exfiles[imesh], m, exactout, names, vecs, "none");*/

		if(imesh > 0) {
			double l2slope = (l2err[imesh]-l2err[imesh-1])/(h[imesh]-h[imesh-1]);
			printf("L2 error slope (%d) = %f\n", imesh, l2slope);
		}
		cout << endl;
	}

	/*ofstream convf(outf);
	for(int i = 0; i < nmesh; i++)
		convf << h[i] << " " << dt[i] << " " << l2err[i] << "\n";
	convf.close();*/

	printf("---\n\n");
	return 0;
}
