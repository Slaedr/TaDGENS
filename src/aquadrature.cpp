/** @file aquadrature.cpp
 * @brief Initialization of quadrature rules
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#include "aquadrature.hpp"

namespace acfd {

/** Note that Gauss-Legendre quadrature (1D) with n quadrature points integrates
 * polynomials upto degree 2n-1 exactly.
 */
void Quadrature1D::initialize(const int n_poly)
{
	using std::sqrt;
	amat::Array2d<a_real> gptemp;

	nPoly = n_poly;
	shape = LINE;

	if(nPoly <= 1) {
		ngauss = 1;
		gweights.resize(ngauss,1);
		gptemp.resize(ngauss,1);
		ggpoints.resize(ngauss,1);
		gptemp(0) = 0.0;
		gweights(0) = 2.0;
	}
	else if(nPoly <= 3){
		ngauss = 2;
		gweights.resize(ngauss,1);
		gptemp.resize(ngauss,1);
		ggpoints.resize(ngauss,1);
		gptemp(0) = -1.0/sqrt(3); gptemp(1) = 1.0/sqrt(3);
		gweights(0) = 1.0; gweights(1) = 1.0;
	}
	else if(nPoly <= 5) {
		ngauss = 3;
		gweights.resize(ngauss,1);
		gptemp.resize(ngauss,1);
		ggpoints.resize(ngauss,1);
		gptemp(0) = -sqrt(3.0/5.0); gptemp(1) = 0.0; gptemp(2) = sqrt(3.0/5.0);
		gweights(0) = 5.0/9.0;  gweights(1) = 8.0/9.0, gweights(2) = 5.0/9.0;
	}
	else if(nPoly <= 7) {
		ngauss = 4;
		gweights.resize(ngauss,1);
		gptemp.resize(ngauss,1);
		ggpoints.resize(ngauss,1);
		gptemp(0) = -sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)); gptemp(1) = -sqrt(3.0/7 - 2.0/7*sqrt(6.0/5));
		gptemp(2) = sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)); gptemp(3) = sqrt(3.0/7 + 2.0/7*sqrt(6.0/5));
		gweights(0) = (18.0-sqrt(30))/36.0; gweights(1) = (18.0+sqrt(30))/36.0;
		gweights(2) = (18.0+sqrt(30))/36.0; gweights(3) = (18.0-sqrt(30))/36.0;
	}
	else if(nPoly <= 9) {
		ngauss = 5;
		gweights.resize(ngauss,1);
		gptemp.resize(ngauss,1);
		ggpoints.resize(ngauss,1);
		gptemp(0) = -1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7)); gptemp(1) = -1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7));
		gptemp(2) = 0.0;
		gptemp(3) = 1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7)); gptemp(4) = 1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7));
		gweights(0) = (322.0-13*sqrt(70.0))/900; gweights(1) = (322.0+13*sqrt(70.0))/900;
		gweights(2) = 128.0/225;
		gweights(3) = (322.0+13*sqrt(70.0))/900; gweights(4) = (322.0-13*sqrt(70.0))/900;
	}
	else {
		nPoly = 15;
		ngauss = 8;
		gweights.resize(ngauss,1);
		gptemp.resize(ngauss,1);
		ggpoints.resize(ngauss,1);
		printf("! Quadrature1D: Quadrature with this strength is not supported! Setting to 15.\n");
		
		a_real gp[][1] = {{-0.960289856497536231684},
							{-0.796666477413626739592},
							{-0.525532409916328985818},
							{-0.183434642495649804939},
							{0.183434642495649804939},
							{0.525532409916328985818},
							{0.796666477413626739592},
							{0.960289856497536231684}};

		a_real gw[][1] = {{0.10122853629037625915},
		                  {0.22238103445337447054},
		                  {0.31370664587788728733},
		                  {0.36268378337836198296},
		                  {0.36268378337836198296},
		                  {0.31370664587788728733},
		                  {0.22238103445337447054},
		                  {0.10122853629037625915}};
		gptemp.initialize(ngauss, 1, (a_real*)gp);
		gweights.initialize(ngauss, 1, (a_real*)gw);
	}

	for(int i = 0; i < ggpoints.rows(); i++)
		ggpoints(i,0) = gptemp(i,0);
}

void Quadrature2DSquare::initialize(const int n_poly)
{
	using std::sqrt;
	amat::Array2d<a_real> gptemp;
	amat::Array2d<a_real> gwtemp;

	nPoly = n_poly;
	shape = QUADRANGLE;
	int ngaussdim;

	if(nPoly <= 1) {
		ngaussdim = 1;
		ngauss = 1;
		a_real gp[] = {0.0};
		a_real gw[] = {2.0};
		gptemp.initialize(ngaussdim, 1, gp);
		gwtemp.initialize(ngaussdim, 1, gw);
	}
	else if(nPoly <= 3){
		ngaussdim = 2;
		ngauss = 4;
		a_real gp[] = {-1.0/sqrt(3), 1.0/sqrt(3) };
		a_real gw[] = { 1.0, 1.0 };
		gptemp.initialize(ngaussdim, 1, gp);
		gwtemp.initialize(ngaussdim, 1, gw);
	}
	else if(nPoly <= 5) {
		ngaussdim = 3;
		ngauss = 9;
		a_real gp[] = { -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0) };
		a_real gw[] = { 5.0/9.0, 8.0/9.0, 5.0/9.0 };
		gptemp.initialize(ngaussdim, 1, gp);
		gwtemp.initialize(ngaussdim, 1, gw);
	}
	else if(nPoly <= 7) {
		ngaussdim = 4;
		ngauss = 16;
		a_real gp[] = {-sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)), -sqrt(3.0/7 - 2.0/7*sqrt(6.0/5)),
		               sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)), sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)) };
		a_real gw[] = { (18.0-sqrt(30))/36.0, (18.0+sqrt(30))/36.0,
		                (18.0+sqrt(30))/36.0, (18.0-sqrt(30))/36.0 };
		gptemp.initialize(ngaussdim, 1, gp);
		gwtemp.initialize(ngaussdim, 1, gw);
	}
	else if(nPoly <= 9) {
		ngaussdim = 5;
		ngauss = 25;
		a_real gp[] = { -1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7)), -1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7)), 0.0,
		                1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7)), 1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7)) };
		a_real gw[] = {(322.0-13*sqrt(70.0))/900, (322.0+13*sqrt(70.0))/900, 128.0/225,
		               (322.0+13*sqrt(70.0))/900, (322.0-13*sqrt(70.0))/900 };
		gptemp.initialize(ngaussdim, 1, gp);
		gwtemp.initialize(ngaussdim, 1, gw);
	}
	else {
		nPoly = 15;
		ngaussdim = 8;
		ngauss = 64;
		printf("! Quadrature1D: Quadrature with this strength is not supported! Setting to 15th degree polynomial.\n");
		
		a_real gp[] =      {-0.960289856497536231684,
							-0.796666477413626739592,
							-0.525532409916328985818,
							-0.183434642495649804939,
							0.183434642495649804939,
							0.525532409916328985818,
							0.796666477413626739592,
							0.960289856497536231684};

		a_real gw[] =    {0.10122853629037625915,
		                  0.22238103445337447054,
		                  0.31370664587788728733,
		                  0.36268378337836198296,
		                  0.36268378337836198296,
		                  0.31370664587788728733,
		                  0.22238103445337447054,
		                  0.10122853629037625915};
		gptemp.initialize(ngaussdim, 1, gp);
		gwtemp.initialize(ngaussdim, 1, gw);
	}

	gweights.resize(ngauss,1);
	ggpoints.resize(ngauss,2);
	for(int i = 0; i < ngaussdim; i++)
	{
		for(int j = 0; j < ngaussdim; j++){
			ggpoints(i*ngaussdim+j,0) = gptemp(i);
			ggpoints(i*ngaussdim+j,1) = gptemp(j);
			gweights(i*ngaussdim+j) = gwtemp(i)*gwtemp(j);
		}
	}
}

void Quadrature2DTriangle::initialize(const int n_poly)
{
	amat::Array2d<a_real> gptemp;
	nPoly = n_poly;
	shape = TRIANGLE;

	if(nPoly == 1) {
		ngauss = 1;
		a_real gp[][2] = {{1.0/3, 1.0/3}};
		a_real gw[][1] = {{0.5}};
		gptemp.initialize(ngauss, 2, (a_real*)gp);
		gweights.initialize(ngauss, 1, (a_real*)gw);
		printf("  Quadrature2DTriangle: Ngauss = 1.\n");
		ggpoints.resize(ngauss,2);
	}
	else if(nPoly == 2) {
		ngauss = 3;
		a_real gp[][2] = {{0.6666666666667,0.1666666666667}, {0.1666666666667,0.6666666666667}, {0.1666666666667,0.1666666666667}};
		a_real gw[][1] = {{0.1666666666667}, {0.1666666666667}, {0.1666666666667}};
		gptemp.initialize(ngauss, 2, (a_real*)gp);
		gweights.initialize(ngauss, 1, (a_real*)gw);
		printf("  Quadrature2DTriangle: Ngauss = 3.\n");
		ggpoints.resize(ngauss,2);
	}
	else if(nPoly == 3) {
		ngauss = 4;
		a_real gp[][2] = {{0.33333333333,0.33333333333}, {0.20000000000,0.20000000000}, {0.20000000000, 0.60000000000}, {0.60000000000, 0.20000000000}};
		a_real gw[][1] = {{-0.28125000000}, {0.26041666667}, {0.26041666667}, {0.26041666667}};
		gptemp.initialize(ngauss, 2, (a_real*)gp);
		gweights.initialize(ngauss, 1, (a_real*)gw);
		printf("  Quadrature2DTriangle: Ngauss = 4.\n");
		ggpoints.resize(ngauss,2);
	}
	else if(nPoly == 4) {
		ngauss = 6;
		a_real gp[][2] = {{0.108103018168070,0.445948490915965},
							{0.445948490915965,0.108103018168070},
							{0.445948490915965,0.445948490915965},
							{0.816847572980459,0.091576213509771},
							{0.091576213509771,0.816847572980459},
							{0.091576213509771,0.091576213509771}};
		a_real gw[][1] = {{0.1116907948390055},
							{0.1116907948390055},
							{0.1116907948390055},
							{0.0549758718276610},
							{0.0549758718276610},
							{0.0549758718276610}};
		gptemp.initialize(ngauss, 2, (a_real*)gp);
		gweights.initialize(ngauss, 1, (a_real*)gw);
		printf("  Quadrature2DTriangle: Ngauss = 6.\n");
		ggpoints.resize(ngauss,2);
	}
	else if(nPoly == 5) {
		ngauss = 7;
		a_real gp[][2] = {
							{0.33333333333333, 0.33333333333333},
							{0.47014206410511, 0.47014206410511},
							{0.47014206410511, 0.05971587178977},
							{0.05971587178977, 0.47014206410511},
							{0.10128650732346, 0.10128650732346},
							{0.10128650732346, 0.79742698535309},
							{0.79742698535309, 0.10128650732346}};
		a_real gw[][1] ={{0.22500000000000},
							{0.13239415278851},
							{0.13239415278851},
							{0.13239415278851},
							{0.12593918054483},
							{0.12593918054483},
							{0.12593918054483}};
		gptemp.initialize(ngauss, 2, (a_real*)gp);
		gweights.initialize(ngauss, 1, (a_real*)gw);
		printf("  Quadrature2DTriangle: Ngauss = 7.\n");
		ggpoints.resize(ngauss,2);
	}
	else { /* npoly == 6 */
		nPoly = 6;
		ngauss = 12;
		a_real gp[][2] = {{0.873821971016996,0.063089014491502},
							{0.063089014491502,0.873821971016996},
							{0.063089014491502,0.063089014491502},
							{0.501426509658179,0.249286745170910},
							{0.249286745170910,0.501426509658179},
							{0.249286745170910,0.249286745170910},
							{0.636502499121399,0.310352451033785},
							{0.636502499121399,0.053145049844816},
							{0.310352451033785,0.636502499121399},
							{0.310352451033785,0.053145049844816},
							{0.053145049844816,0.310352451033785},
							{0.053145049844816,0.636502499121399}};
		a_real gw[][1] = {{0.0254224531851035},
							{0.0254224531851035},
							{0.0254224531851035},
							{0.0583931378631895},
							{0.0583931378631895},
							{0.0583931378631895},
							{0.0414255378091870},
							{0.0414255378091870},
							{0.0414255378091870},
							{0.0414255378091870},
							{0.0414255378091870},
							{0.0414255378091870}};
		gptemp.initialize(ngauss, 2, (a_real*)gp);
		gweights.initialize(ngauss, 1, (a_real*)gw);
		printf("  Quadrature2DTriangle: Ngauss = 12.\n");
		ggpoints.resize(ngauss,2);
	}
	for(int i = 0; i < ggpoints.rows(); i++)
		for(int j = 0; j < 2; j++)
			ggpoints(i,j) = gptemp(i,j);
}

}
