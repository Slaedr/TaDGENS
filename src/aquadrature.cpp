/** @file aquadrature.cpp
 * @brief Initialization of quadrature rules
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#include "aquadrature.hpp"

namespace acfd {

void Quadrature1D::initialize(const int n_poly)
{
	using std::sqrt;

	nPoly = n_poly;
		gweights.resize(ngauss,1);
		gpoints.resize(ngauss,1);

	if(nPoly == 1) {
		ngauss = 1;
		gpoints(0) = 0.0;
		gweights(0) = 2.0;
	}
	else if(nPoly == 2){
		ngauss = 2;
		gpoints(0) = -1.0/sqrt(3); gpoints(1) = 1.0/sqrt(3);
		gweights(0) = 1.0; gweights(1) = 1.0;
	}
	else if(nPoly == 3) {
		ngauss = 3;
		gpoints(0) = -sqrt(3.0/5.0); gpoints(1) = 0.0; gpoints(2) = sqrt(3.0/5.0);
		gweights(0) = 5.0/9.0;  gweights(1) = 8.0/9.0, gweights(2) = 5.0/9.0;
	}
	else if(nPoly == 4) {
		ngauss = 4;
		gpoints(0) = -sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)); gpoints(1) = -sqrt(3.0/7 - 2.0/7*sqrt(6.0/5)); 
		gpoints(2) = sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)); gpoints(3) = sqrt(3.0/7 + 2.0/7*sqrt(6.0/5));
		gweights(0) = (18.0-sqrt(30))/36.0; gweights(1) = (18.0+sqrt(30))/36.0; 
		gweights(2) = (18.0+sqrt(30))/36.0; gweights(3) = (18.0-sqrt(30))/36.0;
	}
	else if(nPoly == 5) {
		ngauss = 5;
		gpoints(0) = -1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7)); gpoints(1) = -1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7));
		gpoints(2) = 0.0;
		gpoints(3) = 1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7)); gpoints(4) = 1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7));
		gweights(0) = (322.0-13*sqrt(70.0))/900; gweights(1) = (322.0+13*sqrt(70.0))/900;
		gweights(2) = 128.0/225;
		gweights(3) = (322.0+13*sqrt(70.0))/900; gweights(4) = (322.0-13*sqrt(70.0))/900;
	}
	else {
		ngauss = 5;
		gweights.resize(ngauss,1);
		gpoints.resize(ngauss,1);
		printf("! Quadrature1D: Quadrature with this number of Gauss points is not supported! Setting to 5.\n");
		gpoints(0) = -1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7)); gpoints(1) = -1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7));
		gpoints(2) = 0.0;
		gpoints(3) = 1.0/3*sqrt(5.0 - 2.0*sqrt(10.0/7)); gpoints(4) = 1.0/3*sqrt(5.0 + 2.0*sqrt(10.0/7));
		gweights(0) = (322.0-13*sqrt(70.0))/900; gweights(1) = (322.0+13*sqrt(70.0))/900;
		gweights(2) = 128.0/225;
		gweights(3) = (322.0+13*sqrt(70.0))/900; gweights(4) = (322.0-13*sqrt(70.0))/900;
	}
}

void Quadrature2DSquare::initialize(const int n_poly)
{
	nPoly = n_poly;
	if(nPoly == 1) {
		ngauss = 1;
		acfd_real gp[][2] = {{0.0, 0.0}};
		acfd_real gw[][1] = {{4.0}};
		/*acfd_real * gp[NDIM] = (acfd_real*[NDIM])malloc(ngauss*sizeof(acfd_real));
		acfd_real ** gw = (acfd_real**)malloc(ngauss*sizeof(acfd_real*));
		for(int i = 0; i < ngauss; i++) {
			gp[i] = (acfd_real*)malloc(NDIM*sizeof(acfd_real));
			gw[i] = (acfd_real*)malloc(sizeof(acfd_real));
		}
		gp = {{0.0, 0.0}};
		gw = {{4.0}};*/
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		std::printf("Quadrature2DSquare: Ngauss = 1.\n");
	}
	else {
		nPoly = 2;	// override argument
		ngauss = 4;
		acfd_real gp[][2] = {{-1.0/SQRT3, -1.0/SQRT3}, {-1.0/SQRT3, 1.0/SQRT3}, {1.0/SQRT3, -1.0/SQRT3}, {1.0/SQRT3, 1.0/SQRT3}};
		acfd_real gw[][1] = {{1.0}, {1.0}, {1.0}, {1.0}};
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		std::printf("Quadrature2DSquare: Ngauss = 4.\n");
	}
}

void Quadrature2DTriangle::initialize(const int n_poly)
{
	nPoly = n_poly;
	if(nPoly == 1) {
		ngauss = 1;
		acfd_real gp[][2] = {{1.0/3, 1.0/3}};
		acfd_real gw[][1] = {{0.5}};
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		printf("Quadrature2DTriangle: Ngauss = 1.\n");
	}
	else if(nPoly == 2) {
		ngauss = 3;
		acfd_real gp[][2] = {{0.6666666666667,0.1666666666667}, {0.1666666666667,0.6666666666667}, {0.1666666666667,0.1666666666667}};
		acfd_real gw[][1] = {{0.1666666666667}, {0.1666666666667}, {0.1666666666667}};
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		printf("Quadrature2DTriangle: Ngauss = 3.\n");
	}
	else if(nPoly == 3) {
		ngauss = 4;
		acfd_real gp[][2] = {{0.33333333333,0.33333333333}, {0.20000000000,0.20000000000}, {0.20000000000, 0.60000000000}, {0.60000000000, 0.20000000000}};
		acfd_real gw[][1] = {{-0.28125000000}, {0.26041666667}, {0.26041666667}, {0.26041666667}};
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		printf("Quadrature2DTriangle: Ngauss = 4.\n");
	}
	else if(nPoly == 4) {
		ngauss = 6;
		acfd_real gp[][2] = {{0.108103018168070,0.445948490915965},
							{0.445948490915965,0.108103018168070},
							{0.445948490915965,0.445948490915965},
							{0.816847572980459,0.091576213509771},
							{0.091576213509771,0.816847572980459},
							{0.091576213509771,0.091576213509771}};
		acfd_real gw[][1] = {{0.1116907948390055},
							{0.1116907948390055},
							{0.1116907948390055},
							{0.0549758718276610},
							{0.0549758718276610},
							{0.0549758718276610}};
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		printf("Quadrature2DTriangle: Ngauss = 6.\n");
	}
	else if(nPoly == 5) {
		ngauss = 7;
		acfd_real gp[][2] = {
							{0.33333333333333, 0.33333333333333},
							{0.47014206410511, 0.47014206410511},
							{0.47014206410511, 0.05971587178977},
							{0.05971587178977, 0.47014206410511},
							{0.10128650732346, 0.10128650732346},
							{0.10128650732346, 0.79742698535309},
							{0.79742698535309, 0.10128650732346}};
		acfd_real gw[][1] ={{0.22500000000000},
							{0.13239415278851},
							{0.13239415278851},
							{0.13239415278851},
							{0.12593918054483},
							{0.12593918054483},
							{0.12593918054483}};
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		printf("Quadrature2DTriangle: Ngauss = 7.\n");
	}
	else { /* npoly == 6 */
		ngauss = 12;
		acfd_real gp[][2] = {{0.873821971016996,0.063089014491502},
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
		acfd_real gw[][1] = {{0.0254224531851035},
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
		gpoints.initialize(ngauss, 2, (acfd_real*)gp);
		gweights.initialize(ngauss, 1, (acfd_real*)gw);
		printf("Quadrature2DTriangle: Ngauss = 12.\n");
	}
}

}
