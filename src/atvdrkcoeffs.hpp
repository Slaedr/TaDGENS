/** @file aode.hpp
 * @brief Provides coefficients in the Butcher tableau of ODE integrators
 * @author Aditya Kashi
 * @date 2017-03-04
 */

#ifndef __AODE_H
#define __AODE_H

#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

namespace acfd {

amat::Array2d<double> tvdrk1(1,3);
tvdrk1[0][0] = 1.0;	tvdrk1[0][1] = 0.0; tvdrk1[0][2] = 1.0;

amat::Array2d<double> tvdrk2(2,3);
tvdrk2[0][0] = 1.0;	tvdrk2[0][1] = 0.0; tvdrk2[0][2] = 1.0;
tvdrk2[1][0] = 0.5;	tvdrk2[1][1] = 0.5;	tvdrk2[1][2] = 0.5;

amat::Array2d<double> tvdrk3(3,3);
tvdrk3[0][0] = 1.0;		tvdrk3[0][1] = 0.0; 	tvdrk3[0][2] = 1.0;
tvdrk3[1][0] = 0.75;	tvdrk3[1][1] = 0.25;	tvdrk3[1][2] = 0.25;
tvdrk3[2][0] = 0.3333333333333333; tvdrk3[2][1] = 0.6666666666666667; tvdrk3[2][2] = 0.6666666666666667;

}
#endif
