/** @file elements.cpp
 * @brief Implementation of finite elements.
 * @author Aditya Kashi
 * @date 2017 March 1
 */

#include "elements.hpp"

namespace acfd {

void LagrangeMapping1D::computeSpeedsAndNormals(const amat::Array2d<acfd_real>& points)
{
	int npoin = points.rows();
	speeds.resize(npoin);
	normals.resize(npoin,NDIM);

	acfd_real vel[NDIM];

	if(order == 1) {
		int ndof = 2;
		for(int i = 0; i < npoin; i++)
		{
			// get sum Lagrange derivatives multiplied by coeffs
			for(int idim = 0; idim < NDIM; idim++)
				vel[idim] = (phyNodes[1][idim] - phyNodes[0][idim])/2.0;
			speeds[i] = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
			// normalize velocity to get unit tangent
			for(int idim = 0; idim < NDIM; idim++)
				vel[idim] /= speeds[i];
			normals[i][0] = vel[1]; normals[i][1] = -vel[0];
		}
	}
	else if(order == 2) {
		int ndof = 3;
		// TODO: Implement P2 tangent and normal
	}
	else
		std::printf("! LagrangeMapping1D: Chosen geometric order not available!\n");
}

}
