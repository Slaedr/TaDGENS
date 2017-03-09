/** @file elements.cpp
 * @brief Implementation of finite elements.
 * @author Aditya Kashi
 * @date 2017 March 1
 */

#include "elements.hpp"

namespace acfd {

/** Currently, Lagrange mappings upto polynomial degree 2 are implemented.
 */
void LagrangeMapping1D::computeAll(const amat::Array2d<acfd_real>& points)
{
	int npoin = points.rows();
	speeds.resize(npoin);
	normals.resize(npoin,NDIM);
	mapping.resize(NDIM);

	acfd_real vel[NDIM];

	if(order == 1) {
		for(int i = 0; i < npoin; i++)
		{
			for(int idim = 0; idim < NDIM; idim++) {
				// mapping
				mapping[i](idim) = phyNodes(0,idim)*(1.0-points(i))*0.5 + phyNodes(1,idim)*(1.0+points(i))*0.5;
				// get sum Lagrange derivatives multiplied by coeffs
				vel[idim] = (phyNodes[1][idim] - phyNodes[0][idim])/2.0;
			}
			
			speeds[i] = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
			
			// normalize velocity to get unit tangent
			for(int idim = 0; idim < NDIM; idim++)
				vel[idim] /= speeds[i];
			normals[i][0] = vel[1]; normals[i][1] = -vel[0];
		}
	}
	else if(order == 2) {
		for(int i = 0; i < npoin; i++)
		{
			for(int idim = 0; idim < NDIM; idim++) {
				mapping[ip](idim) = phyNodes(0,idim)*points(i)*(points(i)-1.0)/2 + phyNodes(1,idim)*points(i)*(points(i)+1.0)/2 + phyNodes(2,idim)*(1.0-points(i)*points(i));
				vel[idim] = phyNodes[0][idim]*(points[i]-0.5) + phyNodes[1][idim]*(points[i]+0.5) + phyNodes[2][idim]*(-2.0*points[i]);
			}
			
			speeds[i] = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
			
			// normalize velocity to get unit tangent
			for(int idim = 0; idim < NDIM; idim++)
				vel[idim] /= speeds[i];
			normals[i][0] = vel[1]; normals[i][1] = -vel[0];
		}
	}
	/* Add another else-if block above this line to add a higher-order map.
	 */
	else
		std::cout << "! LagrangeMapping1D: Chosen geometric order not available!\n";
}

LagrangeMapping2DTriangle::LagrangeMapping2DTriangle() {
}

/** Mappings upto P2 are implemented.
 */
void LagrangeMapping2DTriangle::computeAll(const amat::Array2d<acfd_real>& points)
{
	int npoin = points.rows();
	jaco.resize(npoin);
	jacoinv.resize(npoin);
	jacodet.resize(npoin);
	mapping.resize(npoin);

	if(order == 1) {
		for(int ip = 0; ip < npoin; ip++)
		{
			mapping[ip].resize(NDIM);
			jaco[ip].resize(NDIM,NDIM);
			jacoinv[ip].resize(NDIM,NDIM);
			for(int idim = 0; idim < NDIM; idim++) 
			{
				mapping[ip](idim) = phyNodes(0,idim)*(1.0-points(ip,0)-points(ip,1)) + phyNodes(1,idim)*points(ip,0) + phyNodes(2,idim)*points(ip,1);
				jaco[ip](idim,0) = phyNodes(1,idim)-phyNodes(0,idim);
				jaco[ip](idim,1) = phyNodes(2,idim)-phyNodes(0,idim);
			}
			jacodet[ip] = jaco[ip](0,0)*jaco[ip](1,1) - jaco[ip](0,1)*jaco[ip](1,0)
			jacoinv[ip](0,0) = jaco[ip](1,1)/jacodet[ip]; jacoinv[ip](0,1) = -jaco[ip](0,1)/jacodet[ip]
			jacoinv[ip](1,0) = -jaco[ip](1,0)/jacodet[ip]; jacoinv[ip](1,1) = jaco[ip](0,0)/jacodet[ip]
		}
	}
	else if(order == 2) {
		for(int ip = 0; ip < npoin; ip++)
		{
			mapping[ip].resize(NDIM);
			jaco[ip].resize(NDIM,NDIM);
			jacoinv[ip].resize(NDIM,NDIM);
			for(int idim = 0; idim < NDIM; idim++) 
			{
				mapping[ip](idim) = phynodes(0,idim) * (1.0-3*points(ip,0)-3*points(ip,1)+2*points(ip,0)*points(ip,0)+2*points(ip,1)*points(ip,1)+4*points(ip,0)*points(ip,1)) 
					+ phynodes(1,idim)*(2.0*points(ip,0)*points(ip,0)-points(ip,0)) 
					+ phynodes(2,idim)*(2.0*points(ip,1)*points(ip,1)-points(ip,1)) 
					+ phynodes(3,idim)*4.0*(points(ip,0)-points(ip,0)*points(ip,0)-points(ip,0)*points(ip,1)) 
					+ phynodes(4,idim)*4.0*points(ip,0)*points(ip,1) 
					+ phynodes(5,idim)*4.0*(points(ip,1)-points(ip,1)*points(ip,1)-points(ip,0)*points(ip,1));
				jaco[ip](idim,0) = phyNodes(0,idim)*(-3+4*points(ip,0)+4*points(ip,1)) +phyNodes(1,idim)*(4*points(ip,0)-1) +phyNodes(3,idim)*4*(1-2*points(ip,0)-points(ip,1)) 
					+phyNodes(4,idim)*4*points(ip,1) -phyNodes(5,idim)*4.0*points(ip,1);
				jaco[ip](idim,1) = phyNodes(0,idim)*(-3+4*points(ip,1)+4*points(ip,0)) +phyNodes(2,idim)*(4*points(ip,1)-1) -phyNodes(3,idim)*4*points(ip,0) 
					+ phyNodes(4,idim)*4*points(ip,0) +phyNodes(5,idim)*4*(1-2*points(ip,1)-points(ip,0));
			}
			jacodet[ip] = jaco[ip](0,0)*jaco[ip](1,1) - jaco[ip](0,1)*jaco[ip](1,0)
			jacoinv[ip](0,0) = jaco[ip](1,1)/jacodet[ip]; jacoinv[ip](0,1) = -jaco[ip](0,1)/jacodet[ip]
			jacoinv[ip](1,0) = -jaco[ip](1,0)/jacodet[ip]; jacoinv[ip](1,1) = jaco[ip](0,0)/jacodet[ip]
		}
	}
}

LagrangeMapping2DQuadrangle::LagrangeMapping2DQuadrangle() {
}

}
