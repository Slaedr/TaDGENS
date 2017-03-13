/** @file elements.cpp
 * @brief Implementation of finite elements.
 * @author Aditya Kashi
 * @date 2017 March 1
 */

#include "elements.hpp"

namespace acfd {

using namespace amat;

/** Currently, Lagrange mappings upto polynomial degree 2 are implemented.
 */
void LagrangeMapping1D::computeAll()
{
	const Array2d<acfd_real>& points = quadrature->points();
	int npoin = points.rows();
	//speeds.resize(npoin);
	normals.resize(npoin,NDIM);
	mapping.resize(NDIM);

	acfd_real vel[NDIM];

	if(degree == 1) {
		for(int i = 0; i < npoin; i++)
		{
			for(int idim = 0; idim < NDIM; idim++) {
				// mapping - num of Lagrange shape functions mulitplied by resp coeffs
				mapping[i](idim) = phyNodes(0,idim)*(1.0-points(i))*0.5 + phyNodes(1,idim)*(1.0+points(i))*0.5;
				// get sum Lagrange derivatives multiplied by coeffs
				vel[idim] = (phyNodes[1][idim] - phyNodes[0][idim])/2.0;
			}
			
			//speeds[i] = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
			// (no need to) normalize velocity to get unit tangent
			/*for(int idim = 0; idim < NDIM; idim++)
				vel[idim] /= speeds[i];*/

			normals[i][0] = vel[1]; normals[i][1] = -vel[0];
		}
	}
	else if(degree == 2) {
		for(int i = 0; i < npoin; i++)
		{
			for(int idim = 0; idim < NDIM; idim++) {
				mapping[ip](idim) = phyNodes(0,idim)*points(i)*(points(i)-1.0)/2 + phyNodes(1,idim)*points(i)*(points(i)+1.0)/2 + phyNodes(2,idim)*(1.0-points(i)*points(i));
				vel[idim] = phyNodes[0][idim]*(points[i]-0.5) + phyNodes[1][idim]*(points[i]+0.5) + phyNodes[2][idim]*(-2.0*points[i]);
			}

			normals[i][0] = vel[1]; normals[i][1] = -vel[0];
		}
	}
	/* Add another else-if block above this line to add a higher-order map.
	 */
	else
		std::cout << "! LagrangeMapping1D: Chosen geometric order not available!\n";
}

/** Mappings upto P2 are implemented.
 */
void LagrangeMapping2DTriangle::computeAll()
{
	const Array2d<acfd_real>& points = quadrature->points();
	shape = TRIANGLE;
	int npoin = points.rows();
	jaco.resize(npoin);
	jacoinv.resize(npoin);
	jacodet.resize(npoin);
	mapping.resize(npoin,NDIM);

	if(degree == 1) {
		for(int ip = 0; ip < npoin; ip++)
		{
			jaco[ip].resize(NDIM,NDIM);
			jacoinv[ip].resize(NDIM,NDIM);
			for(int idim = 0; idim < NDIM; idim++) 
			{
				mapping(ip,idim) = phyNodes(0,idim)*(1.0-points(ip,0)-points(ip,1)) + phyNodes(1,idim)*points(ip,0) + phyNodes(2,idim)*points(ip,1);
				jaco[ip](idim,0) = phyNodes(1,idim)-phyNodes(0,idim);
				jaco[ip](idim,1) = phyNodes(2,idim)-phyNodes(0,idim);
			}
			jacodet[ip] = jaco[ip](0,0)*jaco[ip](1,1) - jaco[ip](0,1)*jaco[ip](1,0)
			jacoinv[ip](0,0) = jaco[ip](1,1)/jacodet[ip]; jacoinv[ip](0,1) = -jaco[ip](0,1)/jacodet[ip]
			jacoinv[ip](1,0) = -jaco[ip](1,0)/jacodet[ip]; jacoinv[ip](1,1) = jaco[ip](0,0)/jacodet[ip]
		}
	}
	else if(degree == 2) {
		for(int ip = 0; ip < npoin; ip++)
		{
			jaco[ip].resize(NDIM,NDIM);
			jacoinv[ip].resize(NDIM,NDIM);
			for(int idim = 0; idim < NDIM; idim++) 
			{
				mapping(ip,idim) = phynodes(0,idim) * (1.0-3*points(ip,0)-3*points(ip,1)+2*points(ip,0)*points(ip,0)+2*points(ip,1)*points(ip,1)+4*points(ip,0)*points(ip,1)) 
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

void LagrangeMapping2DTriangle::computeMappingAndJacobianDet()
{
	const Array2d<acfd_real>& points = quadrature->points();
	shape = TRIANGLE;
	int npoin = points.rows();
	jacodet.resize(npoin);
	mapping.resize(npoin,NDIM);

	Matrix jacol.resize(NDIM,NDIM);

	if(degree == 1) {
		for(int ip = 0; ip < npoin; ip++)
		{
			for(int idim = 0; idim < NDIM; idim++) 
			{
				mapping(ip,idim) = phyNodes(0,idim)*(1.0-points(ip,0)-points(ip,1)) + phyNodes(1,idim)*points(ip,0) + phyNodes(2,idim)*points(ip,1);
				jacol(idim,0) = phyNodes(1,idim)-phyNodes(0,idim);
				jacol(idim,1) = phyNodes(2,idim)-phyNodes(0,idim);
			}
			jacodet[ip] = jacol(0,0)*jacol(1,1) - jacol(0,1)*jacol(1,0)
		}
	}
	else if(degree == 2) {
		for(int ip = 0; ip < npoin; ip++)
		{
			for(int idim = 0; idim < NDIM; idim++) 
			{
				mapping(ip,idim) = phynodes(0,idim) * (1.0-3*points(ip,0)-3*points(ip,1)+2*points(ip,0)*points(ip,0)+2*points(ip,1)*points(ip,1)+4*points(ip,0)*points(ip,1)) 
					+ phynodes(1,idim)*(2.0*points(ip,0)*points(ip,0)-points(ip,0)) 
					+ phynodes(2,idim)*(2.0*points(ip,1)*points(ip,1)-points(ip,1)) 
					+ phynodes(3,idim)*4.0*(points(ip,0)-points(ip,0)*points(ip,0)-points(ip,0)*points(ip,1)) 
					+ phynodes(4,idim)*4.0*points(ip,0)*points(ip,1) 
					+ phynodes(5,idim)*4.0*(points(ip,1)-points(ip,1)*points(ip,1)-points(ip,0)*points(ip,1));
				jacol(idim,0) = phyNodes(0,idim)*(-3+4*points(ip,0)+4*points(ip,1)) +phyNodes(1,idim)*(4*points(ip,0)-1) +phyNodes(3,idim)*4*(1-2*points(ip,0)-points(ip,1)) 
					+phyNodes(4,idim)*4*points(ip,1) -phyNodes(5,idim)*4.0*points(ip,1);
				jacol(idim,1) = phyNodes(0,idim)*(-3+4*points(ip,1)+4*points(ip,0)) +phyNodes(2,idim)*(4*points(ip,1)-1) -phyNodes(3,idim)*4*points(ip,0) 
					+ phyNodes(4,idim)*4*points(ip,0) +phyNodes(5,idim)*4*(1-2*points(ip,1)-points(ip,0));
			}
			jacodet[ip] = jacol(0,0)*jacol(1,1) - jacol(0,1)*jacol(1,0)
		}
	}
}

void LagrangeMapping2DQuadrangle::computeAll()
{
	const Array2d<acfd_real>& qp = quadrature->points();
	// TODO: Add geometric computations
}

void LagrangeMapping2DQuadrangle::computeMappingAndJacobianDet()
{
	const Array2d<acfd_real>& qp = quadrature->points();
	// TODO: Add mapping and jaco det
}

/** We currently have upto P2 elements.
 * The number of DOFs is computed as \f$ \sum_{i=1}^{p+1} i \f$ for p = 0,1,2...
 * For computing the element centers and area, we use enough quadrature points to exactly get the element centers.
 * Note that for a quad element of degree bi p (p=1 is bi linear etc), the jacodet is of degree bi 2p-1.
 * For a tri element of degree p, the jacodet is of degree 2p-2.
 */
void TaylorElement::initialize(int degr, const GeomMapping2D* geommap)
{
	degree = degr; gmap = geommap;
	ndof = 0;
	for(int i = 1; i <= degree+1; i++)
		ndof += i;

	int ngauss = gmap->getQuadrature()->numGauss();
	basis.resize(ngauss);
	basisGrad.resize(ngauss);
	for(int i = 0; i < ngauss; i++) {
		basis[i].setup(ndof,1);
		basisGrad[i].resize(ndof,NDIM);
	}
	
	int geomdeg = gmap->getDegree();
	area = 0;

	if(degree > 1)
		basisOffset.resize(degree-1);
	for(int i = 0; i < degree-1; i++) {
		basisOffset[i].resize(i+3);
		for(int j = 0; j < i+3; j++)
			basisOffset[i][j] = 0;
	}

	for(int i = 0; i<NDIM; i++) {
		center[i] = 0;
		delta[i] = 0;
	}

	// Compute max extents of the physical element
	
	const Array2d<acfd_real>& phn = gmap->getPhyNodes();
	for(int i = 0; i < phn.rows(); i++) {
		for(int j = i+1; j < phn.rows(); j++) 
		{
			acfd_real dist[NDIM];
			for(int idim = 0; idim < NDIM; idim++) 
			{
				dist[idim] = phn(i,idim)-phn(j,idim);
				if(dist[idim] > delta[idim])
					delta[idim] = dist[idim];
			}
		}
	}
	
	// Compute element centers and basis offsets

	const Array2d<acfd_real>& gp = gmap->getQuadrature()->points();
	const Array2d<acfd_real>& gw = gmap->getQuadrature()->weights();
	int ng = gp.rows();
	for(int ig = 0; ig < ng; ig++) 
	{
		area += gmap->jacDet(ig) * gw(ig);
		for(int idim = 0; idim < NDIM; idim++)
			center[idim] += gmap->map(ig)(idim) * gmap->jacDet(ig) * gw(ig);
	}
	for(int idim = 0; idim < NDIM; idim++)
		center[idim] /= area;

	if(degree >= 2) {
		for(int ig = 0; ig < ng; ig++)
		{
			basisOffset[0][0] += (gmap->map(ig)(0)-center[0])*(gmap->map(ig)(0)-center[0]) * gmap->jacDet(ig) * gw(ig);
			basisOffset[0][1] += (gmap->map(ig)(1)-center[1])*(gmap->map(ig)(1)-center[1]) * gmap->jacDet(ig) * gw(ig);
			basisOffset[0][2] += (gmap->map(ig)(0)-center[0])*(gmap->map(ig)(1)-center[1]) * gmap->jacDet(ig) * gw(ig);
		}
		basisOffset[0][0] *= 1.0/(area*2*delta[0]*delta[0]);
		basisOffset[0][1] *= 1.0/(area*2*delta[1]*delta[1]);
		basisOffset[0][2] *= 1.0/(area*delta[0]*delta[1]);
	}
	
	// Compute basis functions and gradients
	
	for(int ip = 0; ip < ngauss; ip++)
	{
		// get physical coords of quadrature point
		const Array2d<acfd_real>& gp = gmap->map();
		
		basis[ip](0) = 1.0;
		basisGrad[ip](0,0) = basisGrad[ip](0,1) = 0.0;

		if(degree >= 1) {
			basis[ip](1) = (gp(ip,0)-center[0])/delta[0];
			basis[ip](2) = (gp(ip,1)-center[1])/delta[1];

			basisGrad[ip](1,0) = 1.0/delta[0]; basisGrad[ip](1,1) = 0.0;
			basisGrad[1p](2,0) = 0.0; basisGrad[ip](2,1) = 1.0/delta[1];
		}

		if(degree >= 2) {
			basis[ip](3) = (gp(ip,0)-center[0])*(gp(ip,0)-center[0])/(2.0*delta[0]*delta[0]) - basisOffset[0][0];
			basis[ip](4) = (gp(ip,1)-center[1])*(gp(ip,1)-center[1])/(2.0*delta[1]*delta[1]) - basisOffset[0][1];
			basis[ip](5) = (gp(ip,0)-center[0])*(gp(ip,1)-center[1])/(delta[0]*delta[1]) - basisOffset[0][2];

			basisGrad[ip](3,0) = (gp(ip,0)-center[0])/(delta[0]*delta[0]); basisGrad[ip](3,1) = 0.0;
			basisGrad[ip](4,0) = 0.0; basisGrad[ip](4,1) = (gp(ip,1)-center[1])/(delta[1]*delta[1]);
			basisGrad[ip](5,0) = (gp(ip,1)-center[1])/(delta[0]*delta[1]); basisGrad[ip](5,1) = (gp(ip,0)-center[0])/(delta[0]*delta[1]);
		}
	}
}

void TaylorElement::computeBasis(const acfd_real* gp, acfd_real* basis)
{
	basis[0] = 1.0;

	if(degree >= 1) {
		basis[1] = (gp[0]-center[0])/delta[0];
		basis[2] = (gp[1]-center[1])/delta[1];
	}

	if(degree >= 2) {
		basis[3] = (gp[0]-center[0])*(gp[0]-center[0])/(2.0*delta[0]*delta[0]) - basisOffset[0][0];
		basis[4] = (gp[1]-center[1])*(gp[1]-center[1])/(2.0*delta[1]*delta[1]) - basisOffset[0][1];
		basis[5] = (gp[0]-center[0])*(gp[1]-center[1])/(delta[0]*delta[1]) - basisOffset[0][2];
	}
}

void FaceElement_PhysicalSpace::initialize(int degr, const Element_PhysicalSpace* lelem, const Element_PhysicalSpace* relem, const GeomMapping1D* gmapping)
{
	gmap = gmapping; leftel = lelem; rightel = relem;
	const Array2d<acfd_real>& gp = gmap->getQuadrature()->points();
	int ng = gmap->getQuadrature()->numGauss();
	int ldeg = leftel->getDegree(); int rdeg = rightel->getDegree();

	leftbasis.resize(ng,lelem->getNumDOFs()); rightbasis.resize(ng,relem->getNumDOFs());

	for(int ig = 0; ig < ng; ig++)
	{
		const Array2d<acfd_real>& points = gmap->map();
		lelem->computeBasis(points[ig], leftbasis[ig]);
		relem->computeBasis(points[ig], rightbasis[ig]);
	}
}

}
