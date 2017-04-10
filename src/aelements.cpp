/** @file elements.cpp
 * @brief Implementation of finite elements.
 * @author Aditya Kashi
 * @date 2017 March 1
 */

#include "aelements.hpp"

namespace acfd {

using namespace amat;

/** Currently, Lagrange mappings upto polynomial degree 2 are implemented.
 */
void LagrangeMapping1D::computeAll()
{
	const Array2d<a_real>& points = quadrature->points();
	int npoin = points.rows();
	speeds.resize(npoin);
	normals.resize(npoin);
	for(int i = 0; i < npoin; i++)
		normals[i].resize(NDIM);
	mapping.resize(npoin,NDIM);

	a_real vel[NDIM];

	if(degree == 1) {
		for(int i = 0; i < npoin; i++)
		{
			for(int idim = 0; idim < NDIM; idim++) {
				// mapping - num of Lagrange shape functions mulitplied by resp coeffs
				mapping[i][idim] = phyNodes(0,idim)*(1.0-points(i))*0.5 + phyNodes(1,idim)*(1.0+points(i))*0.5;
				// get sum Lagrange derivatives multiplied by coeffs
				vel[idim] = (phyNodes[1][idim] - phyNodes[0][idim])/2.0;
			}
			
			speeds[i] = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
			
			// normalize tangent to get unit tangent and thus unit normal
			for(int idim = 0; idim < NDIM; idim++)
				vel[idim] /= speeds[i];
			normals[i][0] = vel[1]; normals[i][1] = -vel[0];
		}
	}
	else if(degree == 2) {
		for(int i = 0; i < npoin; i++)
		{
			for(int idim = 0; idim < NDIM; idim++) {
				mapping[i][idim] = phyNodes(0,idim)*points(i)*(points(i)-1.0)/2 + phyNodes(1,idim)*points(i)*(points(i)+1.0)/2 + phyNodes(2,idim)*(1.0-points(i)*points(i));
				vel[idim] = phyNodes[0][idim]*(points(i)-0.5) + phyNodes[1][idim]*(points(i)+0.5) + phyNodes[2][idim]*(-2.0*points(i));
			}
			
			speeds[i] = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);

			// normalize tangent to get unit tangent and thus unit normal
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

void LagrangeMapping2D::computeJacobians(const amat::Array2d<a_real>& __restrict__ po, std::vector<Matrix>& __restrict__ jac, 
		std::vector<Matrix>& __restrict__ jacoi, std::vector<a_real>& __restrict__ jacod)
{
	if (shape == TRIANGLE)
	{
		if(degree == 1) {
			for(int ip = 0; ip < po.rows(); ip++) {
				for(int idim = 0; idim < NDIM; idim++) 
				{
					jac[ip](idim,0) = phyNodes(1,idim)-phyNodes(0,idim);
					jac[ip](idim,1) = phyNodes(2,idim)-phyNodes(0,idim);
				}
			}
		}
		else if(degree == 2) {
			for(int ip = 0; ip < po.rows(); ip++) {
				for(int idim = 0; idim < NDIM; idim++) 
				{
					jac[ip](idim,0) = phyNodes(0,idim)*(-3+4*po(ip,0)+4*po(ip,1)) +phyNodes(1,idim)*(4*po(ip,0)-1) +phyNodes(3,idim)*4*(1-2*po(ip,0)-po(ip,1)) 
						+phyNodes(4,idim)*4*po(ip,1) -phyNodes(5,idim)*4.0*po(ip,1);
					jac[ip](idim,1) = phyNodes(0,idim)*(-3+4*po(ip,1)+4*po(ip,0)) +phyNodes(2,idim)*(4*po(ip,1)-1) -phyNodes(3,idim)*4*po(ip,0) 
						+ phyNodes(4,idim)*4*po(ip,0) +phyNodes(5,idim)*4*(1-2*po(ip,1)-po(ip,0));
				}
			}
		}
	}
	else if (shape == QUADRANGLE)
	{
		//TODO: Add bilinear and biquadratic shape functions
	}

	for(int ip = 0; ip < po.rows(); ip++) {
		jacod[ip] = jac[ip](0,0)*jac[ip](1,1) - jac[ip](0,1)*jac[ip](1,0);
		jacoi[ip](0,0) = jac[ip](1,1)/jacod[ip]; jacoi[ip](0,1) = -jac[ip](0,1)/jacod[ip];
		jacoi[ip](1,0) = -jac[ip](1,0)/jacod[ip]; jacoi[ip](1,1) = jac[ip](0,0)/jacod[ip];
	}
}

/** Mappings upto P2 are implemented.
 */
void LagrangeMapping2D::computeAll()
{
	const Array2d<a_real>& points = quadrature->points();
	shape = quadrature->getShape();
	int npoin = points.rows();
	jaco.resize(npoin);
	jacoinv.resize(npoin);
	jacodet.resize(npoin);
	mapping.resize(npoin,NDIM);

	for(int ip = 0; ip < npoin, ip++) {
		jaco[ip].resize(NDIM,NDIM);
		jacoinv[ip].resize(NDIM,NDIM);
	}

	computeJacobians(points, jaco, jacoinv, jacodet);

	/*if (shape == TRIANGLE)
	{
		if(degree == 1) {
			for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) 
				{
					mapping(ip,idim) = phyNodes(0,idim)*(1.0-points(ip,0)-points(ip,1)) + phyNodes(1,idim)*points(ip,0) + phyNodes(2,idim)*points(ip,1);
				}
			}
		}
		else if(degree == 2) {
			for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) 
				{
					mapping(ip,idim) = phyNodes(0,idim) * (1.0-3*points(ip,0)-3*points(ip,1)+2*points(ip,0)*points(ip,0)+2*points(ip,1)*points(ip,1)+4*points(ip,0)*points(ip,1)) 
						+ phyNodes(1,idim)*(2.0*points(ip,0)*points(ip,0)-points(ip,0)) 
						+ phyNodes(2,idim)*(2.0*points(ip,1)*points(ip,1)-points(ip,1)) 
						+ phyNodes(3,idim)*4.0*(points(ip,0)-points(ip,0)*points(ip,0)-points(ip,0)*points(ip,1)) 
						+ phyNodes(4,idim)*4.0*points(ip,0)*points(ip,1) 
						+ phyNodes(5,idim)*4.0*(points(ip,1)-points(ip,1)*points(ip,1)-points(ip,0)*points(ip,1));
				}
			}
		}
	}
	else if (shape == QUADRANGLE)
	{
		//TODO: Add bilinear and biquadratic shape functions
	}*/
}

void LagrangeMapping2D::computeMappingAndJacobianDet()
{
	const Array2d<a_real>& points = quadrature->points();
	shape = quadrature->getShape();
	int npoin = points.rows();
	jacodet.resize(npoin);
	mapping.resize(npoin,NDIM);

	Matrix jacol; jacol.resize(NDIM,NDIM);

	if(shape == TRIANGLE)
	{
		if(degree == 1) {
			for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) 
				{
					mapping(ip,idim) = phyNodes(0,idim)*(1.0-points(ip,0)-points(ip,1)) + phyNodes(1,idim)*points(ip,0) + phyNodes(2,idim)*points(ip,1);
					jacol(idim,0) = phyNodes(1,idim)-phyNodes(0,idim);
					jacol(idim,1) = phyNodes(2,idim)-phyNodes(0,idim);
				}
			}
		}
		else if(degree == 2) {
			for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) 
				{
					mapping(ip,idim) = phyNodes(0,idim) * (1.0-3*points(ip,0)-3*points(ip,1)+2*points(ip,0)*points(ip,0)+2*points(ip,1)*points(ip,1)+4*points(ip,0)*points(ip,1)) 
						+ phyNodes(1,idim)*(2.0*points(ip,0)*points(ip,0)-points(ip,0)) 
						+ phyNodes(2,idim)*(2.0*points(ip,1)*points(ip,1)-points(ip,1)) 
						+ phyNodes(3,idim)*4.0*(points(ip,0)-points(ip,0)*points(ip,0)-points(ip,0)*points(ip,1)) 
						+ phyNodes(4,idim)*4.0*points(ip,0)*points(ip,1) 
						+ phyNodes(5,idim)*4.0*(points(ip,1)-points(ip,1)*points(ip,1)-points(ip,0)*points(ip,1));
					jacol(idim,0) = phyNodes(0,idim)*(-3+4*points(ip,0)+4*points(ip,1)) +phyNodes(1,idim)*(4*points(ip,0)-1) +phyNodes(3,idim)*4*(1-2*points(ip,0)-points(ip,1)) 
						+phyNodes(4,idim)*4*points(ip,1) -phyNodes(5,idim)*4.0*points(ip,1);
					jacol(idim,1) = phyNodes(0,idim)*(-3+4*points(ip,1)+4*points(ip,0)) +phyNodes(2,idim)*(4*points(ip,1)-1) -phyNodes(3,idim)*4*points(ip,0) 
						+ phyNodes(4,idim)*4*points(ip,0) +phyNodes(5,idim)*4*(1-2*points(ip,1)-points(ip,0));
				}
			}
		}
	}
	else if (shape == QUADRANGLE)
	{
		//TODO: Add bilinear and biquadratic shape functions
	}

	for(int ip = 0; ip < npoin; ip++)
		jacodet[ip] = jacol(0,0)*jacol(1,1) - jacol(0,1)*jacol(1,0);
}

/** We currently have upto P2 elements.
 * The number of DOFs is computed as \f$ \sum_{i=1}^{p+1} i \f$ for p = 0,1,2...
 * For computing the element centers and area, we use the all quadrature points of the quadrature object in the gmap.
 * Note that for a quad element of degree bi p (p=1 is bi linear etc), the jacodet is of degree bi 2p-1.
 * For a tri element of degree p, the jacodet is of degree 2p-2.
 */
void TaylorElement::initialize(int degr, const GeomMapping2D* geommap)
{
	type = PHYSICAL;
	degree = degr; gmap = geommap;
	gmap->computeForPhysicalElement();

	ndof = 0;
	for(int i = 1; i <= degree+1; i++)
		ndof += i;

	int ngauss = gmap->getQuadrature()->numGauss();
	basis.resize(ngauss,ndof);
	basisGrad.resize(ngauss);
	for(int i = 0; i < ngauss; i++) {
		basisGrad[i].resize(ndof,NDIM);
	}
	
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
	
	const Array2d<a_real>& phn = gmap->getPhyNodes();
	for(int i = 0; i < phn.rows(); i++) {
		for(int j = i+1; j < phn.rows(); j++) 
		{
			a_real dist[NDIM];
			for(int idim = 0; idim < NDIM; idim++) 
			{
				dist[idim] = fabs(phn(i,idim)-phn(j,idim));
				if(dist[idim] > delta[idim])
					delta[idim] = dist[idim];
			}
		}
	}
	
	// Compute element centers and basis offsets

	const Array2d<a_real>& gp = gmap->getQuadrature()->points();
	const Array2d<a_real>& gw = gmap->getQuadrature()->weights();
	int ng = gp.rows();
#ifdef DEBUG
	if(ng != ngauss)
		std::cout << "! TaylorElement: initialize(): ngauss does not match size of gauss point array!" << std::endl;
#endif
	for(int ig = 0; ig < ng; ig++) 
	{
		area += gmap->jacDet()[ig] * gw(ig);
		for(int idim = 0; idim < NDIM; idim++)
			center[idim] += gmap->map()(ig,idim) * gmap->jacDet()[ig] * gw(ig);
	}
	for(int idim = 0; idim < NDIM; idim++)
		center[idim] /= area;

	if(degree >= 2) {
		for(int ig = 0; ig < ng; ig++)
		{
			basisOffset[0][0] += (gmap->map()(ig,0)-center[0])*(gmap->map()(ig,0)-center[0]) * gmap->jacDet(ig) * gw(ig);
			basisOffset[0][1] += (gmap->map()(ig,1)-center[1])*(gmap->map()(ig,1)-center[1]) * gmap->jacDet(ig) * gw(ig);
			basisOffset[0][2] += (gmap->map()(ig,0)-center[0])*(gmap->map()(ig,1)-center[1]) * gmap->jacDet(ig) * gw(ig);
		}
		basisOffset[0][0] *= 1.0/(area*2*delta[0]*delta[0]);
		basisOffset[0][1] *= 1.0/(area*2*delta[1]*delta[1]);
		basisOffset[0][2] *= 1.0/(area*delta[0]*delta[1]);
	}
	
	// Compute basis functions and gradients
	const Array2d<a_real>& gp = gmap->map();
	computeBasis(gp, basis);
	computeBasisGrads(gp, basisGrad);
}

void TaylorElement::computeBasis(const amat::Array2d<a_real>& __restrict__ gp, amat::Array2d<a_real>& __restrict__ basiss)
{
	for(int ip = 0; ip < gp.rows(); ip++) 
		basis(ip,0) = 1.0;

	if(degree >= 1) {
		for(int ip = 0; ip < gp.rows(); ip++) 
		{
			basiss(ip,1) = (gp(ip,0)-center[0])/delta[0];
			basiss(ip,2) = (gp(ip,1)-center[1])/delta[1];
		}
	}

	if(degree >= 2) {
		for(int ip = 0; ip < gp.rows(); ip++) 
		{
			basiss(ip,3) = (gp(ip,0)-center[0])*(gp(ip,0)-center[0])/(2.0*delta[0]*delta[0]) - basisOffset[0][0];
			basiss(ip,4) = (gp(ip,1)-center[1])*(gp(ip,1)-center[1])/(2.0*delta[1]*delta[1]) - basisOffset[0][1];
			basiss(ip,5) = (gp(ip,0)-center[0])*(gp(ip,1)-center[1])/(delta[0]*delta[1]) - basisOffset[0][2];
		}
	}
}

void computeBasisGrads(const amat::Array2d<a_real>& __restrict__ gp, std::vector<Matrix>& __restrict__ basisG) const
{
	for(int ip = 0; ip < gp.rows(); ip++)
		basisGrad[ip](0,0) = basisGrad[ip](0,1) = 0.0;

	if(degree >= 1) {
		for(int ip = 0; ip < gp.rows(); ip++) 
		{
			basisG[ip](1,0) = 1.0/delta[0]; basisG[ip](1,1) = 0.0;
			basisG[ip](2,0) = 0.0;          basisG[ip](2,1) = 1.0/delta[1];
		}
	}

	if(degree >= 2) {
		for(int ip = 0; ip < gp.rows(); ip++) 
		{
			basisG[ip](3,0) = (gp(ip,0)-center[0])/(delta[0]*delta[0]); basisG[ip](3,1) = 0.0;
			basisG[ip](4,0) = 0.0;                                      basisG[ip](4,1) = (gp(ip,1)-center[1])/(delta[1]*delta[1]);
			basisG[ip](5,0) = (gp(ip,1)-center[1])/(delta[0]*delta[1]); basisG[ip](5,1) = (gp(ip,0)-center[0])/(delta[0]*delta[1]);
		}
	}
}

void LagrangeElement::initialize(int degr, const GeomMapping2D* geommap)
{
	type = REFERENCE;
	degree = degr; gmap = geommap;
	gmap->computeForReferenceElement();

	if(gmap->getShape() == QUADRANGLE)
		ndof = (degree+1)*(degree+1);
	else {
		ndof = 0;
		for(int i = 1; i <= degree+1; i++)
			ndof += i;
	}	

	int ngauss = gmap->getQuadrature()->numGauss();
	basis.setup(ngauss,ndof);
	basisGrad.resize(ngauss);
	for(int i = 0; i < ngauss; i++) {
		basisGrad[i].resize(ndof,NDIM);
	}
	
	// Compute basis functions and gradients
	const Array2d<a_real>& gp = gmap->getQuadrature()->points();
	computeBasis(gp, basis);
	computeBasisGrads(gp, basisGrad);
}

void LagrangeElement::computeBasis(const Array2d<a_real>& __restrict__ gp, Array2d<a_real>& __restrict__ basisv) const
{
	if(gmap->getShape() == TRIANGLE) {
		if(degree == 1) {
			for(int ip = 0; ip < gp.rows(); ip++)
			{
				basisv(ip,0) = (1.0-gp(ip,0)-gp(ip,1));
				basisv(ip,1) = gp(ip,0);
				basisv(ip,2) = gp(ip,1);
			}
		}
		if(degree == 2) {
			for(int ip = 0; ip < gp.rows(); ip++)
			{
				basisv(ip,0) = 1.0-3*gp(ip,0)-3*gp(ip,1)+2*gp(ip,0)*gp(ip,0)+2*gp(ip,1)*gp(ip,1)+4*gp(ip,0)*gp(ip,1);
				basisv(ip,1) = 2.0*gp(ip,0)*gp(ip,0)-gp(ip,0);
				basisv(ip,2) = 2.0*gp(ip,1)*gp(ip,1)-gp(ip,1);
				basisv(ip,3) = 4.0*(gp(ip,0)-gp(ip,0)*gp(ip,0)-gp(ip,0)*gp(ip,1));
				basisv(ip,4) = 4.0*gp(ip,0)*gp(ip,1);
				basisv(ip,5) = 4.0*(gp(ip,1)-gp(ip,1)*gp(ip,1)-gp(ip,0)*gp(ip,1));
			}
		}
	}
	else {
		//TODO: Add quad lagrange basis
	}
}

void LagrangeElement::computeBasisGrads(const Array2d<a_real>& __restrict__ gp, std::vector<Matrix>& __restrict__ basisG) const
{
	// get ref coords of quadrature points and Jacobians
	const std::vector<Matrix>& __restrict__ jinv = gmap->jacInv();

	if(gmap->getShape() == TRIANGLE) {
		if(degree == 1) {
			for(int ip = 0; ip < ngauss; ip++)
			{
				basisG[ip](0,0) = -1.0; basisG[ip](0,1) = -1.0;
				basisG[ip](1,0) = 1.0;  basisG[ip](1,1) = 0.0;
				basisG[ip](2,0) = 0.0;  basisG[ip](2,1) = 1.0;
			}
		}
		if(degree == 2) {
			for(int ip = 0; ip < ngauss; ip++)
			{
				basisG[ip](0,0) = -3.0+4*gp(ip,0)+4*gp(ip,1);    basisG[ip](0,1) = -3.0+4*gp(ip,0)+4*gp(ip,1);
				basisG[ip](1,0) = 4.0*gp(ip,0)-1.0;              basisG[ip](1,1) = 0;
				basisG[ip](2,0) = 0;                             basisG[ip](2,1) = 4.0*gp(ip,1)-1.0;
				basisG[ip](3,0) = 4.0*(1.0-2*gp(ip,0)-gp(ip,1)); basisG[ip](3,1) = -4.0*gp(ip,0);
				basisG[ip](4,0) = 4.0*gp(ip,1);                  basisG[ip](4,1) = 4.0*gp(ip,0);
				basisG[ip](5,0) = -4.0*gp(ip,1);                 basisG[ip](5,1) = 4.0*(1.0-2*gp(ip,1)-gp(ip,0));
			}
		}	
	}
	else {
		//TODO: Add quad lagrange basis
	}
	
	for(int ip = 0; ip < ngauss; ip++)
	{
		/** To compute gradients in physical space, we use the following.
		 * Let \f$ a := \nabla_x B(x(\xi)) \f$ and \f$ b = \nabla_\xi B(x(\xi)) \f$. Then,
		 * we need \f$ a = J^{-T} b \f$. Instead, we can compute \f$ a^T = b^T J^{-1} \f$,
		 * for efficiency reasons since we have a row-major storage. This latter equation is used.
		 */
		basisG[ip] = basisG[ip]*jinv[ip];
	}
}

/** If the elements' basis functions are defined in physical space, we just compute the physical coordinates of the face quadrature points,
 * and use the physical coordinates to compute basis function values.
 * However, if the elements' basis functions are defined in reference space, we need to compute reference coordinates of the face quadrature points
 * with respect to the elements. This is done as described below.
 *
 * \verbatim
 *   |\              3
 *   | \           |------|
 *  3|   \ 2       |      |
 *   |    \      4 |      |2
 *   |_____\       |______|
 *      1             1                                            
 * \endverbatim
 *
 * Triangle
 * --------
 * The labels indicate local face numbers used (plus 1, as a zero-base is actually used).
 * Suppose \f$ \zeta \in [-1,1] \f$ is the face reference coordinate.
 * Face 1: \f[ \xi = \frac12 (1+\zeta), \, \eta = 0 \f]
 * Face 2: \f[ \xi = \frac12 (1-\zeta), \, \eta = \frac12 (1+\zeta) \f]
 * Face 3: \f[ \xi = 0, \, \eta = \frac12 (1-\zeta) \f]
 *
 * Square
 * ------
 * Face 1: \f[ \xi = \zeta, \, \eta = -1 \f]
 * Face 2: \f[ \xi = 1, \, \eta = \zeta \f]
 * Face 3: \f[ \xi = -\zeta, \, \eta = 1 \f]
 * Face 4: \f[ \xi = -1, \, \eta = -\zeta \f]
 */
void FaceElement::initialize(const Element* lelem, const Element* relem, const GeomMapping1D* gmapping, const int l_lfn, const int r_lfn)
{
	gmap = gmapping; leftel = lelem; rightel = relem; llfn = l_lfn; rlfn = r_lfn;

	int ng = gmap->getQuadrature()->numGauss();

	leftbasis.resize(ng,lelem->getNumDOFs()); rightbasis.resize(ng,relem->getNumDOFs());

	if(lelem->getType() == PHYSICAL) {
		gmap->computeAll();
		const Array2d<a_real>& points = gmap->map();
		lelem->computeBasis(points, leftbasis);
		relem->computeBasis(points, rightbasis);
	}
	else 
	{
		// compute element reference coordinates of quadrature points from their face reference coordinates
		Array2d<a_real> lpoints(ng,NDIM), rpoints(ng,NDIM);
		const Array2d<a_real>& facepoints = gmap->getQuadrature()->points();
		
		if(lelem->getGeometricMapping()->getShape() == TRIANGLE) {
			// TODO: Setup 1d geometric map from reference interval to face of reference element
			if(llfn == 0)
				for(int ig = 0; ig < ng; ig++) {
					lpoints(ig,0) = 0.5*(1.0 + facepoints(ig));
					lpoints(ig,1) = 0;
				}
			else if(llfn == 1)
				for(int ig = 0; ig < ng; ig++) {
					lpoints(ig,0) = 0.5*(1.0 - facepoints(ig));
					lpoints(ig,1) = 0.5*(1.0 + facepoints(ig));
				}
			else
				for(int ig = 0; ig < ng; ig++) {
					lpoints(ig,0) = 0;
					lpoints(ig,1) = 0.5*(1.0 - facepoints(ig));
				}
		}
		else if(lelem->getGeometricMapping()->getShape() == QUADRANGLE) {
			if(llfn == 0)
				for(int ig = 0; ig < ng; ig++) {
					lpoints(ig,0) = facepoints(ig);
					lpoints(ig,1) = -1.0;
				}
			else if(llfn == 1)
				for(int ig = 0; ig < ng; ig++) {
					lpoints(ig,0) = 1.0;
					lpoints(ig,1) = facepoints(ig);
				}
			else if(llfn == 2)
				for(int ig = 0; ig < ng; ig++) {
					lpoints(ig,0) = -facepoints(ig);
					lpoints(ig,1) = 1.0;
				}
			else
				for(int ig = 0; ig < ng; ig++) {
					lpoints(ig,0) = -1.0;
					lpoints(ig,1) = -facepoints(ig);
				}
		}
		
		if(relem->getGeometricMapping()->getShape() == TRIANGLE) {
			if(rlfn == 0)
				for(int ig = 0; ig < ng; ig++) {
					rpoints(ig,0) = 0.5*(1.0 + facepoints(ig));
					rpoints(ig,1) = 0;
				}
			else if(rlfn == 1)
				for(int ig = 0; ig < ng; ig++) {
					rpoints(ig,0) = 0.5*(1.0 - facepoints(ig));
					rpoints(ig,1) = 0.5*(1.0 + facepoints(ig));
				}
			else
				for(int ig = 0; ig < ng; ig++) {
					rpoints(ig,0) = 0;
					rpoints(ig,1) = 0.5*(1.0 - facepoints(ig));
				}
		}
		else if(relem->getGeometricMapping()->getShape() == QUADRANGLE) {
			if(rlfn == 0)
				for(int ig = 0; ig < ng; ig++) {
					rpoints(ig,0) = facepoints(ig);
					rpoints(ig,1) = -1.0;
				}
			else if(rlfn == 1)
				for(int ig = 0; ig < ng; ig++) {
					rpoints(ig,0) = 1.0;
					rpoints(ig,1) = facepoints(ig);
				}
			else if(rlfn == 2)
				for(int ig = 0; ig < ng; ig++) {
					rpoints(ig,0) = -facepoints(ig);
					rpoints(ig,1) = 1.0;
				}
			else
				for(int ig = 0; ig < ng; ig++) {
					rpoints(ig,0) = -1.0;
					rpoints(ig,1) = -facepoints(ig);
				}
		}
	
		// now compute basis function values
		lelem->computeBasis(lpoints, leftbasis);
		relem->computeBasis(rpoints, rightbasis);
	}
}

}
