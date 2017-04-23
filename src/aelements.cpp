/** @file elements.cpp
 * @brief Implementation of finite elements.
 * @author Aditya Kashi
 * @date 2017 March 1
 */

#include "aelements.hpp"

namespace acfd {

using namespace amat;

/** Computes Lagrange basis function values at given points in the reference element.
 * \note NOTE: For efficiency, we would want to able to request computation of only certain basis functions.
 */
void getLagrangeBasis(const Matrix& __restrict__ gp, const Shape shape, const int degree,  Matrix& __restrict__ basisv)
{
	if(shape == TRIANGLE) {
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
				basisv(ip,0) = 1.0 - 3*gp(ip,0) - 3*gp(ip,1) + 2*gp(ip,0)*gp(ip,0) + 2*gp(ip,1)*gp(ip,1) + 4*gp(ip,0)*gp(ip,1);
				basisv(ip,1) = 2.0*gp(ip,0)*gp(ip,0) - gp(ip,0);
				basisv(ip,2) = 2.0*gp(ip,1)*gp(ip,1) - gp(ip,1);
				basisv(ip,3) = 4.0*(gp(ip,0) - gp(ip,0)*gp(ip,0) - gp(ip,0)*gp(ip,1));
				basisv(ip,4) = 4.0*gp(ip,0)*gp(ip,1);
				basisv(ip,5) = 4.0*(gp(ip,1) - gp(ip,1)*gp(ip,1) - gp(ip,0)*gp(ip,1));
			}
		}
	}
	else {
		if(degree == 1) {
			for(int ip = 0; ip < gp.rows(); ip++)
			{
				basisv(ip,0) = (1-gp(ip,0))*(1-gp(ip,1))*0.25;
				basisv(ip,1) = (1+gp(ip,0))*(1-gp(ip,1))*0.25;
				basisv(ip,2) = (1+gp(ip,0))*(1+gp(ip,1))*0.25;
				basisv(ip,3) = (1-gp(ip,0))*(1+gp(ip,1))*0.25;
			}
		}
		if(degree == 2) {
			for(int ip = 0; ip < gp.rows(); ip++)
			{
				basisv(ip,0) = (gp(ip,0)*gp(ip,0)-gp(ip,0)) * (gp(ip,1)*gp(ip,1)-gp(ip,1)) * 0.25;
				basisv(ip,1) = (gp(ip,0)*gp(ip,0)+gp(ip,0)) * (gp(ip,1)*gp(ip,1)-gp(ip,1)) * 0.25;
				basisv(ip,2) = (gp(ip,0)*gp(ip,0)+gp(ip,0)) * (gp(ip,1)*gp(ip,1)+gp(ip,1)) * 0.25;
				basisv(ip,3) = (gp(ip,0)*gp(ip,0)-gp(ip,0)) * (gp(ip,1)*gp(ip,1)+gp(ip,1)) * 0.25;
				basisv(ip,4) = (1-gp(ip,0)*gp(ip,0))        * (gp(ip,1)*gp(ip,1)-gp(ip,1)) * 0.5;
				basisv(ip,5) = (gp(ip,0)*gp(ip,0)+gp(ip,0)) * (1-gp(ip,1)*gp(ip,1))        * 0.5;
				basisv(ip,6) = (1-gp(ip,0)*gp(ip,0))        * (gp(ip,1)*gp(ip,1)+gp(ip,1)) * 0.5;
				basisv(ip,7) = (gp(ip,0)*gp(ip,0)-gp(ip,0)) * (1-gp(ip,1)*gp(ip,1))        * 0.5;
				basisv(ip,8) = (1-gp(ip,0)*gp(ip,0))        * (1-gp(ip,1)*gp(ip,1));
			}
		}
	}
}

/** Computes Lagrange basis function gradients (w.r.t. reference coords) at given points in the reference element.
 * \note NOTE: For efficiency, we would want to able to request computation of gradients of only certain basis functions.
 */
void getLagrangeBasisGrads(const Matrix& __restrict__ gp, const Shape shape, const int degree, std::vector<Matrix>& __restrict__ basisG)
{
	if(shape == TRIANGLE) {
		if(degree == 1) {
			for(int ip = 0; ip < gp.rows(); ip++)
			{
				basisG[ip](0,0) = -1.0; basisG[ip](0,1) = -1.0;
				basisG[ip](1,0) = 1.0;  basisG[ip](1,1) = 0.0;
				basisG[ip](2,0) = 0.0;  basisG[ip](2,1) = 1.0;
			}
		}
		if(degree == 2) {
			for(int ip = 0; ip < gp.rows(); ip++)
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
	else { // QUADRANGLE
		if(degree == 1) {
			for(int ip = 0; ip < gp.rows(); ip++)
			{
				basisG[ip](0,0) = (gp(ip,1)-1)*0.25;    basisG[ip](0,1) = (gp(ip,0)-1)*0.25;
				basisG[ip](1,0) = (1-gp(ip,1))*0.25;    basisG[ip](1,1) = -(1+gp(ip,0))*0.25;
				basisG[ip](2,0) = (1+gp(ip,1))*0.25;    basisG[ip](2,1) = (1+gp(ip,0))*0.25;
				basisG[ip](3,0) = -(1+gp(ip,1))*0.25;   basisG[ip](3,1) = (1-gp(ip,0))*0.25;
			}
		}
		if(degree == 2) {
			for(int ip = 0; ip < gp.rows(); ip++)
			{
				basisG[ip](0,0) = (2*gp(ip,0)-1)*gp(ip,1)*(gp(ip,1)-1)*0.25;  basisG[ip](0,1) = gp(ip,0)*(gp(ip,0)-1)*(2*gp(ip,1)-1)*0.25;
				basisG[ip](1,0) = (2*gp(ip,0)+1)*gp(ip,1)*(gp(ip,1)-1)*0.25;  basisG[ip](1,1) = gp(ip,0)*(gp(ip,0)+1)*(2*gp(ip,1)-1)*0.25;
				basisG[ip](2,0) = (2*gp(ip,0)+1)*gp(ip,1)*(gp(ip,1)+1)*0.25;  basisG[ip](2,1) = gp(ip,0)*(gp(ip,0)+1)*(2*gp(ip,1)+1)*0.25;
				basisG[ip](3,0) = (2*gp(ip,0)-1)*gp(ip,1)*(gp(ip,1)+1)*0.25;  basisG[ip](3,1) = gp(ip,0)*(gp(ip,0)-1)*(2*gp(ip,1)+1)*0.25;
				basisG[ip](4,0) = -gp(ip,0)*gp(ip,1)*(gp(ip,1)-1);            basisG[ip](4,1) = (1-gp(ip,0)*gp(ip,0))*(2*gp(ip,1)-1)*0.5;
				basisG[ip](5,0) = (2*gp(ip,0)+1)*(1-gp(ip,1)*gp(ip,1))*0.5;   basisG[ip](5,1) = -(gp(ip,0)*gp(ip,0)+gp(ip,0))*gp(ip,1);
				basisG[ip](6,0) = -gp(ip,0)*(gp(ip,1)*gp(ip,1)+gp(ip,1));     basisG[ip](6,1) = (1-gp(ip,0)*gp(ip,0))*(2*gp(ip,1)+1)*0.5;
				basisG[ip](7,0) = (2*gp(ip,0)-1)*(1-gp(ip,1)*gp(ip,1))*0.5;   basisG[ip](7,1) = (gp(ip,0)*gp(ip,0)-gp(ip,0))*(-gp(ip,1));
				basisG[ip](8,0) = -2*gp(ip,0)*(1-gp(ip,1)*gp(ip,1));          basisG[ip](8,1) = -2*gp(ip,1)*(1-gp(ip,0)*gp(ip,0));
			}
		}	
	}
}

/// Evaluate the values of functions with a BLAS 3 call
/** \param[in] dofs The matrix of DOFs - each row contains the DOFs of one physical variable (nvars x ndofs)
 * \param[in] basisv Matrix of basis function values at the points at which function values are needed; 
 * each row has all basis function values corresponding to a given point in space (npoints x ndofs)
 * \param[in|out] interp Pre-allocated matrix for storing interpolated values (npoints x nvars)
 */
inline void evaluateFunctions(const Matrix& __restrict__ dofs, const Matrix& __restrict__ basisv, Matrix& __restrict__ interp)
{
	interp.noalias() = basisv * dofs.transpose();
}

/// Evaluate values of gradients of a function using BLAS 3 calls for all points at which evaluation is needed
/** \param[in] dofs The matrix of DOFs - each row contains the DOFs of one physical variable (nvars x ndofs)
 * \param[in] basisg Set of matrices of basis gradient values at the points at which function values are needed; 
 * each row has the gradient corresponding to a given dof (npoints x (ndofs x ndim))
 * \param[in|out] interp Pre-allocated set of matrices for storing interpolated values (npoints x (nvars x ndim))
 */
template <typename T>
inline void evaluateGradients(const Matrix& __restrict__ dofs, const std::vector<Matrix>& __restrict__ basisg, std::vector<T>& __restrict__ interp)
{
	for(size_t ig = 0; ig < basisg.size(); ig++)
		interp[ig].noalias() = dofs * basisg[ig];
}

/// A global function for computing 2D Lagrange mapping derivatives
/** Mappings upto P2 are implemented.
 */
void getLagrangeJacobianDetAndInverse(const Matrix& __restrict__ po, const Shape shape, const int degree, const Matrix& __restrict__ phy,
		std::vector<MatrixDim>& __restrict__ jacoi, std::vector<a_real>& __restrict__ jacod)
{
	int np = po.rows();
	std::vector<MatrixDim> jac(np);

	if (shape == TRIANGLE)
	{
		if(degree == 1) {
			std::vector<Matrix> bgrad(np);
			for(int i = 0; i < np; i++)
				bgrad[i].resize(3,NDIM);
			getLagrangeBasisGrads(po, shape, degree, bgrad);
			evaluateGradients<MatrixDim>(phy, bgrad, jac);

			/*std::vector<MatrixDim> jact(np);
			for(int ip = 0; ip < np; ip++) {
				for(int idim = 0; idim < NDIM; idim++) 
				{
					jact[ip](idim,0) = phy(idim,1)-phy(idim,0);
					jact[ip](idim,1) = phy(idim,2)-phy(idim,0);
				}
				std::cout << jac[ip] << "\n\n";
				std::cout << jact[ip] << "\n**\n";
			}*/
		}
		else if(degree == 2) {
			std::vector<Matrix> bgrad(np);
			for(int i = 0; i < np; i++)
				bgrad[i].resize(6,NDIM);
			getLagrangeBasisGrads(po, shape, degree, bgrad);
			evaluateGradients<MatrixDim>(phy, bgrad, jac);
			
			/*for(int ip = 0; ip < np; ip++) {
				for(int idim = 0; idim < NDIM; idim++) 
				{
					jac[ip](idim,0) = phy(0,idim)*(-3+4*po(ip,0)+4*po(ip,1)) +phy(1,idim)*(4*po(ip,0)-1) +phy(3,idim)*4*(1-2*po(ip,0)-po(ip,1)) 
						+phy(4,idim)*4*po(ip,1) -phy(5,idim)*4.0*po(ip,1);
					jac[ip](idim,1) = phy(0,idim)*(-3+4*po(ip,1)+4*po(ip,0)) +phy(2,idim)*(4*po(ip,1)-1) -phy(3,idim)*4*po(ip,0) 
						+ phy(4,idim)*4*po(ip,0) +phy(5,idim)*4*(1-2*po(ip,1)-po(ip,0));
				}
			}*/
		}
	}
	else if (shape == QUADRANGLE)
	{
		if(degree == 1) {
			std::vector<Matrix> bgrad(np);
			for(int i = 0; i < np; i++)
				bgrad[i].resize(4,NDIM);
			getLagrangeBasisGrads(po, shape, degree, bgrad);
			evaluateGradients<MatrixDim>(phy, bgrad, jac);
			/*for(int ip = 0; ip < np; ip++) {
				for(int idim = 0; idim < NDIM; idim++) 
				{
					jac[ip](idim,0) = ( phy(0,idim)*(po(ip,1)-1) + phy(1,idim)*(1-po(ip,1))
						              + phy(2,idim)*(1+po(ip,1)) + phy(3,idim)*(-1-po(ip,1)) ) * 0.25;
					jac[ip](idim,1) = ( phy(0,idim)*(po(ip,0)-1) + phy(1,idim)*(-1-po(ip,0))
				                      + phy(2,idim)*(1+po(ip,0)) + phy(3,idim)*(1-po(ip,0)) ) * 0.25;
				}
			}*/
		}
		else if(degree == 2) {
			std::vector<Matrix> bgrad(np);
			for(int i = 0; i < np; i++)
				bgrad[i].resize(9,NDIM);
			getLagrangeBasisGrads(po, shape, degree, bgrad);
			evaluateGradients<MatrixDim>(phy, bgrad, jac);
		}
	}

	for(int ip = 0; ip < np; ip++) {
		jacod[ip] = jac[ip](0,0)*jac[ip](1,1) - jac[ip](0,1)*jac[ip](1,0);
		jacoi[ip](0,0) = jac[ip](1,1)/jacod[ip]; jacoi[ip](0,1) = -jac[ip](0,1)/jacod[ip];
		jacoi[ip](1,0) = -jac[ip](1,0)/jacod[ip]; jacoi[ip](1,1) = jac[ip](0,0)/jacod[ip];
	}
}

/// Value of 2D Lagrange mapping at any reference coordinates
void getLagrangeMap(const Matrix& __restrict__ pts, const Shape shape, 
		const int degree, const Matrix& __restrict__ phy, Matrix& __restrict__ mapping)
{
	if(shape == TRIANGLE)
	{
		if(degree == 1) {
			Matrix bas(pts.rows(), 3);
			getLagrangeBasis(pts, shape, degree, bas);
			evaluateFunctions(phy, bas, mapping);
			/*for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) 
				{
					mapping(ip,idim) = phy(0,idim)*(1.0-pts(ip,0)-pts(ip,1)) + phy(1,idim)*pts(ip,0) + phy(2,idim)*pts(ip,1);
				}
			}*/
		}
		else if(degree == 2) {
			Matrix bas(pts.rows(), 6);
			getLagrangeBasis(pts, shape, degree, bas);
			evaluateFunctions(phy, bas, mapping);
			/*for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) 
				{
					mapping(ip,idim) = phy(0,idim) * (1.0-3*pts(ip,0)-3*pts(ip,1)+2*pts(ip,0)*pts(ip,0)+2*pts(ip,1)*pts(ip,1)+4*pts(ip,0)*pts(ip,1)) 
						+ phy(1,idim)*(2.0*pts(ip,0)*pts(ip,0)-pts(ip,0)) 
						+ phy(2,idim)*(2.0*pts(ip,1)*pts(ip,1)-pts(ip,1)) 
						+ phy(3,idim)*4.0*(pts(ip,0)-pts(ip,0)*pts(ip,0)-pts(ip,0)*pts(ip,1)) 
						+ phy(4,idim)*4.0*pts(ip,0)*pts(ip,1) 
						+ phy(5,idim)*4.0*(pts(ip,1)-pts(ip,1)*pts(ip,1)-pts(ip,0)*pts(ip,1));
				}
			}*/
		}
	}
	else if (shape == QUADRANGLE)
	{
		if(degree == 1) {
			Matrix bas(pts.rows(), 4);
			getLagrangeBasis(pts, shape, degree, bas);
			evaluateFunctions(phy, bas, mapping);
			/*for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) {
					mapping(ip,idim) = ( phy(0,idim)*(1-pts(ip,0))*(1-pts(ip,1))
						               + phy(1,idim)*(1+pts(ip,0))*(1-pts(ip,1))
									   + phy(2,idim)*(1+pts(ip,0))*(1+pts(ip,1))
									   + phy(3,idim)*(1-pts(ip,0))*(1+pts(ip,1)) ) * 0.25;
				}
			}*/
		}
		else if(degree == 2) {
			Matrix bas(pts.rows(), 9);
			getLagrangeBasis(pts, shape, degree, bas);
			evaluateFunctions(phy, bas, mapping);
			/*for(int ip = 0; ip < npoin; ip++)
			{
				for(int idim = 0; idim < NDIM; idim++) {
					mapping(ip,idim) =   phy(0,idim)* (pts(ip,0)*pts(ip,0)-pts(ip,0)) * (pts(ip,1)*pts(ip,1)-pts(ip,1)) * 0.25
						               + phy(1,idim)* (pts(ip,0)*pts(ip,0)+pts(ip,0)) * (pts(ip,1)*pts(ip,1)-pts(ip,1)) * 0.25
									   + phy(2,idim)* (pts(ip,0)*pts(ip,0)+pts(ip,0)) * (pts(ip,1)*pts(ip,1)+pts(ip,1)) * 0.25
									   + phy(3,idim)* (pts(ip,0)*pts(ip,0)-pts(ip,0)) * (pts(ip,1)*pts(ip,1)+pts(ip,1)) * 0.25
									   + phy(4,idim)* (1-pts(ip,0)*pts(ip,0))   * (pts(ip,1)*pts(ip,1)-pts(ip,1)) * 0.5
									   + phy(5,idim)* (pts(ip,0)*pts(ip,0)+phi(ip,0)) * (1-pts(ip,1)*pts(ip,1))   * 0.5
									   + phy(6,idim)* (1-pts(ip,0)*pts(ip,0))   * (pts(ip,1)*pts(ip,1)+pts(ip,1)) * 0.5
									   + phy(7,idim)* (pts(ip,0)*pts(ip,0)-phi(ip,0)) * (1-pts(ip,1)*pts(ip,1))   * 0.5
									   + phy(8,idim)* (1-pts(ip,0)*pts(ip,0))*(1-pts(ip,1)*pts(ip,1));
				}
			}*/
		}
	}
}

/** Currently, Lagrange mappings upto polynomial degree 2 are implemented.
 */
void LagrangeMapping1D::computeAll()
{
	const Matrix& points = quadrature->points();
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
				// mapping - sum of Lagrange shape functions multiplied by resp coeffs
				mapping(i,idim) = phyNodes(idim,0)*(1.0-points(i))*0.5 + phyNodes(idim,1)*(1.0+points(i))*0.5;
				// get sum Lagrange derivatives multiplied by coeffs
				vel[idim] = (phyNodes(idim,1) - phyNodes(idim,0))/2.0;
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
				mapping(i,idim) = phyNodes(idim,0)*points(i)*(points(i)-1.0)/2 + phyNodes(idim,1)*points(i)*(points(i)+1.0)/2 + phyNodes(idim,2)*(1.0-points(i)*points(i));
				vel[idim] = phyNodes(idim,0)*(points(i)-0.5) + phyNodes(idim,1)*(points(i)+0.5) + phyNodes(idim,2)*(-2.0*points(i));
			}
			
			speeds[i] = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);

			// normalize tangent to get unit tangent and thus unit normal
			for(int idim = 0; idim < NDIM; idim++)
				vel[idim] /= speeds[i];
			normals[i][0] = vel[1]; normals[i][1] = -vel[0];
		}
	}
	else
		std::cout << "! LagrangeMapping1D: Chosen geometric order not available!\n";
}

void LagrangeMapping2D::calculateMap(const Matrix& __restrict__ points, Matrix& __restrict__ maps) const
{
	getLagrangeMap(points, shape, degree, phyNodes, maps);
}

void LagrangeMapping2D::calculateJacobianDetAndInverse(const Matrix& __restrict__ po, 
		std::vector<MatrixDim>& __restrict__ jacoi, std::vector<a_real>& __restrict__ jacod) const
{
	getLagrangeJacobianDetAndInverse(po, shape, degree, phyNodes, jacoi, jacod);
}

void LagrangeMapping2D::computeForReferenceElement()
{
	const Matrix& points = quadrature->points();
	shape = quadrature->getShape();
	int npoin = points.rows();
	jacoinv.resize(npoin);
	jacodet.resize(npoin);

	getLagrangeJacobianDetAndInverse(points, shape, degree, phyNodes, jacoinv, jacodet);
}

/** We can make this more efficient by not computing the Jacobian inverse below.
 * For that we need (global functions) getLagrangeJacobianDet() and get2by2Inverse() separetely.
 */
void LagrangeMapping2D::computeForPhysicalElement()
{
	const Matrix& points = quadrature->points();
	shape = quadrature->getShape();
	int npoin = points.rows();
	jacodet.resize(npoin);
	mapping.resize(npoin,NDIM);

	std::vector<MatrixDim> jacoi(npoin);
	getLagrangeMap(points, shape, degree, phyNodes, mapping);
	getLagrangeJacobianDetAndInverse(points, shape, degree, phyNodes, jacoi, jacodet);
}

void LagrangeMapping2D::computePhysicalCoordsOfDomainQuadraturePoints()
{
	const Matrix& points = quadrature->points();
	int npoin = points.rows();
	mapping.resize(npoin,NDIM);
	
	getLagrangeMap(points, shape, degree, phyNodes, mapping);
}

/** We currently have upto P2 elements.
 * The number of DOFs is computed as \f$ \sum_{i=1}^{p+1} i \f$ for p = 0,1,2...
 * For computing the element centers and area, we use the all quadrature points of the quadrature object in the gmap.
 * Note that for a quad element of degree bi p (p=1 is bi linear etc), the jacodet is of degree bi 2p-1.
 * For a tri element of degree p, the jacodet is of degree 2p-2.
 */
void TaylorElement::initialize(int degr, GeomMapping2D* geommap)
{
	type = PHYSICAL;
	degree = degr;
	geommap->computeForPhysicalElement();
	gmap = const_cast<const GeomMapping2D*>(geommap);

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
	
	const Matrix& phn = gmap->getPhyNodes();
	for(int i = 0; i < phn.cols(); i++) {
		for(int j = i+1; j < phn.cols(); j++) 
		{
			a_real dist[NDIM];
			for(int idim = 0; idim < NDIM; idim++) 
			{
				dist[idim] = fabs(phn(idim,i)-phn(idim,j));
				if(dist[idim] > delta[idim])
					delta[idim] = dist[idim];
			}
		}
	}
	
	// Compute element centers and basis offsets

	const Array2d<a_real>& gw = gmap->getQuadrature()->weights();
	int ng = gmap->getQuadrature()->numGauss();
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

#ifdef DEBUG
	// check for P1 case
	if(gmap->getDegree() == 1) {
		a_real xc, yc;
		const Matrix& nodes = gmap->getPhyNodes();
		int npo = nodes.cols();
		for(int i = 0; i < npo; i++) {
			xc += nodes(0,i); yc += nodes(1,i);
		}
		xc /= npo; yc /= npo;
		if(std::fabs(center[0]-xc) > SMALL_NUMBER)
			std::printf(" ! TaylorElement: initialize(): Difference between integrated and trivially computed x center is %f!\n", std::fabs(center[0]-xc));
		if(std::fabs(center[1]-yc) > SMALL_NUMBER)
			std::printf(" ! TaylorElement: initialize(): Difference between integrated and trivially computed y center is %f!\n", std::fabs(center[1]-yc));
	}
#endif

	if(degree >= 2) {
		for(int ig = 0; ig < ng; ig++)
		{
			basisOffset[0][0] += (gmap->map()(ig,0)-center[0])*(gmap->map()(ig,0)-center[0]) * gmap->jacDet()[ig] * gw(ig);
			basisOffset[0][1] += (gmap->map()(ig,1)-center[1])*(gmap->map()(ig,1)-center[1]) * gmap->jacDet()[ig] * gw(ig);
			basisOffset[0][2] += (gmap->map()(ig,0)-center[0])*(gmap->map()(ig,1)-center[1]) * gmap->jacDet()[ig] * gw(ig);
		}
		basisOffset[0][0] *= 1.0/(area*2*delta[0]*delta[0]);
		basisOffset[0][1] *= 1.0/(area*2*delta[1]*delta[1]);
		basisOffset[0][2] *= 1.0/(area*delta[0]*delta[1]);
	}
	
	// Compute basis functions and gradients
	const Matrix& gp = gmap->map();
	computeBasis(gp, basis);
	std::vector<MatrixDim> dummy;
	computeBasisGrads(gp, dummy, basisGrad);
}

void TaylorElement::computeBasis(const Matrix& __restrict__ gp, Matrix& __restrict__ basiss) const
{
	for(int ip = 0; ip < gp.rows(); ip++) 
		basiss(ip,0) = 1.0;

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

void TaylorElement::computeBasisGrads(const Matrix& __restrict__ gp, const std::vector<MatrixDim>& jinv, std::vector<Matrix>& __restrict__ basisG) const
{
	for(int ip = 0; ip < gp.rows(); ip++)
		basisG[ip](0,0) = basisG[ip](0,1) = 0.0;

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

void LagrangeElement::initialize(int degr, GeomMapping2D* geommap)
{
	type = REFERENTIAL;
	degree = degr;
	geommap->computeForReferenceElement();
	gmap = const_cast<const GeomMapping2D*>(geommap);

	if(gmap->getShape() == QUADRANGLE)
		ndof = (degree+1)*(degree+1);
	else {
		ndof = 0;
		for(int i = 1; i <= degree+1; i++)
			ndof += i;
	}	

	int ngauss = gmap->getQuadrature()->numGauss();
	basis.resize(ngauss,ndof);
	basisGrad.resize(ngauss);
	for(int i = 0; i < ngauss; i++) {
		basisGrad[i].resize(ndof,NDIM);
	}
	
	// Compute basis functions and gradients
	const Matrix& gp = gmap->getQuadrature()->points();
	const std::vector<MatrixDim>& jinv = gmap->jacInv();
	getLagrangeBasis(gp, gmap->getShape(), degree, basis);
	getLagrangeBasisGrads(gp, gmap->getShape(), degree, basisGrad);
	
	for(int ip = 0; ip < gp.rows(); ip++)
	{
		/** To compute gradients in physical space, we use the following.
		 * Let \f$ a := \nabla_x B(x(\xi)) \f$ and \f$ b = \nabla_\xi B(x(\xi)) \f$. Then,
		 * we need \f$ a = J^{-T} b \f$. Instead, we can compute \f$ a^T = b^T J^{-1} \f$,
		 * for efficiency reasons since we have a row-major storage. This latter equation is used.
		 */
		basisGrad[ip] = (basisGrad[ip]*jinv[ip]).eval();
	}
}

Matrix LagrangeElement::getReferenceNodes() const
{
	Matrix refs(ndof,NDIM);
	if(gmap->getShape() == TRIANGLE)
	{
		if(degree >= 1) {
			refs(0,0) = 0; refs(0,1) = 0;
			refs(1,0) = 1; refs(1,1) = 0;
			refs(2,0) = 0; refs(2,1) = 1;
		}
		if(degree == 2) {
			refs(3,0) = 0.5; refs(3,1) = 0.0;
			refs(4,0) = 0.5; refs(4,1) = 0.5;
			refs(5,0) = 0.0; refs(5,1) = 0.5;
		}
	}
	else if(gmap->getShape() == QUADRANGLE)
	{
		if(degree >= 1) {
			refs(0,0) = -1; refs(0,1) = -1;
			refs(1,0) = 1;  refs(1,1) = -1;
			refs(2,0) = 1;  refs(2,1) = 1;
			refs(3,0) = -1; refs(3,1) = 1;
		}
		if(degree == 2) {
			refs(4,0) = 0;  refs(4,1) = -1;
			refs(5,0) = 1;  refs(5,1) = 0;
			refs(6,0) = 0;  refs(6,1) = 1;
			refs(7,0) = -1; refs(7,1) = 0;
			refs(8,0) = 0;  refs(8,1) = 0;
		}
	}
	return refs;
}

void LagrangeElement::computeBasis(const Matrix& __restrict__ gp, Matrix& __restrict__ basisv) const
{
	getLagrangeBasis(gp, gmap->getShape(), degree, basisv);
}

void LagrangeElement::computeBasisGrads(const Matrix& __restrict__ gp, const std::vector<MatrixDim>& __restrict__ jinv, std::vector<Matrix>& __restrict__ basisG) const
{
	getLagrangeBasisGrads(gp, gmap->getShape(), degree, basisG);
	
	for(int ip = 0; ip < gp.rows(); ip++)
	{
		/** To compute gradients in physical space, we use the following.
		 * Let \f$ a := \nabla_x B(x(\xi)) \f$ and \f$ b = \nabla_\xi B(x(\xi)) \f$. Then,
		 * we need \f$ a = J^{-T} b \f$. Instead, we can compute \f$ a^T = b^T J^{-1} \f$,
		 * for efficiency reasons since we have a row-major storage. This latter equation is used.
		 */
		basisG[ip] = (basisG[ip]*jinv[ip]).eval();
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
void FaceElement::initialize(const Element*const lelem, const Element*const relem, const GeomMapping1D*const gmapping, const int l_lfn, const int r_lfn)
{
	gmap = gmapping; leftel = lelem; rightel = relem; llfn = l_lfn; rlfn = r_lfn;

	const int ng = gmap->getQuadrature()->numGauss();

	if(lelem->getType() == PHYSICAL) {
		leftbasis.resize(ng,lelem->getNumDOFs());
		const Matrix& points = gmap->map();
		lelem->computeBasis(points, leftbasis);
	}
	else if(lelem->getType() == REFERENTIAL)
	{
		// compute element reference coordinates of face quadrature points from their face reference coordinates
		Matrix lpoints(ng,NDIM);
		const Matrix& facepoints = gmap->getQuadrature()->points();
		getElementRefCoords(facepoints, leftel, llfn, 1, lpoints);
	
		// now compute basis function values
		leftbasis.resize(ng,lelem->getNumDOFs());
		lelem->computeBasis(lpoints, leftbasis);
	}

	if(relem->getType() == PHYSICAL) {
		rightbasis.resize(ng,relem->getNumDOFs());
		const Matrix& points = gmap->map();
		relem->computeBasis(points, rightbasis);
	}
	else if(relem->getType() == REFERENTIAL)
	{
		// compute element reference coordinates of face quadrature points from their face reference coordinates
		Matrix rpoints(ng,NDIM);
		const Matrix& facepoints = gmap->getQuadrature()->points();
		a_real fspeed = getElementRefCoords(facepoints, rightel, rlfn, -1, rpoints);

		/*std::vector<MatrixDim> jacinv(ng); std::vector<a_real> jacdet(ng);
		const std::vector<a_real> facedet = gmap->speed();
		rightel->getGeometricMapping()->calculateJacobianDetAndInverse(rpoints, jacinv, jacdet);*/
		/*for(int ig = 0; ig < ng; ig++) {
			a_real diff = std::fabs(jacdet[ig]*fspeed - facedet[ig]);
			std::cout << "  Rightel: Diff = " << diff;
			if(diff > ZERO_TOL) std::cout << "!! --- !!";
			std::cout << std::endl;
		}*/
	
		// now compute basis function values
		rightbasis.resize(ng,relem->getNumDOFs());
		relem->computeBasis(rpoints, rightbasis);

		/*for(int ig = 0; ig < ng; ig++)
			rightbasis.row(ig) *= jacdet[ig]*fspeed/facedet[ig];*/
	}
}

void FaceElement::computeBasisGrads()
{
	int ng = gmap->getQuadrature()->numGauss();

	if(leftel->getType() == PHYSICAL) {
		leftbgrad.resize(ng);
		for(int ip = 0; ip < ng; ip++)
			leftbgrad[ip].resize(leftel->getNumDOFs(), NDIM);

		const Matrix& gp = gmap->map();
		std::vector<MatrixDim> jacinv; // dummy
		leftel->computeBasisGrads(gp, jacinv, leftbgrad);
	}
	else if(leftel->getType() == REFERENTIAL)
	{
		leftbgrad.resize(ng);
		for(int ip = 0; ip < ng; ip++)
			leftbgrad[ip].resize(leftel->getNumDOFs(), NDIM);

		Matrix lpoints(ng,NDIM);
		const Matrix& facepoints = gmap->getQuadrature()->points();
		getElementRefCoords(facepoints, leftel, llfn, 1, lpoints);
		std::vector<MatrixDim> jacinv(ng); std::vector<a_real> jacdet(ng);
		leftel->getGeometricMapping()->calculateJacobianDetAndInverse(lpoints, jacinv, jacdet);
		
		leftel->computeBasisGrads(lpoints, jacinv, leftbgrad);
	}

	if(rightel->getType() == PHYSICAL) {
		rightbgrad.resize(ng);
		for(int ip = 0; ip < ng; ip++)
			rightbgrad[ip].resize(rightel->getNumDOFs(), NDIM);

		const Matrix& gp = gmap->map();
		std::vector<MatrixDim> jacinv; // dummy
		rightel->computeBasisGrads(gp, jacinv, rightbgrad);
	}
	else if(rightel->getType() == REFERENTIAL)
	{
		rightbgrad.resize(ng);
		for(int ip = 0; ip < ng; ip++)
			rightbgrad[ip].resize(rightel->getNumDOFs(), NDIM);

		Matrix rpoints(ng,NDIM);
		const Matrix& facepoints = gmap->getQuadrature()->points();
		getElementRefCoords(facepoints, rightel, rlfn, -1, rpoints);
		std::vector<MatrixDim> jacinv(ng); std::vector<a_real> jacdet(ng);
		rightel->getGeometricMapping()->calculateJacobianDetAndInverse(rpoints, jacinv, jacdet);
		
		rightel->computeBasisGrads(rpoints, jacinv, rightbgrad);
	}
}

/** Note that the order of points has to be reversed for the right element.
 * We use lr for this.
 */
a_real FaceElement::getElementRefCoords(const Matrix& __restrict__ facepoints, const Element *const __restrict__ elem, 
		const int llfn, const int lr, Matrix& __restrict__ dompoints)
{
	int ng = facepoints.rows();
#ifdef DEBUG
	if(ng != dompoints.rows())
		printf("!  FaceElement: getElementRefCoords: Size mismatch!\n");
#endif
	if(elem->getGeometricMapping()->getShape() == TRIANGLE) {
		if(llfn == 0)
		{
			for(int ig = 0; ig < ng; ig++) {
				dompoints(ig,0) = 0.5*(1.0 + lr*facepoints(ig));
				dompoints(ig,1) = 0;
			}
			return 0.5;
		}
		else if(llfn == 1)
		{
			for(int ig = 0; ig < ng; ig++) {
				dompoints(ig,0) = 0.5*(1.0 - lr*facepoints(ig));
				dompoints(ig,1) = 0.5*(1.0 + lr*facepoints(ig));
			}
			return 1/SQRT2;
		}
		else
		{
			for(int ig = 0; ig < ng; ig++) {
				dompoints(ig,0) = 0;
				dompoints(ig,1) = 0.5*(1.0 - lr*facepoints(ig));
			}
			return 0.5;
		}
	}
	else if(elem->getGeometricMapping()->getShape() == QUADRANGLE) {
		if(llfn == 0)
		{
			for(int ig = 0; ig < ng; ig++) {
				dompoints(ig,0) = lr*facepoints(ig);
				dompoints(ig,1) = -1.0;
			}
		}
		else if(llfn == 1)
		{
			for(int ig = 0; ig < ng; ig++) {
				dompoints(ig,0) = 1.0;
				dompoints(ig,1) = lr*facepoints(ig);
			}
		}
		else if(llfn == 2)
		{
			for(int ig = 0; ig < ng; ig++) {
				dompoints(ig,0) = lr*(-facepoints(ig));
				dompoints(ig,1) = 1.0;
			}
		}
		else
		{
			for(int ig = 0; ig < ng; ig++) {
				dompoints(ig,0) = -1.0;
				dompoints(ig,1) = lr*(-facepoints(ig));
			}
		}
		return 1.0;
	}
	
	else return 0;
}

}
