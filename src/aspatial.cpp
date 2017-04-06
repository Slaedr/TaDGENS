/** @file asolverbase.cpp
 * @brief Implements common functionality for solution of Euler/Navier-Stokes equations.
 * @author Aditya Kashi
 * @date 2016-02-24
 */

#include "asolverbase.hpp"

namespace acfd {

SpatialBase::SpatialBase(const UMesh2dh* mesh, const int _p_degree) : m(mesh), p_degree(_p_degree)
{
	std::cout << "SpatialBase: Setting up spaital integrator for FE polynomial degree " << p_degree << std::endl;
	
	// set quadrature strength
	int dom_quaddegree = 2*p_degree;
	if(m->degree() == 2) dom_quaddegree += 1;
	int boun_quaddegree = 2*p_degree;
	if(m->degree() == 2) boun_quaddegree += 1;

	dtquad = new Quadrature2DTriangle();
	dtquad.initialize(dom_quaddegree);
	dsquad = new Quadrature2DSquare();
	dsquad.initialize(dom_quaddegree);
	bquad = new Quadrature1D();
	bquad.initialize(boun_quaddegree);

	map2d = new LagrangeMapping2D[m->gnelem()];
	elems = new TaylorElement[m->gnelem()];

	map1d = new LagrangeMapping1D[m->gnaface()];
	faces = new FaceElement_PhysicalSpace[m->gnaface()-m->gnbface()];
	bfaces = new BFaceElement_PhysicalSpace[m->gnbface()];
}

SpatialBase::~SpatialBase()
{
	delete dtquad;
	delete dsquad;
	delete bquad;
	delete [] map2d;
	delete [] map1d;
	delete [] faces;
	delete [] bfaces;
	delete [] elems;
}

SpatialBase::computeFEData()
{
	std::cout << "SpatialBase: computeFEData(): Computing basis functions, basis gradients and mass matrices for each element" << std::endl;
	minv.resize(m.gnelem());

	// loop over elements to setup and compute maps and elements and compute mass matrices
	for(int iel = 0; iel < m->gnelem(); iel++)
	{
		amat::Array2d<a_real> phynodes(m->gnnode(iel),NDIM);
		for(int i = 0; i < m->gnnode(iel); i++)
			for(int j = 0; j < NDIM; j++)
				phynodes(i,j) = m->gcoords(m->ginpoel(iel,i),j);

		if(m->gnnode(iel) == 4 || m->gnnode(iel) == 9)
			map2d[iel].setAll(m->degree(), phynodes, dsquad);
		else
			map2d[iel].setAll(m->degree(), phynodes, dtquad);

		elems[iel].initialize(p_degree, &map2d[iel]);

		// allocate mass matrix
		minv[iel].resize(elems[iel].getNumDOFs(), elems[iel].getNumDOFs());
		minv[iel] = Matrix::Zero();

		// compute mass matrix
		const Quadrature2D* lquad = map2d[iel].getQuadrature();
		for(int ig = 0; ig < lquad->numGauss(); ig++)
		{
			for(int idof = 0; idof < elems[iel].getNumDOFs(); idof++)
				for(int jdof = 0; jdof < elems[iel].getNumDOFs(); jdof++)
					minv[iel](idof,jdof) += elems[iel].bFunc(ig)(idof)*elems[iel].bFunc(ig)(jdof)*map2d[iel].jacDet(ig);
		}
	}
	
	// loop over faces
	for(int iface = 0; iface < m->gnbface(); iface++)
	{
		int lelem = m->gintfac(iface,0);
		amat::Array2d<a_real> phynodes(m->gnnofa(iface),NDIM);
		for(int i = 0; i < m->gnnofa(iface); i++)
			for(int j = 0; j < NDIM; j++)
				phynodes(i,j) = m->gcoords(m->gintfac(iface,2+i),j);

		map1d[iface].setAll(m->degree(), phynodes, bquad);

		bfaces[iface].initialize(&elems[lelem], &map1d[iface]);
	}

	for(int iface = m->gnbface(); iface < m->gnaface(); iface++)
	{
		int lelem = m->gintfac(iface,0);
		int relem = m->gintfac(iface,1);
		amat::Array2d<a_real> phynodes(m->gnnofa(iface),NDIM);
		for(int i = 0; i < m->gnnofa(iface); i++)
			for(int j = 0; j < NDIM; j++)
				phynodes(i,j) = m->gcoords(m->gintfac(iface,2+i),j);

		map1d[iface].setAll(m->degree(), phynodes, bquad);
		map1d[iface].computeAll();

		faces[iface-m->gnbface()].initialize(&elems[lelem], &map1d[iface]);
	}
	
	std::cout << "SpatialBase: computeFEData(): Done." << std::endl;
}

/*void SpatialBase::compute_ghost_cell_coords_about_midpoint()
{
	int iface, ielem, idim, ip1, ip2;
	std::vector<a_real> midpoint(m->gndim());
	for(iface = 0; iface < m->gnbface(); iface++)
	{
		ielem = m->gintfac(iface,0);
		ip1 = m->gintfac(iface,2);
		ip2 = m->gintfac(iface,3);

		for(idim = 0; idim < m->gndim(); idim++)
		{
			midpoint[idim] = 0.5 * (m->gcoords(ip1,idim) + m->gcoords(ip2,idim));
		}

		for(idim = 0; idim < m->gndim(); idim++)
			rcg(iface,idim) = 2*midpoint[idim] - rc(ielem,idim);
	}
}

void SpatialBase::compute_ghost_cell_coords_about_face()
{
	int ied, ig, ielem;
	a_real x1, y1, x2, y2, xs, ys, xi, yi;

	for(ied = 0; ied < m->gnbface(); ied++)
	{
		ielem = m->gintfac(ied,0); //int lel = ielem;
		//jelem = m->gintfac(ied,1); //int rel = jelem;
		a_real nx = m->ggallfa(ied,0);
		a_real ny = m->ggallfa(ied,1);

		xi = rc.get(ielem,0);
		yi = rc.get(ielem,1);

		// Note: The ghost cell is a direct reflection of the boundary cell about the boundary-face
		//       It is NOT the reflection about the midpoint of the boundary-face
		x1 = m->gcoords(m->gintfac(ied,2),0);
		x2 = m->gcoords(m->gintfac(ied,3),0);
		y1 = m->gcoords(m->gintfac(ied,2),1);
		y2 = m->gcoords(m->gintfac(ied,3),1);

		if(fabs(nx)>SMALL_NUMBER && fabs(ny)>SMALL_NUMBER)		// check if nx != 0 and ny != 0
		{
			xs = ( yi-y1 - ny/nx*xi + (y2-y1)/(x2-x1)*x1 ) / ((y2-y1)/(x2-x1)-ny/nx);
			ys = ny/nx*xs + yi - ny/nx*xi;
		}
		else if(fabs(nx)<=A_SMALL_NUMBER)
		{
			xs = xi;
			ys = y1;
		}
		else
		{
			xs = x1;
			ys = yi;
		}
		rcg(ied,0) = 2*xs-xi;
		rcg(ied,1) = 2*ys-yi;
	}
}*/

EulerFlow::EulerFlow(const UMesh2dh* mesh, const int _p_degree, a_real gamma, Vector& u_inf, Vector& u_in, Vector& u_out, int boun_ids[6])
	: SpatialBase(mesh, _p_degree), uinf(u_inf), uin(u_in), uout(u_out), g(gamma),
	  slipwall_id(boun_ids[0]), inflow_id(boun_ids[1]), outflow_id(boun_ids[2]), farfield_id(boun_ids[3]), periodic_id(boun_id[4]), symmetry_id(boun_ids[5])
{
	Vector condinf = u_inf;
	uinf(0) = condinf(0);
	uinf(1) = uinf(0)*condinf(1)*std::cos(condinf(2));
	uinf(2) = uinf(0)*condinf(1)*std::sin(condinf(2));
	a_real p = uinf(0)*condinf(1)*condinf(1)/(g*condinf(3)*condinf(3));
	uinf(3) = p/(g-1.0) + 0.5*uinf(0)*condinf(1)*condinf(1);
}

void EulerFlow::compute_boundary_states(const a_real ins[NVARS], const Vector& n, int iface, a_real bs[NVARS])
{
	a_real vni = (ins[1]*n[0] + ins[2]*n[1])/ins[0];
	a_real pi = (g-1.0)*(ins[3] - 0.5*(pow(ins[1],2)+pow(ins[2],2))/ins[0]);
	a_real pinf = (g-1.0)*(uinf[3] - 0.5*(pow(uinf[1],2)+pow(uinf[2],2))/uinf[0]);
	a_real ci = sqrt(g*pi/ins[0]);
	a_real Mni = vni/ci;

	if(m->gintfacbtags(iface,0) == slipwall_id)
	{
		bs[0] = ins[0];
		bs[1] = ins[1] - 2*vni*n[0]*bs[0];
		bs[2] = ins[2] - 2*vni*n[1]*bs[0];
		bs[3] = ins[3];
	}
	
	if(m->gintfacbtags(iface,0) == periodic_id) {
		// TODO: Implement periodic boundary here //
	}

	if(m->gintfacbtags(iface,0) == freestream_id)
	{
		//if(Mni <= -1.0)
		{
			for(int i = 0; i < NVARS; i++)
				bs[i] = uinf[i];
		}
		/*else if(Mni > -1.0 && Mni < 0)
		{
			// subsonic inflow, specify rho and u according to FUN3D BCs paper
			for(i = 0; i < NVARS-1; i++)
				bs(ied,i) = uinf.get(0,i);
			bs(ied,3) = pi/(g-1.0) + 0.5*( uinf.get(0,1)*uinf.get(0,1) + uinf.get(0,2)*uinf.get(0,2) )/uinf.get(0,0);
		}
		else if(Mni >= 0 && Mni < 1.0)
		{
			// subsonic ourflow, specify p accoording FUN3D BCs paper
			for(i = 0; i < NVARS-1; i++)
				bs(ied,i) = ins.get(ied,i);
			bs(ied,3) = pinf/(g-1.0) + 0.5*( ins.get(ied,1)*ins.get(ied,1) + ins.get(ied,2)*ins.get(ied,2) )/ins.get(ied,0);
		}
		else
			for(i = 0; i < NVARS; i++)
				bs(ied,i) = ins.get(ied,i);*/
	}
}

void EulerFlow::compute_RHS()
{
	//std::cout << "Computing res ---\n";
#pragma omp parallel default(shared)
	{
#pragma omp for simd
		for(int iel = 0; iel < m->gnelem(); iel++)
		{
			for(int i = 0; i < NVARS; i++)
				residual(iel,i) = 0.0;
			integ(iel) = 0.0;
		}

		// first, set cell-centered values of boundary cells as left-side values of boundary faces
#pragma omp for
		for(a_int ied = 0; ied < m->gnbface(); ied++)
		{
			a_int ielem = m->gintfac(ied,0);
			for(int ivar = 0; ivar < NVARS; ivar++)
				uleft(ied,ivar) = u.get(ielem,ivar);
		}
	}

	if(order == 2)
	{
		// get cell average values at ghost cells using BCs
		compute_boundary_states(uleft, ug);

		rec->compute_gradients();
		lim->compute_face_values();
	}
	else
	{
		// if order is 1, set the face data same as cell-centred data for all faces

		// set both left and right states for all interior faces
#pragma omp parallel for
		for(a_int ied = m->gnbface(); ied < m->gnaface(); ied++)
		{
			a_int ielem = m->gintfac(ied,0);
			a_int jelem = m->gintfac(ied,1);
			for(int ivar = 0; ivar < NVARS; ivar++)
			{
				uleft(ied,ivar) = u.get(ielem,ivar);
				uright(ied,ivar) = u.get(jelem,ivar);
			}
		}
	}

	// set right (ghost) state for boundary faces
	compute_boundary_states(uleft,uright);

	/** Compute fluxes.
	 * The integral of the maximum magnitude of eigenvalue over each face is also computed:
	 * \f[
	 * \int_{f_i} (|v_n| + c) \mathrm{d}l
	 * \f]
	 * so that time steps can be calculated for explicit time stepping.
	 */

	std::vector<a_real> ci(m->gnaface()), vni(m->gnaface()), cj(m->gnaface()), vnj(m->gnaface());

#pragma omp parallel default(shared)
	{
#pragma omp for
		for(a_int ied = 0; ied < m->gnaface(); ied++)
		{
			a_real n[NDIM];
			n[0] = m->ggallfa(ied,0);
			n[1] = m->ggallfa(ied,1);
			a_real len = m->ggallfa(ied,2);

			const a_real* ulp = uleft.const_row_pointer(ied);
			const a_real* urp = uright.const_row_pointer(ied);
			a_real* fluxp = fluxes.row_pointer(ied);

			// compute flux
			inviflux->get_flux(ulp, urp, n, fluxp);

			// integrate over the face
			for(int ivar = 0; ivar < NVARS; ivar++)
					fluxp[ivar] *= len;

			//calculate presures from u
			a_real pi = (g-1)*(uleft.get(ied,3) - 0.5*(pow(uleft.get(ied,1),2)+pow(uleft.get(ied,2),2))/uleft.get(ied,0));
			a_real pj = (g-1)*(uright.get(ied,3) - 0.5*(pow(uright.get(ied,1),2)+pow(uright.get(ied,2),2))/uright.get(ied,0));
			//calculate speeds of sound
			ci[ied] = sqrt(g*pi/uleft.get(ied,0));
			cj[ied] = sqrt(g*pj/uright.get(ied,0));
			//calculate normal velocities
			vni[ied] = (uleft.get(ied,1)*n[0] +uleft.get(ied,2)*n[1])/uleft.get(ied,0);
			vnj[ied] = (uright.get(ied,1)*n[0] + uright.get(ied,2)*n[1])/uright.get(ied,0);
		}

		// update residual and integ
		//std::cout << "Beginning new loop --- \n";
#pragma omp for
		for(a_int iel = 0; iel < m->gnelem(); iel++)
		{
			for(int ifael = 0; ifael < m->gnfael(iel); ifael++)
			{
				a_int ied = m->gelemface(iel,ifael);
				a_real len = m->ggallfa(ied,2);
				a_int nbdelem = m->gesuel(iel,ifael);

				if(nbdelem > iel) {
					for(int ivar = 0; ivar < NVARS; ivar++)
						residual(iel,ivar) -= fluxes(ied,ivar);
					integ(iel) += (fabs(vni[ied]) + ci[ied])*len;
				}
				else {
					for(int ivar = 0; ivar < NVARS; ivar++)
						residual(iel,ivar) += fluxes(ied,ivar);
					integ(iel) += (fabs(vnj[ied]) + cj[ied])*len;
				}
			}
		}
	} // end parallel region
}

a_real EulerFlow::compute_entropy()
{
	postprocess_cell();
	a_real vmaginf2 = uinf(0,1)/uinf(0,0)*uinf(0,1)/uinf(0,0) + uinf(0,2)/uinf(0,0)*uinf(0,2)/uinf(0,0);
	a_real sinf = ( uinf(0,0)*(g-1) * (uinf(0,3)/uinf(0,0) - 0.5*vmaginf2) ) / pow(uinf(0,0),g);

	amat::Array2d<a_real> s_err(m->gnelem(),1);
	a_real error = 0;
	for(int iel = 0; iel < m->gnelem(); iel++)
	{
		s_err(iel) = (scalars(iel,2)/pow(scalars(iel,0),g) - sinf)/sinf;
		error += s_err(iel)*s_err(iel)*m->gjacobians(iel)/2.0;
	}
	error = sqrt(error);

	//a_real h = sqrt((m->jacobians).max());
	a_real h = 1.0/sqrt(m->gnelem());

	std::cout << "EulerFV:   " << log10(h) << "  " << std::setprecision(10) << log10(error) << std::endl;

	return error;
}

void EulerFlow::postprocess_point()
{
	std::cout << "SpatialBase: postprocess_point(): Creating output arrays...\n";
	scalars.setup(m->gnpoin(),3);
	velocities.setup(m->gnpoin(),2);
	amat::Array2d<a_real> c(m->gnpoin(),1);

	amat::Array2d<a_real> areasum(m->gnpoin(),1);
	amat::Array2d<a_real> up(m->gnpoin(), NVARS);
	up.zeros();
	areasum.zeros();

	int inode, ivar;
	a_int ielem, iface, ip1, ip2, ipoin;

	for(ielem = 0; ielem < m->gnelem(); ielem++)
	{
		for(inode = 0; inode < m->gnnode(ielem); inode++)
			for(ivar = 0; ivar < NVARS; ivar++)
			{
				up(m->ginpoel(ielem,inode),ivar) += u.get(ielem,ivar)*m->garea(ielem);
				areasum(m->ginpoel(ielem,inode)) += m->garea(ielem);
			}
	}
	for(iface = 0; iface < m->gnbface(); iface++)
	{
		ielem = m->gintfac(iface,0);
		ip1 = m->gintfac(iface,2);
		ip2 = m->gintfac(iface,3);
		for(ivar = 0; ivar < NVARS; ivar++)
		{
			up(ip1,ivar) += ug.get(iface,ivar)*m->garea(ielem);
			up(ip2,ivar) += ug.get(iface,ivar)*m->garea(ielem);
			areasum(ip1) += m->garea(ielem);
			areasum(ip2) += m->garea(ielem);
		}
	}

	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		for(ivar = 0; ivar < NVARS; ivar++)
			up(ipoin,ivar) /= areasum(ipoin);

	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
	{
		scalars(ipoin,0) = up.get(ipoin,0);
		velocities(ipoin,0) = up.get(ipoin,1)/up.get(ipoin,0);
		velocities(ipoin,1) = up.get(ipoin,2)/up.get(ipoin,0);
		//velocities(ipoin,0) = dudx(ipoin,1);
		//velocities(ipoin,1) = dudy(ipoin,1);
		a_real vmag2 = pow(velocities(ipoin,0), 2) + pow(velocities(ipoin,1), 2);
		scalars(ipoin,2) = up.get(ipoin,0)*(g-1) * (up.get(ipoin,3)/up.get(ipoin,0) - 0.5*vmag2);		// pressure
		c(ipoin) = sqrt(g*scalars(ipoin,2)/up.get(ipoin,0));
		scalars(ipoin,1) = sqrt(vmag2)/c(ipoin);
	}
	std::cout << "EulerFV: postprocess_point(): Done.\n";
}

amat::Array2d<a_real> EulerFlow::getOutput() const
{
	return scalars;
}

}	// end namespace
