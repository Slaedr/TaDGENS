/** @file aspatialadvection.cpp
 * @brief Implementatio of spatial discretization for linear advection
 * @author Aditya Kashi
 * @date 2017 April 14
 */

#include "aspatialadvection.hpp"

namespace acfd {

LinearAdvection::LinearAdvection(const UMesh2dh* mesh, const int _p_degree, const char basis, const Vector vel, const a_real b_val, const int inoutflag, const int extrapflag)
	: SpatialBase(mesh, _p_degree, basis), a(vel), bval(b_val), inoutflow_flag(inoutflag), extrapolation_flag(extrapflag)
{
	computeFEData();
	u.resize(m->gnelem());
	res.resize(m->gnelem());
	for(int iel = 0; iel < m->gnelem(); iel++) {
		u[iel].resize(1, elems[iel].getNumDOFs());
		res[iel].resize(1, elems[iel].getNumDOFs());
	}
}

void computeBoundaryState(const int iface, const Vector& instate, Vector& bstate)
{
	if(m->gintfacbtags(iface, 0) == extrapolation_flag)
	{
		bstate = instate;
	}
	else if(m->gintfacbtags(iface, 0) == extrapolation_flag)
	{
		// TODO: compute normal velocity and decide whether to extrapolate or impose specified boundary value
	}
}

void LinearAdvection::computeFaceTerms()
{
	for(int iel = 0; iel < m->gnbface(); iel++)
	{
		;
	}
}

}
