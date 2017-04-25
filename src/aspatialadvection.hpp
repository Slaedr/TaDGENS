/** @file aspatialadvection.hpp
 * @brief Spatial discretization for linear advection
 * @author Aditya Kashi
 * @date 2017 April 14
 */

#ifndef __ASPATIALADVECTION_H
#define __ASPATIALADVECTION_H

#define NVARS 1

#include "aspatial.hpp"

namespace acfd {

/// Residual computation for linear advection
/** \note Make sure to call both compute_topological and compute_boundary_maps on the mesh object before using an object of this class!
 *
 * If the ODE is \f$ \frac{du}{dt} + R(u) = 0 \f$, [res](@ref res) holds \f$ R \f$.
 */
class LinearAdvection : public SpatialBase
{
protected:
	Vector a;								///< Advection velocity
	a_real amag;							///< Magnitude of advection velocity
	a_real bval;							///< Value to be imposed at inflow
	int inoutflow_flag;						///< Boundary flag at faces where inflow or outflow is required
	int extrapolation_flag;					///< Boundary flag for extrapolation condition
	amat::Array2d<a_real> output;			///< Pointwise values for output

	/// Computes upwind flux
	void computeNumericalFlux(const a_real* const uleft, const a_real* const uright, const a_real* const n, a_real* const flux);

	/// Computes face integrals from flow state described by the parameter
	void computeFaceTerms(const std::vector<Matrix>& ustage);

	/// Computes boundary (ghost) states depending on face marker for the face denoted by the first argument
	void computeBoundaryState(const int iface, const Matrix& instate, Matrix& bstate);

public:
	LinearAdvection(const UMesh2dh* mesh, const int _p_degree, const char basis, const Vector vel, const a_real b_val, const int inoutflag, const int extrapflag);
	
	/// Adds face contributions and computes domain contribution to the [right hand side](@ref residual) 
	void update_residual(const std::vector<Matrix>& ustage);
	
	/// Adds source term contribution to residual
	void add_source( a_real (*const rhs)(a_real, a_real, a_real), a_real t);

	/// Compute quantities to export
	void postprocess();

	/// Read-only access to output quantities
	const amat::Array2d<a_real>& getOutput() const {
		return output;
	}

	a_real computeL2Error(double (*const exact)(double,double,double), const double time) const;
};

}
#endif