/** @file aspatialadvection.hpp
 * @brief Spatial discretization for linear advection
 * @author Aditya Kashi
 * @date 2017 April 14
 */

#ifndef ASPATIALADVECTION_H
#define ASPATIALADVECTION_H

#include "aspatial.hpp"

namespace acfd {

/// Residual computation for linear advection
/** \note Make sure to call both compute_topological and compute_boundary_maps on the mesh object
 * before using an object of this class!
 *
 * If the ODE is \f$ \frac{du}{dt} + R(u) = 0 \f$, [res](@ref res) holds \f$ R \f$.
 */
class LinearAdvection : public SpatialBase
{
public:
	LinearAdvection(const UMesh2dh* mesh, const int _p_degree, const char basis,
	                const int inoutflag, const int extrapflag);

	/// Adds face contributions and computes domain contribution to the [right hand side](@ref residual)
	void update_residual(const std::vector<Matrix>& u, std::vector<Matrix>& res, std::vector<a_real>& mets);

	/// Compute quantities to export
	void postprocess(const std::vector<Matrix>& u);

	/// Read-only access to output quantities
	const amat::Array2d<a_real>& getOutput() const {
		return output;
	}

	/// provide the exact solution for a verification case
	a_real exact_solution(const a_real position[NDIM], const a_real time) const;

protected:
	std::array<a_real,NDIM> a;              ///< Advection velocity
	a_real amag;							///< Magnitude of advection velocity

	int inoutflow_flag;						///< Boundary flag at faces where inflow or outflow is required
	int extrapolation_flag;					///< Boundary flag for extrapolation condition
	amat::Array2d<a_real> output;			///< Pointwise values for output

	const int nvars;                        ///< Number of advected variables - currently 1

	const a_real aa,bb,dd,ee;               ///< Coeffs of exact solution

	/// Computes upwind flux
	void computeNumericalFlux(const a_real* const uleft, const a_real* const uright, const a_real* const n,
	                          a_real* const flux);

	/// Computes face integrals from flow state described by the parameter
	void computeFaceTerms(const std::vector<Matrix>& u);

	/// Computes boundary (ghost) states depending on face marker for the face denoted by the first argument
	void computeBoundaryState(const int iface, const Matrix& instate, Matrix& bstate);

	/// provide a test source term for a verification case
	a_real source_term(const a_real pos[NDIM], const a_real time) const
	{
		return a[0]*aa*std::cos(aa*pos[0]+dd)*std::sin(bb*pos[1]+ee)
			+ a[1]*bb*std::sin(aa*pos[0]+dd)*std::cos(bb*pos[1]*ee);
		//return 0;
	}

	/// State at inflow boundaries
	a_real bcfunc(const a_real position[NDIM]) const
	{
		return std::sin(aa*position[0]+dd)*std::sin(bb*position[1]+ee);
		//return sin(2*PI*position[1]);
	}
};

}
#endif
