/** \file
 * \brief Specification of PDE properties
 * \author Aditya Kashi
 */

#ifndef TADGENS_PDE_H
#define TADGENS_PDE_H

#include "aconstants.hpp"

namespace tadgens {

/// An abstract PDE context
/** Every PDE is assumed to have a convective term, a diffusive term and a source term.
 * Note that all functions operate not on one point (in the mesh) but a set of points.
 */
template <typename Derived>
class PDE
{
public:
	/// Computes the convective flux along a given direction
	virtual void convFlux(const ArrayBase<Derived>& u, const a_real *const n, a_real *const flux) const = 0;

	/// Computes the Jacobian of the convective flux along a direction w.r.t. the state variable(s)
	virtual void dConvFlux_du(const a_real *const u, const a_real *const n, a_real *const dfdu) const = 0;

	/// Computes the eigenvalues of the convective flux Jacobian along a given direction
	virtual void convFluxJacEigenvalues(const a_real *const u, const a_real *const n,
	                                    a_real *const eigs) const = 0;

	virtual void viscFlux(const a_real *const u, const a_real *const gradu, const a_real *const n,
	                      a_real *const flux) const = 0;

	virtual void dViscFlux_du(const a_real *const u, const a_real *const gradu, const a_real *const n,
	                          a_real *const dfdu) const = 0;

	virtual void dViscFlux_dgradu(const a_real *const u, const a_real *const gradu, const a_real *const n,
	                              a_real *const dfdgu) const = 0;

	virtual void source(const a_real *const u, a_real *const sourceterm) const = 0;

	virtual void dSource_du(const a_real *const u, a_real *const dsourcedu) const = 0;
};

}

#endif
