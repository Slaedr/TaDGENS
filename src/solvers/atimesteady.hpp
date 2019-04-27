/** @file atimesteady.hpp
 * @brief Steady-state pesudo-time stepping schemes
 * @author Aditya Kashi
 * @date 2017 April 18
 */

#ifndef ATIMESTEADY_H
#define ATIMESTEADY_H

#include "spatial/aspatial.hpp"

namespace acfd {

struct SteadySolverConfig {
	double tol;
	int maxiter;
	double cflbeg;
	double cflend;
};
	
/// TVD RK explicit time stepping
/** The initial condition must be specified in the spatial discretization context
 * [elsewhere](@ref SpatialBase::setInitialConditionModal).
 */
class SteadyBase
{
protected:
	const UMesh2dh *const m;						///< Mesh context
	SpatialBase *const spatial;                     ///< Spatial discretization context

	std::vector<Matrix> R;							///< Residuals

	/// vector of unknowns
	/** Each Eigen3 (E3) Matrix contains the DOF values of all physical variables for an element.
	 */
	std::vector<Matrix> u;

	/// Maximum allowable explicit time step for each element
	/** For Euler, stores (for each elem i)
	 *   Vol(i) / \f$ \sum_{j \in \partial\Omega_I} \int_j( |v_n| + c) d \Gamma \f$, 
	 * where v_n and c are average values of the cell faces
	 */
	std::vector<a_real> tsl;

	double cfl;										///< CFL number
	double tol;										///< Tolerance for residual
	int maxiter;									///< Max number of iterations

public:
	SteadyBase(const UMesh2dh*const mesh, SpatialBase* s,
	           a_real cflnumber, double toler, int max_iter);
	
	/// Read-only access to solution
	const std::vector<Matrix>& solution() const {
		return u;
	}

	/// Carries out the time stepping process
	virtual void solve() = 0;
};

/// Explicit forward-Euler scheme with local time stepping
class SteadyExplicit : public SteadyBase
{
public:
	/** \param[in] mesh The mesh context
	 * \param[in] s The spatial discretization context
	 * \param[in] cflnumber
	 * \param[in] toler Tolerance for the relative residual
	 * \param[in] max_iter Maximum number of iterations
	 */
	SteadyExplicit(const UMesh2dh *const mesh, SpatialBase *const s,
	               a_real cflnumber, double toler, int max_iter);

	/// Carries out the time stepping process
	void solve();
};

/// Implicit backward-Euler pseudo-time scheme with local time stepping
class SteadyImplicit : public SteadyBase
{
public:
	/** \param[in] mesh The mesh context
	 * \param[in] s The spatial discretization context
	 * \param[in] cflnumber
	 * \param[in] toler Tolerance for the relative residual
	 * \param[in] max_iter Maximum number of iterations
	 */
	SteadyImplicit(const UMesh2dh *const mesh, SpatialBase *const s,
	               const a_real cflnumber, double toler, int max_iter);

	/// Carries out the time stepping process
	void solve();
};

}
#endif
