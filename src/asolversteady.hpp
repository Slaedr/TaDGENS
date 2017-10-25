/** @file asolversteady.hpp
 * @brief Solver class for time-stepping for steady-state problems
 * @author Aditya Kashi
 * @date 2017-03-12
 */

#ifndef ASOLVERSTEADY_H
#define ASOLVERSTEADY_H

#include "aspatial.hpp"

namespace acfd{

/// Abstract interface for solvers that iterate to steady-state
class SolverSteady
{
protected:
	SolverBase* space;				///< Spatial discretization context
	acfd_real icfl;
	acfd_real fcfl;
	acfd_real tol;
	int maxits;
	int rampstart;
	int rampend;

public:
	void loadFlowData(acfd_real Minf, acfd_real vinf, acfd_real a, acfd_real rhoinf);
	void loadTimeParameters(acfd_real toler, int maxiter, acfd_real icfl, acfd_real fcfl, int rampStartIter, int rampEndIter);
	virtual void solve() = 0;
};

/// Explicit pseudo-time forward Euler time-stepping with local time steps
class SolverSteadyExplicit : public SolverSteady
{
public:
	void solve();
};

}

#endif
