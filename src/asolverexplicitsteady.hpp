/** @file asolverexplicitsteady.hpp
 * @brief Solver class for explicit time-stepping for steady-state problems
 * @author Aditya Kashi
 * @date 2017-03-12
 */

#ifndef __ASOLVEREXPLICITSTEADY_H
#define __ASOLVEREXPLICITSTEADY_H

#ifndef __ASOLVERBASE_H
#include "asolverbase.hpp"
#endif

namespace acfd{

class SolverExplicitSteady : public SolverBase
{
public:
	void loaddata(acfd_real Minf, acfd_real vinf, acfd_real a, acfd_real rhoinf);
};

}

#endif
