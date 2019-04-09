/** @file aspatialpoisson.hpp
 * @brief Poisson equation solver
 * @author Aditya Kashi
 * @date 2017-04-10
 */

#ifndef ASPATIALPOISSON_H
#define ASPATIALPOISSON_H 1

#include "aspatial.hpp"

#include <Eigen/Sparse>

namespace acfd {

/// Symmetric interior penalty scheme for Laplace operator
/** \note Strong boundary conditions, hence only nodal basis!
 * Currently only Dirichlet boundaries.
 */
class LaplaceSIP : public SpatialBase<1>
{
protected:
	a_real nu;											///< Diffusivity
	a_real eta;											///< Penalty
	a_real (*const rhs)(a_real, a_real);				///< forcing function
	a_real (*const exact)(a_real, a_real);              ///< Exact solution
	a_real (*const exactgradx)(a_real, a_real);			///< Exact x-derivative of exact solution
	a_real (*const exactgrady)(a_real, a_real);			///< Exact y-derivative of exact solution
	int dirichlet_id;									///< Boundary marker for Dirichlet boundary
	int neumann_id;										///< Boundary marker for homogeneous Neumann boundary
	a_real dirichlet_value;								///< Dirichlet boundary value
	std::vector<a_int> dirdofflags;						///< Binary flag for each DOF, identifying as lying on a Dirichlet boundary or not
	a_int ndirdofs;										///< Number of Dirichlet DOFs
	a_real cbig;										///< Penalty for Dirichlet condition

	Eigen::SparseMatrix<a_real> Ag;						///< Global left hand side matrix
	Vector bg;											///< Global load vector
	Vector ug;											///< 'Global' solution vector
	amat::Array2d<a_real> output;						///< Output array for plotting

public:
	LaplaceSIP(const UMesh2dh* mesh, const int _p_degree, a_real stab_param,
			a_real(*const f)(a_real,a_real), a_real(*const exact_sol)(a_real,a_real), 
			a_real(*const exact_gradx)(a_real,a_real), a_real(*const exact_grady)(a_real,a_real));

	void assemble();
	void solve();

	/// Computes errors in L2 and SIP norms
	void computeErrors(a_real& l2error, a_real& siperror) const;

	void postprocess(const std::vector<Matrix>& u);

	const amat::Array2d<a_real>& getOutput() const {
		return output;
	}

	void update_residual(const std::vector<Matrix>& u,
	                     std::vector<Matrix>& res, 
	                     std::vector<a_real>& mets)
	{ }
};

}	// end namespace
#endif
