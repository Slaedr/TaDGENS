/** @file aspatialpoisson_continuous.hpp
 * @brief Poisson equation solver
 * @author Aditya Kashi
 * @date 2017-04-10
 */

#ifndef ASPATIALPOISSON_CONTINUOUS_H
#define ASPATIALPOISSON_CONTINUOUS_H 1

#include "spatial/aspatial.hpp"

#include <Eigen/Sparse>

namespace tadgens {

/// Continuous FEM for Poisson equation
class LaplaceC : public SpatialBase
{
public:
	LaplaceC(const UMesh2dh *const mesh, const int _p_degree,
	         const int bc_dirichlet_id, const int bc_neumann_id);

	void assemble();
	void solve();

	/// Computes errors in L2 and H1 norms
	void computeErrors(a_real& l2error, a_real& h1error) const;

	void postprocess(const std::vector<Matrix>& u);

	const amat::Array2d<a_real>& getOutput() const {
		return output;
	}

	void update_residual(const std::vector<Matrix>& u, std::vector<Matrix>& res, std::vector<a_real>& mets) { }

	/// provide the exact solution for a verification case
	/** \warning This is an eigenfunction of the Laplace operator!
	 */
	a_real exact_solution(const a_real pos[NDIM], const a_real time) const
	{
		return std::sin(aa*pos[0]+dd)*std::sin(bb*pos[1]+ee);
	}

protected:
	const a_real nu;                        ///< Diffusivity (currently always 1)
	const a_real cbig;                      ///< Penalty for Dirichlet condition
	const int dirichlet_id;                 ///< Boundary marker for Dirichlet boundary
	const int neumann_id;                   ///< Boundary marker for homogeneous Neumann boundary
	const a_int ntotaldofs;                 ///< Total number of DOFs
	const a_real aa, bb, dd, ee;            ///< Params of exact solution

	//amat::Array2d<int> dofmap;              ///< Identifies global dof index with local dofs and element index
	std::vector<int> bflag;                   ///< Marks whether a DOF lies on a boundary or not

	Eigen::SparseMatrix<a_real> Ag;         ///< Global left hand side matrix
	Vector bg;                              ///< Global load vector
	Vector ug;                              ///< 'Global' solution vector
	amat::Array2d<a_real> output;           ///< Output array for plotting


	a_int getGlobalDofIdx(const a_int elem, const int localnode) const;

	a_real source_term(const a_real pos[NDIM], const a_real time) const
	{
		return (aa*aa+bb*bb)*exact_solution(pos,time);
	}

	std::array<a_real,NDIM> exact_gradient(const a_real pos[NDIM], const a_real time) const
	{
		std::array<a_real,2> grad;
		grad[0] = aa*std::cos(aa*pos[0]+dd)*std::sin(bb*pos[1]+ee);
		grad[1] = std::sin(aa*pos[0]+dd)*bb*std::cos(bb*pos[1]+ee);
		return grad;
	}
};

}

#endif
