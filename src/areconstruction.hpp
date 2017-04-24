/** @file areconstruction.hpp
 * @brief Classes for higher-order derivative reconstruction schemes
 * @author Aditya Kashi
 * @date April 23, 2017
 */

#ifndef __ARECONSTRUCTION_H
#define __ARECONSTRUCTION_H 1

#ifndef __AELEMENTS_H
#include "aelements.hpp"
#endif

namespace acfd
{

/// Abstract class for a erivative reconstruction context for a single element
class Reconstructor
{
protected:
	/// Number of face-neighboring elements
	int nsur;

	/// Surrounding elements
	const Element** selms;

	/// DOFs of surrounding elements
	/** The first matrix in the list should be the DOFs of this element, ie,
	 * the element whose derivatives this object is to reconstruct.
	 */
	const Matrix* sdofs;

	/// Reconstructed derivatives
	Matrix rder;

public:
	Reconstructor(const Element** sur_elems, const Matrix* sur_dofs);
	virtual ~Reconstruction();

	/// Computes reconstructed derivatives corresponding to the set [DOFs](@ref sdofs)
	virtual void compute() = 0;

	const Matrix& recDOFs() const {
		return rder;
	}
};

/// Least-squares reconstruction, currently only for RDG-P0P1 and P1P2.
/** Only works for DOFs corresponding to Taylor basis.
 * The Element pointers must point to TaylorElement objects or Bad Things (TM) will happen.
 *
 * The least-squares problem is solved by normal equations. The left-hand side is pre-inverted and stored.
 */
class LeastSquaresReconstructor : public Reconstructor
{
	int ntotaldofs;
	int numknowndofs;
	int numders;				///< numders = ntotaldofs - numknowndofs
	int pdegree;				///< Polynomial degree of underlying DG discretization
	Matrix R;
	Matrix A;
	/// center element's basis function values at center of each surrounding element
	std::vector<Matrix> cbasis;
public:
	/// Computes least-squares reconstruction operator and inverts it
	LeastSquaresReconstructor(const Element** sur_elems, const Matrix* sur_dofs);

	/// Assemly of RHS and matrix-vector multiplication to compute reconstructed DOFs
	void compute();
};


} // end namespace
#endif
