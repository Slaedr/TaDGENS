/** @file areconstruction.cpp
 * @brief Implementations for reconstruction schemes
 * @author Aditya Kashi
 * @date April 23, 2017
 */

#include <areconstruction.hpp>

namespace acfd
{

Reconstructor::Reconstructor(int n_sur, const Element** sur_elems, const Matrix* sur_dofs) : nsur(n_sur), selms(sur_elems), sdofs(sur_dofs)
{ }

LeastSquaresReconstructor::LeastSquaresReconstructor(int n_sur, const Element** sur_elems, const Matrix* sur_dofs) : Reconstructor(n_sur, sur_elems, sur_dofs)
{
	const int reclevel = 1;
	// get degree before reconstruction
	int pdegree = selms[0]->getDegree() - 1;
	
	numknowndofs = 0; ntotaldofs = 0
	for(int i = 0; i <= pdegree; i++)
		numknowndofs += i+1;			// in 2D
	for(int i = 0; i <= pdegree+1; i++)	
		ntotaldofs += i+1;				// in 2D

	numders = 0;
	for(int i = 0; i < reclevel; i++)
		numders += pdegree + i+2;		// in 2D
	std::cout << " LeastSquaresReconstructor: Num known DOFs = " << numknowndofs << ", num total DOFs = " << ntotaldofs << ", num DOFs to reconstruct = " << numders << std::endl;
	if(numknowndofs != sdofs[0].cols()) std::cout << "! LeastSquaresReconstructor: No. of known DOFs and size of DOF input don't match!!\n";
	
	int nvars = sdofs[0].rows();
	
	if(pdegree == 0)
	A.resize(numknowndofs*nsur, numders);
	R.resize(numders,numders);

	// get center element's basis function values at center of each element
	// we need all basis functions, corresponding to both known and to-be-reconstructed DOFs
	cbasis.resize(nsur);
	const a_real* center = reinterpret_cast<TaylorElement*>(selms[0])->getCenter();
	const a_real* delta = reinterpret_cast<TaylorElement*>(selms[0])->getDelta();
	const std::vector<std::vector<a_real>>& bo = reinterpret_cast<TaylorElement*>(selms[0])->getBasisOffsets();
	for(int i = 0; i < nsur; i++) 
	{
		cbasis[i].resize(1, ntotaldofs);
		Matrix cc(1,NDIM);
		const a_real* centeri = reinterpret_cast<TaylorElement*>(selms[i+1])->getCenter();
		cc(0,0) = centeri[0]; cc(0,1) = centeri[1];
		getTaylorBasis(cc, pdegree+1, center, delta, bo, cbasis[i]);
	}

	// assemble
	for(int isur = 0; isur < nsur; isur++)
	{
		for(int i = 0; i < numders; i++)
			A(isur*numknowndofs,i) = cbasis[isur](0,i+numknowndofs);
		if(pdegree == 1) {
			A(isur*numknowndofs+1, 0) = cbasis[isur](0,1);
			A(isur*numknowndofs+1, 1) = cbasis[isur](0,2);
			A(isur*numknowndofs+1, 2) = 0;
			A(isur*numknowndofs+2, 0) = 0;
			A(isur*numknowndofs+2, 1) = cbasis[isur](0,1);
			A(isur*numknowndofs+2, 2) = cbasis[isur](0,2);
		}
	}

	R = A.transpose()*A;
	R = (R.inverse()).eval();
}

void LeastSquaresReconstructor::compute()
{
	int nvars = sdofs[0].rows();
	Matrix b = Matrix::Zero(numknowndofs*nsur, nvars);
	for(int isur = 0; isur < nsur; isur++)
	{
		for(int ivar = 0; ivar < nvars; ivar++) {
			// first row of the block
			for(int j = 0; j < pdegree+1; j++)
				b(isur*numknowndofs, ivar) -= sdofs[0].row(ivar) * cbasis[isur].block(0,0,1,numknowndofs).transpose();
			b(isur*numknowndofs, ivar) += sdofs[isur+1](ivar,0);

			// other rows
			if(pdegree == 1) {
				b(isur*numknowndofs+1, ivar) = selms[0]->getDelta()[0]/selms[isur+1]->getDelta()[0]*sdofs[isur+1](1,ivar) - sdofs[0](1,ivar);
				b(isur*numknowndofs+2. ivar) = selms[0]->getDelta()[1]/selms[isur+1]->getDelta()[1]*sdofs[isur+1](2,ivar) - sdofs[0](2,ivar);
			}
		}
	}

	b = (A.transpose()*b).eval();
	rder = R*b;
}

} // end namespace
