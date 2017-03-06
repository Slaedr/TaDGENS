#include "amatrix.hpp"

using namespace amat;

void printmatrix(const Matrix<double>& A)
{
    std::cout << A[1][2] << std::endl;
}

int matrix()
{
    Matrix<double> A(3,3);
    A.zeros();
    A(0,2) = 6.6;
    A[1][2] = 3;
    std::cout << A[0][2] << " " << A(0,2) << std::endl;
    printmatrix(A);
    
    return 0;
}

int main()
{
    matrix();
    
    return 0;
}
