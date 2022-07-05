// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef NNLS_INCLUDE_H
#define NNLS_INCLUDE_H

#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::matrix<double> MatrixType;
typedef boost::numeric::ublas::vector<double> VectorType;

namespace NNLS {

// Wrapper for nnls solver
double nnls(MatrixType& A, VectorType& b, VectorType& x);

} // End Namespace NNLS

#endif // NNLS_INCLUDE_H