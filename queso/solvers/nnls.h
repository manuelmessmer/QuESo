// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef NNLS_INCLUDE_H
#define NNLS_INCLUDE_H

//// External includes
#include <boost/numeric/ublas/matrix.hpp>

namespace queso {
namespace nnls {

typedef boost::numeric::ublas::matrix<double> MatrixType;
typedef boost::numeric::ublas::vector<double> VectorType;

// Wrapper for nnls solver
double nnls(MatrixType& A, const VectorType& b, VectorType& x);

} // End Namespace nnls
} // End namespace queso

#endif // NNLS_INCLUDE_H