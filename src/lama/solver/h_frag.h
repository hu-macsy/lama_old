/*
 * h_frag.h
 *
 *  Created on: 20.07.2011
 *      Author: rrehrman
 */

#ifndef LAMA_H_FRAG_H_
#define LAMA_H_FRAG_H_

#include <lama/matrix/DenseMatrix.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>

template<typename T>
inline bool matrixIs( const lama::Matrix& matrix )
{
    return typeid(lama::DenseMatrix<T>) == typeid(matrix) || typeid(lama::CSRSparseMatrix<T>) == typeid(matrix)
                    || typeid(lama::ELLSparseMatrix<T>) == typeid(matrix)
                    || typeid(lama::COOSparseMatrix<T>) == typeid(matrix)
                    || typeid(lama::DIASparseMatrix<T>) == typeid(matrix)
                    || typeid(lama::JDSSparseMatrix<T>) == typeid(matrix);
}

inline bool matrixIsDense( const lama::Matrix& matrix )
{
    return typeid(lama::DenseMatrix<float>) == typeid(matrix) || typeid(lama::DenseMatrix<double>) == typeid(matrix);
}

#endif /* LAMA_H_FRAG_H_ */
