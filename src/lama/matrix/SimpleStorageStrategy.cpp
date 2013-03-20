/**
 * @file SimpleStorageStrategy.cpp
 *
 * @license
 * Copyright (c) 2012
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.using phoenix::val;
 * @endlicense
 *
 * @brief SimpleStorageStrategy.cpp
 * @author Kai Buschulte
 * @date 07.05.2012
 * $Id$
 */

// hpp
#include <lama/matrix/SimpleStorageStrategy.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/ContextManager.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/CommunicatorManager.hpp>

// boost
#include <boost/config/warning_disable.hpp>
// spirit
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/object/new.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/statement/if.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace lama
{
LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, SimpleStorageStrategy<ValueType>::logger,
                              "Matrix.SimpleStorageStrategy" );

template<typename ValueType>
SimpleStorageStrategy<ValueType>::SimpleStorageStrategy( Matrix& other ) :
    Matrix( other )
{
}

template<typename ValueType>
SimpleStorageStrategy<ValueType>::~SimpleStorageStrategy()
{
}

/***** Strategy decision method *****/

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::applyStrategy( const Matrix )
{
    if ( getContextPtr()->getType() == Context::CUDA )
    {
        mInnerMatrix.reset( new ELLSparseMatrix( getDistributionPtr(), getColDistributionPtr() ) );
    }
    else
    {
        mInnerMatrix.reset( new CSRSparseMatrix( getDistributionPtr(), getColDistributionPtr() ) );
    }

    LAMA_LOG_INFO( logger, "assign " << matrix << " to " << *this );

    const _SparseMatrix* sparseMatrix = dynamic_cast<const _SparseMatrix*>( &matrix );

    if ( sparseMatrix )
    {
        // for a sparse matrix local + halo part can be assigned

        Matrix::setDistributedMatrix( matrix.getDistributionPtr(), matrix.getColDistributionPtr() );

        // ToDo: allow flexibility regarding the context, e.g. format conversion should be done on GPU

        mLocalData->assign( matrix.getLocalStorage() );
        mHaloData->assign( matrix.getHaloStorage() );
        mHalo = matrix.getHalo();
    }
    else
    {
        // convert dense matrix to a sparse matrix

        DistributionPtr colDist = matrix.getColDistributionPtr();
        DistributionPtr rowDist = matrix.getDistributionPtr();

        Matrix::setDistributedMatrix( rowDist, colDist );

        matrix.buildLocalStorage( *mLocalData ); // local storage with all columns

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }
}

/***** Separate configuration methods *****/

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::setContext( const ContextPtr context )
{
    setContext( context );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    setDistributedMatrix( rowDistribution, colDistribution );
}

/***** Forward Methods SimpleStorageStrategy *****/

template<typename ValueType>
const char* SimpleStorageStrategy<ValueType>::getTypeName() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getTypeName();
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::clear()
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->clear();
}
template<typename ValueType>
void SimpleStorageStrategy<ValueType>::allocate( const IndexType numRows, const IndexType numColumns )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->allocate( numRows, numColumns );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->allocate( rowDistribution, colDistribution );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::setIdentity()
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->setIdentity();
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::assign( const Matrix& other )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->assign( other );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::assign( const _MatrixStorage& other )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->assign( other );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::assign(
    const _MatrixStorage& storage,
    DistributionPtr rowDist,
    DistributionPtr colDist )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->assign( storage, rowDist, colDist );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::buildLocalStorage( _MatrixStorage& storage ) const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->buildLocalStorage( storage );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::getRow( Vector& row, const IndexType globalRowIndex ) const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->getRow( row, globalRowIndex );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::getDiagonal( Vector& diagonal ) const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->getDiagonal( diagonal );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::setDiagonal( const Vector& diagonal )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->setDiagonal( diagonal );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::setDiagonal( const Scalar scalar )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->setDiagonal( scalar );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::scale( const Vector& scaling )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->scale( scaling );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::scale( const Scalar scaling )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->scale( scaling );
}

template<typename ValueType>
Scalar SimpleStorageStrategy<ValueType>::getValue( IndexType i, IndexType j ) const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getValue( i, j );
}

template<typename ValueType>
IndexType SimpleStorageStrategy<ValueType>::getNumValues() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getNumValues();
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::matrixTimesVector(
    Vector& result,
    const Scalar alpha,
    const Vector& x,
    const Scalar beta,
    const Vector& y ) const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->matrixTimesVector( result, alpha, x, beta, y );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::matrixTimesScalar( const Matrix& other, const Scalar alpha )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->matrixTimesScalar( other, alpha );
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::matrixTimesMatrix(
    Matrix& result,
    const Scalar alpha,
    const Matrix& B,
    const Scalar beta,
    const Matrix& C ) const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->matrixTimesMatrix( result, alpha, B, beta, C );
}

template<typename ValueType>
IndexType SimpleStorageStrategy<ValueType>::getLocalNumValues() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getLocalNumValues();
}

template<typename ValueType>
IndexType SimpleStorageStrategy<ValueType>::getLocalNumRows() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getLocalNumRows();
}

template<typename ValueType>
IndexType SimpleStorageStrategy<ValueType>::getLocalNumColumns() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getLocalNumColumns();
}

template<typename ValueType>
ContextPtr SimpleStorageStrategy<ValueType>::getContextPtr() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getContextPtr();
}

template<typename ValueType>
Matrix::MatrixKind SimpleStorageStrategy<ValueType>::getMatrixKind() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getMatrixKind();
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::prefetch() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->prefetch();
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::wait() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->wait();
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::invert( const Matrix& other )
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->invert( other );
}

template<typename ValueType>
Scalar SimpleStorageStrategy<ValueType>::maxDiffNorm( const Matrix& other ) const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->maxDiffNorm( other );
}

template<typename ValueType>
std::auto_ptr<Matrix> SimpleStorageStrategy<ValueType>::create() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->create();
}

template<typename ValueType>
std::auto_ptr<Matrix> SimpleStorageStrategy<ValueType>::copy() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->copy();
}

template<typename ValueType>
Scalar::ScalarType SimpleStorageStrategy<ValueType>::getValueType() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getValueType();
}

template<typename ValueType>
bool SimpleStorageStrategy<ValueType>::hasDiagonalProperty() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->hasDiagonalProperty();
}

template<typename ValueType>
void SimpleStorageStrategy<ValueType>::resetDiagonalProperty()
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    mInnerMatrix->resetDiagonalProperty();
}

template<typename ValueType>
size_t SimpleStorageStrategy<ValueType>::getMemoryUsage() const
{
    LAMA_ASSERT( mInnerMatrix, "Matrix is not created. Already assigned a matrix to the SimpleStorageStrategy?" );
    return mInnerMatrix->getMemoryUsage();
}

template class LAMA_DLL_IMPORTEXPORT SimpleStorageStrategy<float> ;
template class LAMA_DLL_IMPORTEXPORT SimpleStorageStrategy<double> ;

} //namespace lama
