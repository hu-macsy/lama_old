/**
 * @file Matrix.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Abstract base class for all matrices of a certain value type.
 * @author Thomas Brandes
 * @date 31.10.2017
 */
#pragma once

#include <scai/lama/matrix/_Matrix.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Abstract base class for a matrix of a certain value type.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Matrix:

    public _Matrix

{

public:

    /** Create a new matrix of a certain format but with same value type */

    static Matrix<ValueType>* getMatrix( const Format format );

    virtual ~Matrix();

    /** Overwrite _Matrix::newMatrix to get the covariant return type */

    virtual Matrix<ValueType>* newMatrix( void ) const = 0;

    /** Overwrite _Matrix::copy to get the covariant return type */

    virtual Matrix<ValueType>* copy( void ) const = 0;

    // Important: pass assign operator so it can also be used here

    using _Matrix::operator=;

    /** 
     *  Implementation of pure method of _Matrix::getValueType 
     *
     *  The following code demonstrates how this routine might be used to
     *  make a safe reinterpret_cast from the base class _Matrix.
     *
     *  \code
     *    _Matrix& m = ...   
     *    if ( m.getValueType() == TypeTraits<T>::stype )
     *    {
     *        Matrix<T>& mt = reinterpret_cast<Matrix<T>&>( m );
     *    }
     *  \endcode
     */
    virtual common::scalar::ScalarType getValueType() const;

    /** Implementation of pure method _Matrix::getValueTypeSize */

    virtual size_t getValueTypeSize() const;

    /** Implementation of _Matrix::matrixTimesVector */

    void matrixTimesVector(
        _Vector& result,
        const Scalar alpha,
        const _Vector& x,
        const Scalar beta,
        const _Vector& y ) const;

    /** Implementation of _Matrix::vectorTimesMatrix */

    void vectorTimesMatrix(
        _Vector& result,
        const Scalar alpha,
        const _Vector& x,
        const Scalar beta,
        const _Vector& y ) const;

    /** Implementation of _Matrix::setRow for all typed matrices
     *
     *  Note: all derived classes must provide setLocalRow( rowArray, localRowIndex, op )
     */
    void setRow( const _Vector& row, const IndexType globalRowIndex,
                 const common::BinaryOp op );

    /** Implementation of _Matrix::setColumn for all typed matrices
     *
     *  The method is implemented by setting the local part of the column on each partition.
     *  All derived classes must provide setLocalColum( colArray, colIndex, op )
     */
    void setColumn( const _Vector& column,
                    const IndexType colIndex,
                    const common::BinaryOp op );

    /**
     * @brief Returns the L1 norm of this matrix.
     *
     * @return the L1 norm of this as real value 
     *
     * l1Norm computes the sum of the absolute values of all entries
     */
    virtual NormType<ValueType> l1Norm( void ) const = 0;
    virtual NormType<ValueType> l2Norm( void ) const = 0;
    virtual NormType<ValueType> maxNorm( void ) const = 0;

    /** 
     * we provide here an implementation that works for any kind of
     * matrices.
     */
    virtual NormType<ValueType> maxDiffNorm( const _Matrix& other ) const = 0;

protected:

    Matrix();

    /**
     * @brief Constructs a matrix of a given size replicated on each partition.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     *
     * Same as _Matrix( NoDistribution(numRows), NoDistribution(numColumns) )
     */
    Matrix( const IndexType numRows, const IndexType numColumns );

    Matrix( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    Matrix( const _Matrix& other, dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    Matrix( const _Matrix& other );

    Matrix( const Matrix<ValueType>& other );

    /** typed version of matrixTimesVector must be implemented by derived classes */

    virtual void matrixTimesVectorImpl(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>& denseY ) const = 0;

    virtual void vectorTimesMatrixImpl(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>& denseY ) const = 0;

    virtual void setLocalRow( const hmemo::HArray<ValueType>& row,
                              const IndexType localRowIndex,
                              const common::BinaryOp op  ) = 0;

    virtual void setLocalColumn( const hmemo::HArray<ValueType>& column,
                                 const IndexType colIndex,
                                 const common::BinaryOp op  ) = 0;

    /** This method is the same for dense/sparse matrices as column distribution is replicated */

    void vectorTimesMatrixRepCols(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>& denseY ) const;

    // Implementations of pure _Matrix methods to guarantee upward compatibilty

    Scalar _l1Norm()  const;
    Scalar _l2Norm()  const;
    Scalar _maxNorm() const;
    Scalar _maxDiffNorm( const _Matrix& other ) const;
};

/** 
 * Definiton of corresponding shared pointer type for the class Matrix<ValueType> by a type alias.
 *
 *  \code
 *      MatrixPtr<ValueType> x( Matrix<ValueType>::getMatrix( Format::COO ) );
 *      std::shared_ptr<Matrix<ValueType> > x( Matrix<ValueType>::getMatrix( Format::COO ) );
 *  \endcode
*/
template<typename ValueType>
using MatrixPtr = std::shared_ptr<Matrix<ValueType> >;

} /* end namespace lama */

} /* end namespace scai */

