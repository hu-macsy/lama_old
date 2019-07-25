/**
 * @file HybridMatrix.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Matrix class that stands for A1 + A2 but without building it explicitly
 * @author Thomas Brandes
 * @date 28.06.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/OperatorMatrix.hpp>

namespace scai
{

namespace lama
{

/** Operator matrix class that stands for A1 + A2 without building it explicitly
 *
 *  The above matrix is not built explicitly and only some methods are implemented so
 *  this class can be used in solvers that exploit matrix-free methods.
 *
 *  This class might be very useful to combine a stencil matrix with a sparse matrix that
 *  contains only the boundary conditions.
 */
template<typename ValueType>
class HybridMatrix : public OperatorMatrix<ValueType>
{

public:

    /** Constructor that builds the matrix A1 + A2 symbollically
     *
     *  @param[in] A1, A2 are the two matrices that build the hybrid one.
     *
     *  The distributions of A1 and A2 must be the same and are the corresponding
     *  distributions for this matrix.
     *
     *  Note: some updates operations like setContextPtr, scale will update the original
     *        matrices A1, A2 that are only kept as reference.
     */
    HybridMatrix( Matrix<ValueType>& A1, Matrix<ValueType>& A2 ) :

        OperatorMatrix<ValueType>( A1.getRowDistributionPtr(), A1.getColDistributionPtr() ),
        mA1( A1 ),
        mA2( A2 )

    {
        SCAI_ASSERT_EQ_ERROR( A1.getRowDistribution(), A2.getRowDistribution(), "HybridMatrix: size mismatch" )
        SCAI_ASSERT_EQ_ERROR( A1.getColDistribution(), A2.getColDistribution(), "HybridMatrix: size mismatch" )

        tmpTarget.reset( A1.newTargetVector() );
        tmpSource.reset( A1.newSourceVector() );
    }

    /** 
     *  Implementation of the linear operator
     */
    virtual void matrixTimesVector(
        Vector<ValueType>& result,
        const ValueType alpha,
        const Vector<ValueType>& x,
        const ValueType beta,
        const Vector<ValueType>* y,
        const common::MatrixOp op ) const
    {
        // result = alpha * ( A1 + A2 ) * x +  beta * y 
        //  ->  result = alpha * A1 * x + [ alpha * A2 * x + beta * y ]

        Vector<ValueType>& tmp = common::isTranspose( op ) ? *tmpSource : *tmpTarget;

        SCAI_LOG_DEBUG( logger, "gemv 1, op = " << op << ", A2 = " << mA2 )
        mA2.matrixTimesVector( tmp, alpha, x, beta, y, op );
        SCAI_LOG_DEBUG( logger, "gemv 2, op = " << op << ", A1 = " << mA1 )
        mA1.matrixTimesVector( result, alpha, x, ValueType( 1 ), &tmp, op );
    }

    /**
     *  @brief Implementation of pure method Matrix<ValueType>::scale 
     *
     *  Be careful, this operation has side effects for the both matrices used here.
     */
    virtual void scale( const ValueType& alpha )
    {
        mA1.scale( alpha );
        mA2.scale( alpha );
    }

    /**
     * @brief Provide implementation for _Matrix::setContextPtr 
     *
     * The better approach might be to set the context for the 
     * original two matrices.
     */
    virtual void setContextPtr( const hmemo::ContextPtr context ) 
    {
        mA1.setContextPtr( context );
        mA2.setContextPtr( context );
    }

    virtual void matrixTimesVectorDense(
        DenseVector<ValueType>&,
        const ValueType,
        const DenseVector<ValueType>&,
        const ValueType,
        const DenseVector<ValueType>*,
        const common::MatrixOp ) const 
    {
        COMMON_THROWEXCEPTION( "NEVER CALLED" )
    }

    /** Provide the context where linear operator is executed 
     *
     *  This seems a bit ambiguous, we could take one of the matrices
     */
    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA1.getContextPtr();
    }

    virtual size_t getMemoryUsage() const
    {
        return mA1.getMemoryUsage() + mA2.getMemoryUsage();
    }

    /* Implementation of method writeAt for hybrid matrix. */

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "HybridMatrix[ " << mA1 << " + " << mA2 << " ]";
    }

    /** Just use here the logger of the base class. */

    using OperatorMatrix<ValueType>::logger;

private:

    std::unique_ptr<Vector<ValueType>> tmpTarget;  // temporary vector, will be reused
    std::unique_ptr<Vector<ValueType>> tmpSource;  // temporary vector, will be reused

    Matrix<ValueType>& mA1;     // This class keeps only a reference
    Matrix<ValueType>& mA2;     // This class keeps only a reference
};

}

}
