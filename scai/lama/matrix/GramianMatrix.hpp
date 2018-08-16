/**
 * @file GramianMatrix.hpp
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
 * @brief Operator matrix class that stands for transpose( A ) * A but without building it explicitly
 * @author Thomas Brandes
 * @date 28.06.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/OperatorMatrix.hpp>
#include <scai/lama/_Vector.hpp>

namespace scai
{

namespace lama
{

/** Operator matrix class that stands for transpose( A ) * A without building it explicitly
 *
 *  The above matrix is not built explicitly and only some methods are implemented so
 *  this class can be used in solvers that exploit matrix-free methods.
 */
template<typename ValueType>
class GramianMatrix : public OperatorMatrix<ValueType>
{

public:

    /** Constructor that builds A' * A
     *
     *  @param[in] A is the rectangular matrix 
     *
     *  The size of the Gramian matrix is n x n if A has the size m x n
     */
    GramianMatrix( const Matrix<ValueType>& A ) :

        OperatorMatrix<ValueType>( A.getColDistributionPtr(), A.getColDistributionPtr() ),
        mA( A )

    {
        SCAI_LOG_INFO( logger, "GramianMatrix( A = " << A << " )" )
    }

    /** 
     *  Implementation of the linear operator
     */
    virtual void matrixTimesVectorDense(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>* y,
        const common::MatrixOp op ) const
    {
        SCAI_ASSERT_ERROR( !common::isConj( op ), "conj matrix operation not supported here." )

        // Note: transpose flag does not matter as marix is symmetric

        SCAI_LOG_INFO( logger, "matrixTimesVector, mA = " << mA )

        mAx    = mA * x;
        result = alpha * transpose( mA ) * mAx;  

        if ( y != nullptr )
        {
            result += beta * *y;
        }
    }

    /** Provide the context where linear operator is executed */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }

    /** Just use here the logger of the base class. */

    using OperatorMatrix<ValueType>::logger;

private:

    mutable DenseVector<ValueType> mAx;      // help vector for A * x
    const Matrix<ValueType>& mA;     // This class keeps only a reference
};

}

}
