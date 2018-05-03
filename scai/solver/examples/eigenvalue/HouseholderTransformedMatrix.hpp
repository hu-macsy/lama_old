/**
 * @file HouseholderTransformedMatrix.cpp
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
 * @brief _Matrix class that builds HLH for a Householder matrix H = I - u * u' / alpha
 * @author Thomas Brandes
 * @date 28.06.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/OperatorMatrix.hpp>

namespace scai
{

namespace lama
{

/** Matrix class that stands for H * L * H with H = I - u * u' / alpha 
 *
 *  The above matrix is not built explicitly and only some methods are implemented so
 *  this class can be used in solvers that exploit matrix-free methods.
 */
template<typename ValueType>
class HouseholderTransformedMatrix : public OperatorMatrix<ValueType>
{

public:

    /** Constructor that builds HLH with H = 1 - u * u' * alpha 
     *
     *  @param[in] L is the matrix that to which the Householder matrix is applied
     *  @param[in] u is the vector of the Householder matrix
     *  @param[in] alpha is the scaling factor.
     */

    HouseholderTransformedMatrix( const Matrix<ValueType>& L, const Vector<ValueType>& u, const ValueType alpha ) : 

        OperatorMatrix<ValueType>( L ),
        mL( L )

    {
        SCAI_LOG_INFO( logger, "HouseholderTransformedMatrix( L = " << L << ", u = " << u << ", alpha = " << alpha << " )" )

        // build the help vectors mR and mS that are used within the matrix * vector operation

        DenseVector<ValueType> h( L.getContextPtr() );

        h = L * u;
        h /= alpha;

        ValueType gamma = u.dotProduct( h ) / alpha * 0.5;
        
        mR = u;
        mS = h - gamma * u;

        mS[0] = ValueType( 0 );
        mR[0] = ValueType( 0 );
    }

    /** Reimplement the matrix * vector operation
     *
     *  H L H * x = L * x - s' * x * r - r' * x * s
     */
    virtual void matrixTimesVectorDense(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>* y,
        common::MatrixOp op ) const
    {
        if ( op != common::MatrixOp::NORMAL )
        {
            COMMON_THROWEXCEPTION( op << " for matrixTimesVector not supported here" )
        }

        SCAI_LOG_INFO( logger, "matrixTimesVector, mL = " << mL )

        result = mL * x;
        result -= mS.dotProduct( x ) * mR;
        result -= mR.dotProduct( x ) * mS;
        result[0] = ValueType( 0 );
        result *= alpha;
 
        if ( y != nullptr )
        {
            result += beta * *y;
        }
    }

    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mL.getContextPtr();
    }

    using OperatorMatrix<ValueType>::logger;

private:

    DenseVector<ValueType> mR;   // help vector to make matrix * vector more efficient
    DenseVector<ValueType> mS;   // help vector to make matrix * vector more efficient

    const Matrix<ValueType>& mL;
};

}

}
