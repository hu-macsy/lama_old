/**
 * @file GramianMatrix.cpp
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
 * @brief _Matrix class that stands for A' * A but without building it explicitly
 * @author Thomas Brandes
 * @date 28.06.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/AbstractMatrix.hpp>
#include <scai/lama/_Vector.hpp>

namespace scai
{

namespace lama
{

/** Abstract matrix class that stands for A' * A without building it explicitly
 *
 *  The above matrix is not built explicitly and only some methods are implemented so
 *  this class can be used in solvers that exploit matrix-free methods.
 */
class GramianMatrix : public AbstractMatrix
{

public:

    /** Constructor that builds A' * A
     *
     *  @param[in] A is the rectangular matrix 
     *
     *  The size of the Gramian matrix is n x n if A has the size m x n
     */
    GramianMatrix( const _Matrix& A ) :

        AbstractMatrix( A.getColDistributionPtr(), A.getColDistributionPtr() ),
        mA( A )

    {
        SCAI_LOG_INFO( logger, "GramianMatrix( A = " << A << " )" )

        mAx.reset( _Vector::getVector( VectorKind::DENSE, A.getValueType() ) );
    }

    /** Reimplement the matrix * vector operation
     *
     *  
     */
    virtual void matrixTimesVector(
        _Vector& result,
        const Scalar alpha,
        const _Vector& x,
        const Scalar beta,
        const _Vector& y ) const
    {
        SCAI_LOG_INFO( logger, "matrixTimesVector, mA = " << mA )

        *mAx= mA * x;
        result = *mAx * mA;  // same as mA' * mAx
        result *= alpha;
        result += beta * y;
    }

    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }

    /** This method must be provided so that solvers can decide about the type of additional runtime vectors. */

    virtual common::scalar::ScalarType getValueType() const
    {
        return mA.getValueType();
    }

private:

    _VectorPtr mAx;      // help vector for A * x
    const _Matrix& mA;
};

}

}
