/**
 * @file Kaczmarz.cpp
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
 * @brief Implementation of methods for the Kaczmarz solver.
 * @author Thomas Brandes
 * @date 24.08.2011
 */

// hpp
#include <scai/solver/Kaczmarz.hpp>

// internal scai libraries
#include <scai/lama/Vector.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, Kaczmarz<ValueType>::logger, "Solver.IterativeSolver.Kaczmarz" )

using lama::Matrix;
using lama::MatrixKind;
using lama::Vector;
using lama::VectorKind;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* Kaczmarz<ValueType>::create()
{
    return new Kaczmarz<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType Kaczmarz<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "Kaczmarz" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
Kaczmarz<ValueType>::Kaczmarz( const std::string& id ) : 

    IterativeSolver<ValueType>( id )
{
}

template<typename ValueType>
Kaczmarz<ValueType>::Kaczmarz( const std::string& id, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id, logger )
{
}

template<typename ValueType>
Kaczmarz<ValueType>::Kaczmarz( const Kaczmarz& other ) : 

    IterativeSolver<ValueType>( other )
{
}

template<typename ValueType>
Kaczmarz<ValueType>::~Kaczmarz()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void Kaczmarz<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.Kaczmarz.initialize" )

    IterativeSolver<ValueType>::initialize( coefficients );

    KaczmarzRuntime& runtime = getRuntime();

    // matrix type dense or sparse decides whether vector mRow is dense or sparse

    VectorKind kind = VectorKind::DENSE;

    if ( coefficients.getMatrixKind() == MatrixKind::SPARSE )
    {
        kind = VectorKind::SPARSE;
    }

    runtime.mRow.reset( Vector<ValueType>::getVector( kind ) );
    runtime.mRow->allocate( coefficients.getRowDistributionPtr() );
}

/* ========================================================================= */
/*    IterativeSolver: one iteration                                         */
/* ========================================================================= */

template<typename ValueType>
void Kaczmarz<ValueType>::iterate()
{
    SCAI_REGION( "Solver.Kaczmarz.iterate" )

    KaczmarzRuntime& runtime = getRuntime();

    IndexType iter = this->getIterationCount();

    if ( iter == 0 )
    {
        this->getResidual();
    }

    const Matrix<ValueType>& A = *runtime.mCoefficients;

    const Vector<ValueType>& b = *runtime.mRhs;

    Vector<ValueType>& x = runtime.mSolution.getReference(); // -> dirty

    Vector<ValueType>& z = *runtime.mRow;

    SCAI_LOG_INFO( logger, "Iteration " << iter )

    // Kaczmarz implementation start

    for ( IndexType iRow = 0; iRow < A.getNumRows(); ++iRow )
    {
        SCAI_LOG_DEBUG( logger, "Update irow " << iRow )

        A.getRow( z, iRow );

        z.redistribute( x.getDistributionPtr() );

        ValueType bi = b[ iRow ];

        ValueType p = ( bi - z.dotProduct( x ) ) / z.dotProduct( z );

        x = x + p * z;
    }


    // Kaczmarz implementation end

    mKaczmarzRuntime.mSolution.setDirty( true );
}

/* ========================================================================= */
/*       Runtime getter                                                      */
/* ========================================================================= */

template<typename ValueType>
typename Kaczmarz<ValueType>::KaczmarzRuntime& Kaczmarz<ValueType>::getRuntime()
{
    return mKaczmarzRuntime;
}

template<typename ValueType>
const typename Kaczmarz<ValueType>::KaczmarzRuntime& Kaczmarz<ValueType>::getRuntime() const
{
    return mKaczmarzRuntime;
}

/* ========================================================================= */
/*       Virtual methods                                                     */
/* ========================================================================= */

template<typename ValueType>
Kaczmarz<ValueType>* Kaczmarz<ValueType>::copy()
{
    return new Kaczmarz( *this );
}

template<typename ValueType>
void Kaczmarz<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "Kaczmarz<" << common::TypeTraits<ValueType>::id() << "> ( id = " << Solver<ValueType>::getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Kaczmarz, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
