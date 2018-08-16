/**
 * @file Jacobi.cpp
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
 * @brief Implementation of class Jacobi.
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/Jacobi.hpp>

// local library
#include <scai/utilskernel/HArrayUtils.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/instantiate.hpp>

#include <functional>

using std::function;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;

using std::bind;

namespace scai
{

using tasking::SyncToken;

using utilskernel::HArrayUtils;
using lama::Matrix;
using lama::MatrixKind;
using lama::SparseMatrix;
using lama::VectorKind;
using lama::DenseVector;
using lama::Vector;

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, Jacobi<ValueType>::logger, "Solver.IterativeSolver.Jacobi" )

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
Jacobi<ValueType>::Jacobi( const std::string& id ) : 

    OmegaSolver<ValueType>( id )
{
}

template<typename ValueType>
Jacobi<ValueType>::Jacobi( const std::string& id, LoggerPtr logger ) : 

    OmegaSolver<ValueType>( id, logger )
{
}

template<typename ValueType>
Jacobi<ValueType>::Jacobi( const std::string& id, ValueType omega ) : 

    OmegaSolver<ValueType>( id, omega )
{
}

template<typename ValueType>
Jacobi<ValueType>::Jacobi( const std::string& id, ValueType omega, LoggerPtr logger ) : 
 
    OmegaSolver<ValueType>( id, omega, logger )
{
}

template<typename ValueType>
Jacobi<ValueType>::Jacobi( const Jacobi& other ) : 

    OmegaSolver<ValueType>( other )
{
}

template<typename ValueType>
Jacobi<ValueType>::~Jacobi()
{
    SCAI_LOG_INFO( logger, "~Jacobi" )
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void Jacobi<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    OmegaSolver<ValueType>::initialize( coefficients );

    if ( coefficients.getMatrixKind() != MatrixKind::SPARSE )
    {
        COMMON_THROWEXCEPTION( "ERROR: Jacobi can only be applied for sparse matrices, but here: " << coefficients )
    }

    // Additional stuff: inherit context of matrix to runtime data, extract diagonal 

    SCAI_LOG_DEBUG( logger, "Initialization started" )

    JacobiRuntime& runtime = getRuntime();

    DenseVector<ValueType>& diagonal    = runtime.mDiagonal;
    DenseVector<ValueType>& oldSolution = runtime.mOldSolution;

    hmemo::ContextPtr ctx = coefficients.getContextPtr();

    diagonal.setContextPtr( ctx );
    oldSolution.setContextPtr( ctx );

    coefficients.getDiagonal( diagonal );  // diagonal has correct distribution

    SCAI_LOG_INFO( logger, "Jacobi initialized, diagonal = " << runtime.mDiagonal
                            << ", oldSolution = " << runtime.mOldSolution )
}

template<typename ValueType>
void Jacobi<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    IterativeSolver<ValueType>::solveInit( solution, rhs );

    if ( VectorKind::DENSE != solution.getVectorKind() )
    {
        COMMON_THROWEXCEPTION( "Jacobi solver requires dense vector for solution." )
    }
    if ( VectorKind::DENSE != rhs.getVectorKind() )
    {
        COMMON_THROWEXCEPTION( "Jacobi solver requires dense vector for rhs." )
    }

    // Importan: allocate old solution with same dist as solution

    getRuntime().mOldSolution.allocate( solution.getDistributionPtr() );

    SCAI_LOG_INFO( logger, "solveInit, solution = " << getRuntime().mSolution.getConstReference() 
                          << ", rhs = " << *getRuntime().mRhs )
}

template<typename ValueType>
void Jacobi<ValueType>::iterate()
{
    SCAI_REGION( "Solver.SpJacobi.iterate" )

    using hmemo::HArray;

    JacobiRuntime& runtime = getRuntime();  

    const Matrix<ValueType>& m = *getRuntime().mCoefficients;

    if ( m.getNumRows() == 0 )
    {
        SCAI_LOG_WARN( logger, "Zero sized matrix given. Won't execute any calculations in this iteration. " )
        return;
    }

    Vector<ValueType>& solutionV = runtime.mSolution.getReference();   // mark solution as dirty

    // Note: start casts are safe as already verified in initialize, solveInit

    const SparseMatrix<ValueType>& coefficients = static_cast<const SparseMatrix<ValueType>&>( m );
  
    DenseVector<ValueType>& oldSolution = runtime.mOldSolution;
    DenseVector<ValueType>& solution    = static_cast<DenseVector<ValueType>&>( solutionV );
    const DenseVector<ValueType>& rhs   = static_cast<const DenseVector<ValueType>&>( *runtime.mRhs );

    // Swap solution and old solution

    SCAI_ASSERT_EQ_ERROR( solution.getDistribution(), oldSolution.getDistribution(), "mismatch" )

    solution.swap( oldSolution );

    const ValueType omega = OmegaSolver<ValueType>::mOmega;

    // from rhs and solution we need only the local parts as LAMA arrays

    const HArray<ValueType>& localRhs = rhs.getLocalValues();
    HArray<ValueType>& localSolution  = solution.getLocalValues();

    HArray<ValueType>& localOldSolution = oldSolution.getLocalValues();
    HArray<ValueType>& haloOldSolution  = oldSolution.getHaloValues();

    const HArray<ValueType>& localDiagonal = runtime.mDiagonal.getLocalValues();

    // Here we use a general compute/communicate scheme that is exactly the same
    // as for matrix-vector multiplication but with other local/halo computations

    using std::bind;
    using std::cref;
    using namespace std::placeholders;  // placeholders are also needed

    void ( scai::lama::MatrixStorage<ValueType>::*jacobiIterateHalo )(
        HArray<ValueType>& localSolution,
        const HArray<ValueType>& localDiagonal,
        const HArray<ValueType>& oldHaloSolution,
        const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterateHalo;

    // will call jacobiIterateHalo( haloMatrix, localSolution, diagonal, haloOldSolution, omega )
    function <
    void(
        const lama::MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX ) > haloF =
            bind( jacobiIterateHalo, _1, _2, cref( localDiagonal ), _3, omega );

    if ( lama::SyncKind::SYNCHRONOUS == coefficients.getCommunicationKind() )
    {
        // For the local operation a jacobi step is done
        void ( lama::MatrixStorage<ValueType>::*jacobiIterate )(
            HArray<ValueType>& solution,
            const HArray<ValueType>& oldSolution,
            const HArray<ValueType>& rhs,
            const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterate;
        // Bind the additional arguments like localRhs and omega
        function <
        void(
            const lama::MatrixStorage<ValueType>* haloMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX ) > localF =
                bind( jacobiIterate, _1, _2, _3, cref( localRhs ), omega );

        coefficients.haloOperationSync( localSolution, localOldSolution, haloOldSolution, localF, haloF );
    }
    else
    {
        SyncToken* ( lama::MatrixStorage<ValueType>::*jacobiIterateAsync )(
            HArray<ValueType>& solution,
            const HArray<ValueType>& oldSolution,
            const HArray<ValueType>& rhs,
            const ValueType omega ) const = &lama::MatrixStorage<ValueType>::jacobiIterateAsync;

        function <
        SyncToken*(
            const lama::MatrixStorage<ValueType>* haloMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX ) > localAsyncF =
                bind( jacobiIterateAsync, _1, _2, _3, cref( localRhs ), omega );

        coefficients.haloOperationAsyncLocal( localSolution, localOldSolution, haloOldSolution, localAsyncF, haloF );
    }

    SCAI_LOG_INFO( logger, "Jacobi iterate done, local sol = " << localSolution 
                           << ", local old sol = " << localOldSolution )
}

template<typename ValueType>
typename Jacobi<ValueType>::JacobiRuntime& Jacobi<ValueType>::getRuntime()
{
    return mJacobiRuntime;
}

template<typename ValueType>
const typename Jacobi<ValueType>::JacobiRuntime& Jacobi<ValueType>::getRuntime() const
{
    return mJacobiRuntime;
}

template<typename ValueType>
Jacobi<ValueType>* Jacobi<ValueType>::copy()
{
    return new Jacobi( *this );
}

template<typename ValueType>
void Jacobi<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "Jacobi<" << common::TypeTraits<ValueType>::id() << "> ( id = " << this->getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
SolverCreateKeyType Jacobi<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "Jacobi" );
}


template<typename ValueType>
_Solver* Jacobi<ValueType>::create()
{
    return new Jacobi<ValueType>( "_genByFactory" );
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Jacobi, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

// template solver::_Solver::Register<solver::Jacobi<ValueType>>::RegisterGuard
//          solver::_Solver::Register<solver::Jacobi<ValueType>>::registerGuard;

} /* end namespace scai */
