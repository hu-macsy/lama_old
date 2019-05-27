/**
 * @file AMGSetup.cpp
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
 * @brief AMGSetup.cpp
 * @author Jiri Kraus
 * @date 28.10.2011
 */

// hpp
#include <scai/solver/AMGSetup.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using lama::Matrix;
using lama::Vector;

namespace solver
{

SCAI_LOG_DEF_LOGGER( _AMGSetup::logger, "AMGSetup" )

/* ========================================================================= */
/*    static methods (for _AMGSetup - Factory, untyped )                     */
/* ========================================================================= */

_AMGSetup* _AMGSetup::getAMGSetup( const common::ScalarType scalarType, const std::string& setupType )
{
    return create( AMGSetupCreateKeyType( scalarType, setupType ) );
}

/* ========================================================================= */
/*    static methods (for AMGSetup<ValueType> - Factory                      */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::getCreateValues( std::vector<std::string>& values )
{
    std::vector<AMGSetupCreateKeyType> createValues;

    _AMGSetup::getCreateValues( createValues );  // all solvers ( valueType, solvertype )

    values.clear();

    for ( size_t i = 0; i < createValues.size(); ++i )
    {
        if ( createValues[i].first == common::TypeTraits<ValueType>::stype )
        {
            // AMGSetup for this value type
            values.push_back( createValues[i].second );
        }
    }
}

template<typename ValueType>
AMGSetup<ValueType>* AMGSetup<ValueType>::getAMGSetup( const std::string& setupType )
{
    _AMGSetup* setup = _AMGSetup::getAMGSetup( common::TypeTraits<ValueType>::stype, setupType );

    SCAI_ASSERT_DEBUG( dynamic_cast<AMGSetup<ValueType>*>( setup ), "Illegal setup" )

    return reinterpret_cast<AMGSetup<ValueType>*>( setup );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
AMGSetup<ValueType>::AMGSetup() : 

    mHostOnlyLevel( std::numeric_limits<IndexType>::max() ), 
    mHostOnlyVars( 0 ), 
    mReplicatedLevel( std::numeric_limits<IndexType>::max() ),
    mMainGalerkinMatrix( nullptr )
{
}

template<typename ValueType>
AMGSetup<ValueType>::~AMGSetup()
{
}

/* ========================================================================= */
/*    Methods                                                                */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::setHostOnlyLevel( IndexType hostOnlyLevel )
{
    mHostOnlyLevel = hostOnlyLevel;
}

template<typename ValueType>
void AMGSetup<ValueType>::setHostOnlyVars( IndexType hostOnlyVars )
{
    mHostOnlyLevel = hostOnlyVars;
}

template<typename ValueType>
void AMGSetup<ValueType>::setReplicatedLevel( IndexType replicatedLevel )
{
    mReplicatedLevel = replicatedLevel;
}

template<typename ValueType>
void AMGSetup<ValueType>::setMaxLevels( IndexType maxLevels )
{
    mMaxLevels = maxLevels;
}

template<typename ValueType>
void AMGSetup<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "AMGSetup( ... )";
}

/* ========================================================================= */
/*     initialize ( common part for all AMG Setup classes )                  */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::initialize( const Matrix<ValueType>& mainSystemMatrix )
{
    SCAI_REGION( "AMGSetup.initialize" )

    SCAI_ASSERT_EQ_ERROR( 0, getNumLevels(), "AMGSetup already initialized" )

    // pushback main system matrix to storage vector

    mMainGalerkinMatrix = &mainSystemMatrix;

    // clear the last matrix hierarchy

    mGalerkinMatrices.clear();
    mInterpolationMatrices.clear();
    mRestrictionMatrices.clear();

    createMatrixHierarchy();   // virtual function implemented individually by derived AMGSetup classes

    createVectorHierarchy();   // create objects for rhs, solution, residual on each level

    createSolverHierarchy();   // create smoothers and coarse level solver
}

/* ========================================================================= */
/*      Management of vectors for all AMG levels                             */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::createVectorHierarchy()
{
    SCAI_ASSERT_ERROR( mMainGalerkinMatrix, "null pointer for main matrix, setup not not initialized yet" )

    SCAI_REGION( "AMGSetup.createVectorHierarchy" );

    IndexType numLevels = getNumLevels();

    SCAI_LOG_INFO( logger, "Creating Vectors for hierarchy, #levels = " << numLevels );

    mSolutionHierarchy.resize( numLevels );
    mRhsHierarchy.resize( numLevels );
    mTmpResHierarchy.resize( numLevels );
}

template<typename ValueType>
Vector<ValueType>& AMGSetup<ValueType>::getSolutionVector( const IndexType level )
{
    IndexType numLevels = mSolutionHierarchy.size();

    SCAI_ASSERT_ERROR( level != 0, "SolutionVector on level 0 is not stored in hierarachy." );
    SCAI_ASSERT_VALID_INDEX_ERROR( level, numLevels, "SolutionVector on Level " << level << " does not exist" );

    return mSolutionHierarchy[level];
}

template<typename ValueType>
Vector<ValueType>& AMGSetup<ValueType>::getRhsVector( const IndexType level )
{
    IndexType numRhsVectors = mRhsHierarchy.size();

    SCAI_ASSERT_ERROR( level != 0, "RhsVector on level 0 is not stored in hierarachy." );
    SCAI_ASSERT_VALID_INDEX_ERROR( level, numRhsVectors, "illegal level for rhs vecotr" )
    return mRhsHierarchy[level];
}

template<typename ValueType>
Vector<ValueType>& AMGSetup<ValueType>::getTmpResVector( const IndexType level )
{
    IndexType numResVectors = mTmpResHierarchy.size();
    SCAI_ASSERT_VALID_INDEX_ERROR( level, numResVectors, "TmpResVector on Level " << level << " does not exist" );
    return mTmpResHierarchy[level];
}

/* ========================================================================= */
/*      Management of the matrices on / between all AMG levels               */
/* ========================================================================= */

template<typename ValueType>
IndexType AMGSetup<ValueType>::getNumLevels() const
{
    // Attention: main matrix counts as an additional level

    if ( mMainGalerkinMatrix == nullptr )
    {
        return 0;
    }
    else
    {
        return static_cast<IndexType>( mGalerkinMatrices.size() + 1 );
    }
}

template<typename ValueType>
void AMGSetup<ValueType>::addNextLevel( 
    std::unique_ptr<lama::Matrix<ValueType>> interpolationMatrix,
    std::unique_ptr<lama::Matrix<ValueType>> galerkinMatrix,
    std::unique_ptr<lama::Matrix<ValueType>> restrictionMatrix )
{
    SCAI_ASSERT_ERROR( interpolationMatrix.get(), "null pointer for interpolation matrix" )

    // check for a correct interpolation matrix, #rows must be same as size of last galerkin matrix

    const Matrix<ValueType>& lastGalerkin = getGalerkin( getNumLevels() - 1 );

    SCAI_ASSERT_EQ_ERROR( interpolationMatrix->getNumRows(), lastGalerkin.getNumRows(), 
                          "interpolation matrix for level does not match size of last Galerkin matrix." )

    if ( restrictionMatrix.get() )
    {
        // check for a correct restriction matrix
 
        SCAI_ASSERT_EQ_ERROR( restrictionMatrix->getNumRows(), interpolationMatrix->getNumColumns(), "serious mismatch" )
        SCAI_ASSERT_EQ_ERROR( restrictionMatrix->getNumColumns(), interpolationMatrix->getNumRows(), "serious mismatch" )

        mRestrictionMatrices.push_back( std::move( restrictionMatrix ) );
    }
    else
    {
        std::unique_ptr<Matrix<ValueType>> ownRestrictionMatrix( interpolationMatrix->newMatrix() );
        ownRestrictionMatrix->assignTranspose( *interpolationMatrix );
        mRestrictionMatrices.push_back( std::move( ownRestrictionMatrix ) );
    }

    if ( galerkinMatrix.get() )
    {
        // check galerkin matrix

        SCAI_ASSERT_EQ_ERROR( galerkinMatrix->getNumRows(), galerkinMatrix->getNumColumns(), "galerkin matrix not square" )
        SCAI_ASSERT_EQ_ERROR( galerkinMatrix->getNumRows(), interpolationMatrix->getNumColumns(), "galerkin matrix not square" )

        mGalerkinMatrices.push_back( std::move( galerkinMatrix ) );
    }
    else
    {
        // compute own galerkin matrix
        COMMON_THROWEXCEPTION( "compute of restriction matrix not supported yet" )
    }

    mInterpolationMatrices.push_back( std::move( interpolationMatrix ) );
}

/* ========================================================================= */

template<typename ValueType>
const Matrix<ValueType>& AMGSetup<ValueType>::getGalerkin( const IndexType level )
{
    IndexType numLevels = getNumLevels();

    SCAI_ASSERT_VALID_INDEX_ERROR( level, numLevels, "Galerkin on Level " << level << " does not exist" );

    if ( level == 0 )
    {
        return *mMainGalerkinMatrix;
    }
    else
    {
        return *mGalerkinMatrices[level-1];
    }
}

template<typename ValueType>
const Matrix<ValueType>& AMGSetup<ValueType>::getRestriction( const IndexType level )
{
    IndexType numRestrictionMatrices = mRestrictionMatrices.size();   // type conversion to LAMA index type

    SCAI_ASSERT_VALID_INDEX_ERROR( level, numRestrictionMatrices, "illegal get for restriction matrix" )

    return *mRestrictionMatrices[level];
}

template<typename ValueType>
const Matrix<ValueType>& AMGSetup<ValueType>::getInterpolation( const IndexType level )
{
    IndexType numInterpolationMatrices = mInterpolationMatrices.size();

    SCAI_ASSERT_VALID_INDEX_ERROR( level, numInterpolationMatrices, "illegal level for interpolatin" )
    return *mInterpolationMatrices[level];
}

/* ========================================================================= */
/*       Managment of Smoothers on all levels                                */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::setCoarseLevelSolver( SolverPtr<ValueType> solver )
{
    mCoarseLevelSolver = solver;
}

template<typename ValueType>
void AMGSetup<ValueType>::setSmoother( SolverPtr<ValueType> solver )
{
    SCAI_LOG_INFO( logger, "Set own smoother" )
    mSmoother = solver;
}

template<typename ValueType>
Solver<ValueType>& AMGSetup<ValueType>::getSmoother( const IndexType level )
{
    IndexType numLevels = mSolverHierarchy.size();

    SCAI_ASSERT_EQ_DEBUG( numLevels, getNumLevels(), "serious mismatch at smoother query" )

    SCAI_ASSERT_VALID_INDEX_ERROR( level, numLevels, "Smoother on Level " << level << " does not exist" );

    return *mSolverHierarchy[level];
}

template<typename ValueType>
Solver<ValueType>& AMGSetup<ValueType>::getCoarseLevelSolver()
{
    IndexType numLevels = mSolverHierarchy.size();

    SCAI_ASSERT_GT_ERROR( numLevels, 0, "no solvers set." )

    return *mSolverHierarchy[numLevels - 1];
}

template<typename ValueType>
std::string AMGSetup<ValueType>::getCoarseLevelSolverInfo() const
{
    IndexType numLevels = mSolverHierarchy.size();

    if ( numLevels > 0 )
    {
        return mSolverHierarchy.back()->getId();
    }
    else
    {
        return "No Coarse Solver";
    }
}

template<typename ValueType>
std::string AMGSetup<ValueType>::getSmootherInfo() const
{
    if ( getNumLevels() <= 1 )
    {
        return "No Smoother / 1-Level Approach";
    }
    else
    {
        return mSolverHierarchy[0]->getId();
    }
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::createSolverHierarchy()
{
    mSolverHierarchy.clear();    // just in case 

    SCAI_REGION( "AMGSetup.createSolverHierarchy" );

    IndexType numLevels = getNumLevels();

    SCAI_ASSERT_GE_ERROR( numLevels, 1, "no levels, cannot create solver" )

    SCAI_LOG_INFO( logger, "Initializing solver-hierarchy, numLevels = " << numLevels );

    for ( IndexType i = 0; i < numLevels; ++i )
    {
        bool isCoarseLevel = ( i + 1 ) == numLevels;  

        SolverPtr<ValueType> solver;

        if ( mSmoother && !isCoarseLevel )
        {
            solver.reset( mSmoother->copy() );
        }
        else if ( isCoarseLevel && mCoarseLevelSolver )
        {
            solver.reset( mCoarseLevelSolver->copy() );
        }
        else
        {
            // get the default solver as provided by the derived AMG setup class

            solver = createSolver( isCoarseLevel );
        }

        solver->initialize( getGalerkin( i ) );

        SCAI_LOG_INFO( logger, "\tLevel " << i << *solver )

        mSolverHierarchy.push_back( solver );
    }
}

/* ========================================================================= */

template<typename ValueType>
void AMGSetup<ValueType>::setMinVarsCoarseLevel( const IndexType vars )
{
    mMinVarsCoarseLevel = vars;
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( AMGSetup, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
