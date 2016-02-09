/**
 * @file SOR.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * SOFTWARE.
 * @endlicense
 *
 * @brief SOR.cpp
 * @author Martin Schenk
 * @date 29.08.2011
 * @since 1.0.0
 */

// hpp
#include <scai/solver/SOR.hpp>

// local library
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

// tracing
#include <scai/tracing.hpp>

using namespace scai::hmemo;

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( SOR::logger, "Solver.IterativeSolver.OmegaSolver.SOR" )

using lama::Matrix;
using lama::SparseMatrix;
using lama::CSRSparseMatrix;
using lama::Vector;
using lama::DenseVector;
using lama::Scalar;
using lama::CSRStorage;

SOR::SOR( const std::string& id )
    :
//constructor for GaussSeidel (SOR with omega==1)
    OmegaSolver( id, 1.0 )
{
}

//TODO check for 0 < omega < 2
SOR::SOR( const std::string& id, const Scalar omega )
    :
//constructor for SOR with relaxation factor omega
    OmegaSolver( id, omega )
{
}

SOR::SOR( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id, 1.0, logger )
{
}

SOR::SOR( const std::string& id, const Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger )
{
}

SOR::SOR( const SOR& other )
    : OmegaSolver( other )
{
}

SOR::SORRuntime::SORRuntime()
    : OmegaSolverRuntime()
{
}

SOR::~SOR()
{
}

SOR::SORRuntime::~SORRuntime()
{
}

void SOR::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "Solver.SOR.initialize" )

    Solver::initialize( coefficients );

    // ToDo: handling with other matrix types

    if( coefficients.getValueType() == common::scalar::FLOAT )
    {
        const lama::SparseMatrix<float>* A = dynamic_cast<const lama::SparseMatrix<float>*>( &coefficients );

        if( A && A->getLocalStorage().getFormat() == lama::Format::CSR && A->getHaloStorage().getFormat() == lama::Format::CSR )
        {
            return;
        }
        else
        {
            mIterationMatrix.reset( new CSRSparseMatrix<float>( coefficients ) );
            SCAI_LOG_INFO( logger, "conversion of iteration matrix to CSR: " << *mIterationMatrix )
        }
    }
    else if( coefficients.getValueType() == common::scalar::DOUBLE )
    {
        const SparseMatrix<double>* A = dynamic_cast<const SparseMatrix<double>*>( &coefficients );

        if( A && A->getLocalStorage().getFormat() == lama::Format::CSR && A->getHaloStorage().getFormat() == lama::Format::CSR )
        {
            return;
        }
        else
        {
            mIterationMatrix.reset( new CSRSparseMatrix<double>( coefficients ) );
            SCAI_LOG_INFO( logger, "conversion of iteration matrix to CSR: " << *mIterationMatrix )
        }
    }
    else
    {
        COMMON_THROWEXCEPTION(
            "Coefficients matrix " << typeid( coefficients ).name() << "(" << coefficients << ") is of unsupported type for SOR." );
    }
}

void SOR::iterate()
{
    SCAI_REGION( "Solver.SOR.iterate" )

    switch( getRuntime().mCoefficients->getValueType() )
    {
        case common::scalar::FLOAT:
        {
            iterateImpl<float>();
            break;
        }

        case common::scalar::DOUBLE:
        {
            iterateImpl<double>();
            break;
        }

        default:
        {
            COMMON_THROWEXCEPTION( "Value type " << getRuntime().mCoefficients->getValueType() << " is not implement." )
        }
    }
}

template<typename ValueType>
void SOR::iterateImpl()
{
    DenseVector<ValueType> & x = dynamic_cast<DenseVector<ValueType>&>( *getRuntime().mSolution );
    const DenseVector<ValueType> & b = dynamic_cast<const DenseVector<ValueType>&>( *getRuntime().mRhs );
    const Matrix* matrixPtr = getRuntime().mCoefficients;

    if( mIterationMatrix.get() )
    {
        SCAI_LOG_DEBUG( logger,
                        "Taking CSR converted matrix " << *mIterationMatrix << " instead of " << *getRuntime().mCoefficients )
        matrixPtr = mIterationMatrix.get();
    }

    const SparseMatrix<ValueType> & A = dynamic_cast<const SparseMatrix<ValueType>&>( *matrixPtr );
    const CSRStorage<ValueType> & csrA = dynamic_cast<const CSRStorage<ValueType>&>( A.getLocalStorage() );
    SCAI_ASSERT_ERROR( csrA.hasDiagonalProperty(), "csrA = " << csrA << " has not diagonal property" )
    const CSRStorage<ValueType> & csrAHalo = dynamic_cast<const CSRStorage<ValueType>&>( A.getHaloStorage() );
    const dmemo::Halo& Ahalo = A.getHalo();

    if( !Ahalo.isEmpty() )
    {
        x.updateHalo( Ahalo );
    }

    ContextPtr hostContext = Context::getHostPtr();

    WriteAccess<ValueType> xAcc( x.getLocalValues(), hostContext );
    ReadAccess<ValueType> bAcc( b.getLocalValues(), hostContext );
    ReadAccess<ValueType> aValues( csrA.getValues(), hostContext );
    ReadAccess<IndexType> ia( csrA.getIA(), hostContext );
    ReadAccess<IndexType> ja( csrA.getJA(), hostContext );

    // Getting access to arrays is safe even if not needed

    ReadAccess<ValueType> aHaloValues( csrAHalo.getValues(), hostContext );
    ReadAccess<IndexType> haloIa( csrAHalo.getIA(), hostContext );
    ReadAccess<IndexType> haloJa( csrAHalo.getJA(), hostContext );
    ReadAccess<ValueType> xHalo( x.getHaloValues(), hostContext );

    ValueType omega = mOmega.getValue<ValueType>();

    SCAI_LOG_DEBUG( logger, "iterate, omega = " << omega << ", #rows = " << csrA.getNumRows() )

    //SOR with relaxation factor omega
    if( omega != 1.0 )
    {
        SCAI_REGION( "Solver.SOR.iterate:Relaxation" )

        const ValueType oneMinusOmega = 1 - omega;

        for( IndexType i = 0; i < csrA.getNumRows(); i++ )
        {
            ValueType sum = 0;
            const ValueType diag = aValues[ia[i]];

            //Local loop
            for( IndexType pos = ia[i] + 1; pos < ia[i + 1]; pos++ )
            {
                sum = sum + xAcc[ja[pos]] * aValues[pos];
            }

            // Halo loop

            if( Ahalo.getHaloSize() > 0 )
            {

                for( IndexType pos = haloIa[i]; pos < haloIa[i + 1]; pos++ )
                {
                    sum = sum + xHalo[haloJa[pos]] * aHaloValues[pos];
                }
            }

            xAcc[i] = oneMinusOmega * xAcc[i] + omega * ( bAcc[i] - sum ) / diag;
        }

    }
    else //GaussSeidel (SOR without relaxation factor omega)
    {
        SCAI_REGION( "Solver.SOR.iterate:GaussSeidel" )

        for( IndexType i = 0; i < csrA.getNumRows(); i++ )
        {
            ValueType sum = 0;
            const ValueType diag = aValues[ia[i]];

            for( IndexType pos = ia[i] + 1; pos < ia[i + 1]; pos++ )
            {
                sum = sum + xAcc[ja[pos]] * aValues[pos];
            }

            // Halo loop

            if( Ahalo.getHaloSize() > 0 )
            {
                for( IndexType pos = haloIa[i]; pos < haloIa[i + 1]; pos++ )
                {
                    sum = sum + xHalo[haloJa[pos]] * aHaloValues[pos];
                }
            }

            xAcc[i] = ( bAcc[i] - sum ) / diag;
        }
    }

}

SOR::SORRuntime& SOR::getRuntime()
{
    return mSORRuntime;
}

const SOR::SORRuntime& SOR::getConstRuntime() const
{
    return mSORRuntime;
}

SolverPtr SOR::copy()
{
    return SolverPtr( new SOR( *this ) );
}

void SOR::writeAt( std::ostream& stream ) const
{
    stream << "SOR ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string SOR::createValue()
{
	return "SOR";
}

Solver* SOR::create( const std::string name )
{
	return new SOR( name );
}

} /* end namespace solver */

} /* end namespace scai */
