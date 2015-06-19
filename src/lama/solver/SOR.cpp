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
#include <lama/solver/SOR.hpp>

// others
#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>

// tracing
#include <tracing/tracing.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( SOR::logger, "Solver.IterativeSolver.OmegaSolver.SOR" )

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
    LAMA_REGION( "Solver.SOR.initialize" )

    Solver::initialize( coefficients );

    // ToDo: handling with other matrix types

    if( coefficients.getValueType() == Scalar::FLOAT )
    {
        const SparseMatrix<float>* A = dynamic_cast<const SparseMatrix<float>*>( &coefficients );

        if( A && A->getLocalStorage().getFormat() == Format::CSR && A->getHaloStorage().getFormat() == Format::CSR )
        {
            return;
        }
        else
        {
            mIterationMatrix.reset( new CSRSparseMatrix<float>( coefficients ) );
            LAMA_LOG_INFO( logger, "conversion of iteration matrix to CSR: " << *mIterationMatrix )
        }
    }
    else if( coefficients.getValueType() == Scalar::DOUBLE )
    {
        const SparseMatrix<double>* A = dynamic_cast<const SparseMatrix<double>*>( &coefficients );

        if( A && A->getLocalStorage().getFormat() == Format::CSR && A->getHaloStorage().getFormat() == Format::CSR )
        {
            return;
        }
        else
        {
            mIterationMatrix.reset( new CSRSparseMatrix<double>( coefficients ) );
            LAMA_LOG_INFO( logger, "conversion of iteration matrix to CSR: " << *mIterationMatrix )
        }
    }
    else
    {
        LAMA_THROWEXCEPTION(
            "Coefficients matrix " << typeid( coefficients ).name() << "(" << coefficients << ") is of unsupported type for SOR." );
    }
}

void SOR::iterate()
{
    LAMA_REGION( "Solver.SOR.iterate" )

    switch( getRuntime().mCoefficients->getValueType() )
    {
        case Scalar::FLOAT:
        {
            iterateImpl<float>();
            break;
        }

        case Scalar::DOUBLE:
        {
            iterateImpl<double>();
            break;
        }

        default:
        {
            LAMA_THROWEXCEPTION( "Value type " << getRuntime().mCoefficients->getValueType() << " is not implement." )
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
        LAMA_LOG_DEBUG( logger,
                        "Taking CSR converted matrix " << *mIterationMatrix << " instead of " << *getRuntime().mCoefficients )
        matrixPtr = mIterationMatrix.get();
    }

    const SparseMatrix<ValueType> & A = dynamic_cast<const SparseMatrix<ValueType>&>( *matrixPtr );
    const CSRStorage<ValueType> & csrA = dynamic_cast<const CSRStorage<ValueType>&>( A.getLocalStorage() );
    LAMA_ASSERT_ERROR( csrA.hasDiagonalProperty(), "csrA = " << csrA << " has not diagonal property" )
    const CSRStorage<ValueType> & csrAHalo = dynamic_cast<const CSRStorage<ValueType>&>( A.getHaloStorage() );
    const Halo& Ahalo = A.getHalo();

    if( !Ahalo.isEmpty() )
    {
        x.updateHalo( Ahalo );
    }

    HostWriteAccess<ValueType> xAcc( x.getLocalValues() );
    HostReadAccess<ValueType> bAcc( b.getLocalValues() );
    HostReadAccess<ValueType> aValues( csrA.getValues() );
    HostReadAccess<IndexType> ia( csrA.getIA() );
    HostReadAccess<IndexType> ja( csrA.getJA() );

    // Getting access to arrays is safe even if not needed

    HostReadAccess<ValueType> aHaloValues( csrAHalo.getValues() );
    HostReadAccess<IndexType> haloIa( csrAHalo.getIA() );
    HostReadAccess<IndexType> haloJa( csrAHalo.getJA() );
    HostReadAccess<ValueType> xHalo( x.getHaloValues() );

    ValueType omega = mOmega.getValue<ValueType>();

    LAMA_LOG_DEBUG( logger, "iterate, omega = " << omega << ", #rows = " << csrA.getNumRows() )

    //SOR with relaxation factor omega
    if( omega != 1.0 )
    {
        LAMA_REGION( "Solver.SOR.iterate:Relaxation" )

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
        LAMA_REGION( "Solver.SOR.iterate:GaussSeidel" )

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

} /* namespace lama */
