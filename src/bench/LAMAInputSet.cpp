/**
 * @file LAMAInputSet.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief LAMAInputSet.cpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 12.09.2011
 * $Id$
 */

#include <bench/LAMAInputSet.hpp>

LAMA_LOG_DEF_LOGGER( LAMAInputSet::logger, "InputSet.LAMAInputSet" );

LAMAInputSet::LAMAInputSet( const std::string& id )
    : InputSet( id ), mAlpha( 0.0 ), mBeta( 0.0 ), mX( 0 ), mY( 0 ), mA( 0 ), mB( 0 ), mC( 0 ), mDenseA(
        0 ), mDenseB( 0 ), mDenseC( 0 )
{
    LAMA_LOG_INFO( logger, "LAMA Input set " << id << " created, all elems are NULL" );
}

LAMAInputSet::LAMAInputSet(
    const std::string& id,
    double alpha,
    double beta,
    std::auto_ptr<lama::DenseVector<double> > x,
    std::auto_ptr<lama::CSRSparseMatrix<double> > A )
    : InputSet( id ), mAlpha( alpha ), mBeta( beta ), mX( x.release() ), mY( mX ), mA( A.release() ), mB(
        mA ), mC( mA ), mDenseA( 0 ), mDenseB( 0 ), mDenseC( 0 )
{
    if( mX && mA )
    {
        LAMA_LOG_INFO( logger,
                       "LAMA Input set " << id << " created, " << ", alpha = " << alpha << ", beta = " << beta << ", x = " << *mX << ", A = " << *mA );
    }
    else
    {
        if( !mX )
        {
            LAMA_LOG_WARN( logger, "LAMA Input set with x = NULL created" );
        }
        if( !mA )
        {
            LAMA_LOG_WARN( logger, "LAMA Input set with A = NULL created" );
        }
    }
}

LAMAInputSet::LAMAInputSet(
    const std::string& id,
    double alpha,
    double beta,
    std::auto_ptr<lama::DenseVector<double> > x,
    std::auto_ptr<lama::DenseVector<double> > y,
    std::auto_ptr<lama::CSRSparseMatrix<double> > A )
    : InputSet( id ), mAlpha( alpha ), mBeta( beta ), mX( x.release() ), mY( y.release() ), mA(
        A.release() ), mB( mA ), mC( mA ), mDenseA( 0 ), mDenseB( 0 ), mDenseC( 0 )
{
    LAMA_LOG_INFO( logger,
                   "LAMA Input set " << id << " created: " << "alpha = " << alpha << ", beta = " << beta << ", x = " << *mX << ", y = " << *mY << ", A = " << *mA );
}

LAMAInputSet::LAMAInputSet(
    const std::string& id,
    double alpha,
    double beta,
    std::auto_ptr<lama::DenseVector<double> > x,
    std::auto_ptr<lama::DenseVector<double> > y,
    std::auto_ptr<lama::DenseMatrix<double> > A )
    : InputSet( id ), mAlpha( alpha ), mBeta( beta ), mX( x.release() ), mY( y.release() ), mA( 0 ), mB(
        0 ), mC( 0 ), mDenseA( A.release() ), mDenseB( mDenseA ), mDenseC( mDenseA )
{
    LAMA_LOG_INFO( logger,
                   "LAMA Input set " << id << " created: " << "alpha = " << alpha << ", beta = " << beta << ", x = " << *mX << ", y = " << *mY << ", DenseA = " << mDenseA );
}

LAMAInputSet::~LAMAInputSet()
{
    LAMA_LOG_INFO( logger, "~LAMAInputSet, will free input data" );

    if( mX == mY )
    {
        mY = 0;
    }
    if( mX != 0 )
    {
        delete mX;
    }
    mX = 0;
    if( mY != 0 )
    {
        delete mY;
    }
    mY = 0;
    if( mC == mA || mC == mB )
    {
        mC = 0;
    }
    if( mB == mA )
    {
        mB = 0;
    }
    if( mA != 0 )
    {
        delete mA;
    }
    mA = 0;
    if( mB != 0 )
    {
        delete mB;
    }
    mB = 0;
    if( mC != 0 )
    {
        delete mC;
    }
    mC = 0;

    if( mDenseC == mDenseA || mDenseC == mDenseB )
    {
        mDenseC = 0;
    }
    if( mDenseB == mDenseA )
    {
        mDenseB = 0;
    }
    if( mDenseA != 0 )
    {
        delete mDenseA;
    }
    mDenseA = 0;
    if( mDenseB != 0 )
    {
        delete mDenseB;
    }
    mDenseB = 0;
    if( mDenseC != 0 )
    {
        delete mDenseC;
    }
    mDenseC = 0;
}

double LAMAInputSet::getAlpha() const
{
    return mAlpha;
}

double LAMAInputSet::getBeta() const
{
    return mBeta;
}

const lama::DenseVector<double>& LAMAInputSet::getX() const
{
    return *mX;
}

const lama::DenseVector<double>& LAMAInputSet::getY() const
{
    return *mY;
}

const lama::CSRSparseMatrix<double>& LAMAInputSet::getA() const
{
    if( mA == 0 )
    {
        if( mDenseA == mDenseB && mB != 0 )
        {
            mA = mB;
        }
        else if( mDenseA == mDenseC && mC != 0 )
        {
            mA = mC;
        }
        else
        {
            mA = new lama::CSRSparseMatrix<double>( *mDenseA );
        }
    }
    return *mA;
}

const lama::CSRSparseMatrix<double>& LAMAInputSet::getB() const
{
    if( mB == 0 )
    {
        if( mDenseB == mDenseA && mA != 0 )
        {
            mB = mA;
        }
        else if( mDenseB == mDenseC && mC != 0 )
        {
            mB = mC;
        }
        else
        {
            mB = new lama::CSRSparseMatrix<double>( *mDenseB );
        }
    }
    return *mB;
}

const lama::CSRSparseMatrix<double>& LAMAInputSet::getC() const
{
    if( mC == 0 )
    {
        if( mDenseC == mDenseA && mA != 0 )
        {
            mC = mA;
        }
        else if( mDenseC == mDenseB && mB != 0 )
        {
            mC = mB;
        }
        else
        {
            mC = new lama::CSRSparseMatrix<double>( *mDenseC );
        }
    }
    return *mC;
}

const lama::DenseMatrix<double>& LAMAInputSet::getDenseA() const
{
    if( mDenseA == 0 )
    {
        if( mA == mB && mDenseB != 0 )
        {
            mDenseA = mDenseB;
        }
        else if( mA == mC && mDenseC != 0 )
        {
            mDenseA = mDenseC;
        }
        else
        {
            mDenseA = new lama::DenseMatrix<double>( *mA );
        }
    }
    return *mDenseA;
}

const lama::DenseMatrix<double>& LAMAInputSet::getDenseB() const
{
    if( mDenseB == 0 )
    {
        if( mB == mA && mDenseA != 0 )
        {
            mDenseB = mDenseA;
        }
        else if( mB == mC && mDenseC != 0 )
        {
            mDenseB = mDenseC;
        }
        else
        {
            mDenseB = new lama::DenseMatrix<double>( *mB );
        }
    }
    return *mDenseB;
}

const lama::DenseMatrix<double>& LAMAInputSet::getDenseC() const
{
    if( mDenseC == 0 )
    {
        if( mC == mA && mDenseA != 0 )
        {
            mDenseC = mDenseA;
        }
        else if( mC == mB && mDenseB != 0 )
        {
            mDenseC = mDenseB;
        }
        else
        {
            mDenseC = new lama::DenseMatrix<double>( *mC );
        }
    }
    return *mDenseC;
}
