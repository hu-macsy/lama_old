/**
 * @file PGASMaxLocFunctor.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief PGASMaxLocFunctor.hpp
 * @author Michael Drost
 * @date 02.02.2012
 * $Id$
 */
#ifndef LAMA_PGASMAXLOCFUNCTOR_HPP_
#define LAMA_PGASMAXLOCFUNCTOR_HPP_

#include <lama/pgas/PGASInterface.hpp>
#include <lama/pgas/PGASFunctor.hpp>

namespace lama
{
template<typename T>
class PGASMaxLocFunctor: public PGASFunctor
{
    T* mResult;
    T* mWork;
    int* mWorkRank;
    int* mResultRank;
    const PGASInterface* const mInterface;

public:
    PGASMaxLocFunctor( T val );
    virtual ~PGASMaxLocFunctor();
    virtual void iteration( int partner, bool active );
    T getResult();
    int getLoc();
};

template<typename T>
PGASMaxLocFunctor<T>::PGASMaxLocFunctor( T val )
    : mInterface( PGASInterface::getInstance() )
{
    mResult = static_cast<T*>( mInterface->allocate( sizeof(T) ) );
    mWork = static_cast<T*>( mInterface->allocate( sizeof(T) ) );
    mResultRank = static_cast<int*>( mInterface->allocate( sizeof(int) ) );
    *mWork = val;
    *mResult = val;
    *mResultRank = mInterface->getRank();
}
template<typename T>
PGASMaxLocFunctor<T>::~PGASMaxLocFunctor()
{
    mInterface->free( mWork, sizeof(T) );
    mInterface->free( mResult, sizeof(T) );
    mInterface->free( mResultRank, sizeof(int) );
}

template<typename T>
T PGASMaxLocFunctor<T>::getResult()
{
    return *mResult;
}

template<typename T>
int PGASMaxLocFunctor<T>::getLoc()
{
    return *mResultRank;
}

template<typename T>
void PGASMaxLocFunctor<T>::iteration( int partner, bool active )
{
    mInterface->syncronizeAll();
    *mWork = *mResult;
    mInterface->syncronizeAll();
    if ( active )
    {
        mInterface->get( mWork, mResult, sizeof(T), partner );

        if ( ( ( *mWork ) > ( *mResult ) ) )
        {
            *mResult = *mWork;
            mInterface->get( mResultRank, mResultRank, sizeof(int), partner );
        }
    }
}
} //namespace lama

#endif // LAMA_PGASMAXLOCFUNCTOR_HPP_
