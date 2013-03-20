/**
 * @file PGASMaxFunctor.hpp
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
 * @brief PGASMaxFunctor.hpp
 * @author Michael Drost
 * @date 02.02.2012
 * $Id$
 */
#ifndef LAMA_PGASMAXFUNCTOR_HPP_
#define LAMA_PGASMAXFUNCTOR_HPP_

#include <lama/pgas/PGASInterface.hpp>
#include <lama/pgas/PGASFunctor.hpp>

namespace lama
{
template<typename T>
class PGASMaxFunctor: public PGASFunctor
{
    T* mResult;
    T* mWork;
    const PGASInterface* const mInterface;

public:
    PGASMaxFunctor( T val );
    virtual ~PGASMaxFunctor();
    virtual void iteration( int partner, bool active );
    T getResult();
};

template<typename T>
PGASMaxFunctor<T>::PGASMaxFunctor( T val )
    : mInterface( PGASInterface::getInstance() )
{
    mResult = static_cast<T*>( mInterface->allocate( sizeof(T) ) );
    mWork = static_cast<T*>( mInterface->allocate( sizeof(T) ) );
    *mWork = val;
    *mResult = val;
}
template<typename T>
PGASMaxFunctor<T>::~PGASMaxFunctor()
{
    mInterface->free( mWork, sizeof(T) );
    mInterface->free( mResult, sizeof(T) );
}

template<typename T>
T PGASMaxFunctor<T>::getResult()
{
    return *mResult;
}

template<typename T>
void PGASMaxFunctor<T>::iteration( int partner, bool active )
{
    mInterface->syncronizeAll();
    *mWork = *mResult;
    mInterface->syncronizeAll();
    if( active )
    {
        mInterface->get( mResult, mWork, sizeof(T), partner );
        *mResult = ( *mWork < *mResult ) ? *mResult : *mWork;
    }
}
} //namespace lama

#endif // LAMA_PGASMAXFUNCTOR_HPP_
