/**
 * @file NoSyncToken.hpp
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
 * @brief NoSyncToken.hpp
 * @author Thomas Brandes
 * @date 25.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_NOSYNCTOKEN_HPP_
#define LAMA_NOSYNCTOKEN_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/SyncToken.hpp>

namespace lama
{

class LAMA_DLL_IMPORTEXPORT NoSyncToken: public SyncToken
{
public:

    inline NoSyncToken();

    inline virtual ~NoSyncToken();

    inline virtual void wait();

    inline virtual bool probe() const;

    inline void writeAt( std::ostream& stream ) const;
};

NoSyncToken::NoSyncToken()
{
    LAMA_LOG_DEBUG( logger, "NoSyncToken constructed" )
}

NoSyncToken::~NoSyncToken()
{
    LAMA_LOG_DEBUG( logger, "~NoSyncToken, synchronized = " << isSynchronized() )

    if( !isSynchronized() )
    {
        setSynchronized(); // Important: accesses should be freed
    }
}

void NoSyncToken::writeAt( std::ostream& stream ) const
{
    stream << "NoSyncToken( synchronized = " << isSynchronized() << " )";
}

void NoSyncToken::wait()
{
    if( !isSynchronized() )
    {
        setSynchronized(); // Important: accesses should be freed
    }
}

bool NoSyncToken::probe() const
{
    return true; // always ready
}

}

#endif // LAMA_NOSYNCTOKEN_HPP_
