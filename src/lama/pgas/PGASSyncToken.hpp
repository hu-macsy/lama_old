/**
 * @file PGASSyncToken.hpp
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
 * @brief PGASSyncToken.hpp
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */
#ifndef LAMA_PGAS_SYNC_TOKEN_HPP_
#define LAMA_PGAS_SYNC_TOKEN_HPP_

#include <lama/SyncToken.hpp>
#include <lama/LAMATypes.hpp>
#include <boost/scoped_array.hpp>

namespace lama
{

/** Class for PGAS synchronization that waits on pending messages. */

class PGASSyncToken: public SyncToken
{
public:

    /** Constructor for an PGAS synchronization token. */

    PGASSyncToken();

    /** Destructor, will also wait for synchronization and cleanup. */

    virtual ~PGASSyncToken();

    /** This method waits for all requests. */

    virtual void wait();

    /** Method to check whether communications have already been finished. */

    virtual bool probe() const;

    virtual void writeAt( std::ostream& stream ) const;

};

}

#endif // LAMA_PGAS_SYNC_TOKEN_HPP_
