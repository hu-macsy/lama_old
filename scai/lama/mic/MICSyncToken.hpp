/**
 * @file MICSyncToken.hpp
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
 * @brief MICSyncToken.hpp
 * @author Thomas Brandes
 * @date 02.07.2013
 * @since 1.1.0
 */

#pragma once

// base classes
#include <scai/lama/SyncToken.hpp>
#include <scai/lama/mic/MICContext.hpp>

// boost
#include <scai/common/shared_ptr.hpp>

namespace scai
{

namespace lama
{

typedef common::shared_ptr<const MICContext> MICContextPtr;

/** Class that sycnchronizes with a MIC offload transfer or computation. */

class COMMON_DLL_IMPORTEXPORT MICSyncToken: public scai::lama::SyncToken

{
public:

    /** Constructor for a MIC sychronization token.
     *
     *  @param[in]  context  is the MICcontext where asynchronous operation takes place
     *
     *  A pointer to the MIC context is required to enable/disable it.
     */

    MICSyncToken( MICContextPtr context );

    virtual ~MICSyncToken();

    /** After starting the offload computation/transfer with a signal this signal is set here. */

    void setSignal( int signal );

    virtual void wait();

    virtual bool probe() const;

private:

    MICContextPtr mMICContext; // needed for synchronization

    int mSignal; // set by an offload computation
};

} /* end namespace lama */

} /* end namespace scai */
