/**
 * @file PoissonInputSetCreator.hpp
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
 * @brief PoissonInputSetCreator.hpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 06.05.2010, revised 13.09.2011
 * $Id$
 */

#pragma once

#include <scai/benchmark.hpp>

#include <scai/lama/benchmark/PoissonInputSetCreator.hpp>
#include <scai/lama/benchmark/LAMAInputSet.hpp>

#include <scai/logging.hpp>

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>

/** This class creates 'distributed' input sets for poisson solvers.
 *
 *  \code
 *     common::unique_ptr<LAMAInputSet> inputSet( PoissonInputSetCreator::createSet( "2D_9P_4_4" ) );
 *  \endcode
 *
 *  The matrix A will have a general block distribution and is in CSR format.
 */

class PoissonInputSetCreator: public bf::InputSetCreator<LAMAInputSet>
{
public:

    typedef bf::InputSetCreator<LAMAInputSet>::InputSetType InputSetType;

    static const std::string& id();

    static InputSetType* createSet( const std::string& params );

    PoissonInputSetCreator();

    virtual ~PoissonInputSetCreator();

    virtual InputSetType* create() const;

    /** Implements pure method by using static method createSet. */

    virtual InputSetType* create( const std::string& params ) const;

    /** Implements pure method by using static method id. */

    virtual const std::string& getId() const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger );
};
