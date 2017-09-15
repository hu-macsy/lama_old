/**
 * @file PoissonInputSetCreator.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief PoissonInputSetCreator.hpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 06.05.2010
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

namespace scai
{

namespace lama
{

/** This class creates 'distributed' input sets for poisson solvers.
 *
 *  \code
 *     common::unique_ptr<LAMAInputSet> inputSet( PoissonInputSetCreator::createSet( "2D_9P_4_4" ) );
 *  \endcode
 *
 *  The matrix A will have a general block distribution and is in CSR format.
 */

class PoissonInputSetCreator: public benchmark::InputSetCreator<LAMAInputSet>
{
public:

    typedef benchmark::InputSetCreator<LAMAInputSet>::InputSetType InputSetType;

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

}

}

