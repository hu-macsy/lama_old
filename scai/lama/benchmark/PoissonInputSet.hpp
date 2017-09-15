/**
 * @file PoissonInputSet.hpp
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
 * @brief PoissonInputSet.hpp
 * @author Thomas Brandes
 * @date 15.09.2017
 */

#pragma once

#include <scai/benchmark.hpp>
#include <scai/lama/benchmark/LAMAInputSet.hpp>

#include <string>

namespace scai
{

namespace lama
{

/** Derived class of LAMAInputSet that creates the data by a stencil matrix */

class PoissonInputSet: public  LAMAInputSet, 
                       private benchmark::InputSet::Register<PoissonInputSet>
{
public:

    /** Constructor of an input set with poisson stencil matrix */

    PoissonInputSet( const std::string argument );

    /** Unique id for this derived class, used for InputSet::Factory to create objects */

    static std::string createValue();

    /** Static method required for create, called by InputSet::create( PoissonInputSet::createValue(), argument ) */

    static InputSet* create( const std::string argument );
};

}

}
