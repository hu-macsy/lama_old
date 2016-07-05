/**
 * @file solver/mepr/SolverEps.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Wrapper for templated calls in Solver
 * @author Eric Schricker
 * @date 14.04.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/ScalarType.hpp>

#include <scai/lama/Scalar.hpp>

namespace scai
{

namespace solver
{

namespace mepr
{

template<typename TList>
struct SolverEps;

template<>
struct SolverEps<common::mepr::NullType>
{
    static lama::Scalar eps1( const common::scalar::ScalarType& )
    {
        return lama::Scalar( 0 );
    }
 
    static lama::Scalar eps0( const common::scalar::ScalarType& )
    {
        return lama::Scalar( 0 );
    }
};

template<typename H, typename T>
struct SolverEps<common::mepr::TypeList<H, T> >
{
    static lama::Scalar eps1( const common::scalar::ScalarType& type )
    {
        if ( common::TypeTraits<H>::stype == type )
        {
            return lama::Scalar( common::TypeTraits<H>::eps1() );
        }
        else
        {
            return SolverEps<T>::eps1( type );
        }
    }

    static lama::Scalar eps0( const common::scalar::ScalarType& type )
    {
        if ( common::TypeTraits<H>::stype == type )
        {
            return lama::Scalar( common::TypeTraits<H>::eps0() );
        }
        else
        {
            return SolverEps<T>::eps0( type );
        }
    }
};

} /* end namespace mepr */
} /* end namespace solver */
} /* end namepsace scai */
