/**
 * @file utility.hpp
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
 * @brief Some utilities for parallel I/O
 * @author Thomas Brandes
 * @date 13.07.2016
 */

#include <scai/common/Settings.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>

using namespace scai;

/** Help routine to read a scalar value type --SCAI_TYPE=... */

static common::scalar::ScalarType getType() 
{
    common::scalar::ScalarType type = common::TypeTraits<double>::stype;
    
    std::string val;
    
    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {   
        scai::common::scalar::ScalarType env_type = scai::common::str2ScalarType( val.c_str() );
        
        if ( env_type == scai::common::scalar::UNKNOWN )
        {   
            std::cout << "SCAI_TYPE=" << val << " illegal, is not a scalar type" << std::endl;
        }
        
        type = env_type;
    }

    return type;
}
