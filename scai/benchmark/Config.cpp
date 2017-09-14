/**
 * @file Config.cpp
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
 * @brief Config.cpp
 * @author : robin
 * @date 06.04.2011
 */
/*
 * @file Config.cpp
 * @author: robin
 * Created on: 03.08.2010
 */

#include <scai/benchmark/Config.hpp>

namespace scai
{

namespace bf
{

Config& Config::getInstance()
{
    static Config config;
    return config;
}

const std::string& Config::getValueOf( const std::string& param )
{
    return m_params.operator[]( param );
}

void Config::setValueFor( std::string param, std::string value )
{
    m_params[param] = value;
}

Config::~Config()
{
}

Config::Config()
{
}

}

}
