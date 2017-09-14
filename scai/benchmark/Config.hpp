/**
 * @file Config.hpp
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
 * @brief Config.h
 * @author : robin
 * @date 06.04.2011
 */
/*
 * @file Config.h
 * @author: robin
 * Created on: 03.08.2010
 */

#pragma once

#include <string>
#include <map>

#include <scai/common/config.hpp>

namespace bf
{

/**
 * @brief This class holds static parameters for the benchmarks, which can be
 *        set by a main program, and asked by each benchmark, in case it is not
 *        possible to give this values as parameter values to the benchmarks.
 *        For example, this class could hold the path to input data.
 *        This class is a singleton, so there does exist only one instance of
 *        this class.
 */
class COMMON_DLL_IMPORTEXPORT Config
{
public:
    virtual ~Config();
    /**
     * @brief returns the onliest instance of this class.
     *
     * @return the onliest instance of this class.
     */
    static Config& getInstance();

    /**
     * @brief returns the value of the given parameter.
     *
     * @param[in] param the parameter, to get the value of.
     *
     * @return A constant reference to the value of the parameter.
     */
    const std::string& getValueOf( const std::string& param );

    /**
     * @brief sets the given value for the given parameter.
     *
     * Sets the given value for the given parameter, in case the given
     * parameter does not exist, yet.
     *
     * @param[in]   param   the parameter to be set.
     * @param[in]   value   the value of the parameter.
     */
    void setValueFor( std::string param, std::string value );

private:
    Config();
    Config( const Config& cc );
    std::map<std::string,std::string> m_params;
};

} // namespace bf
