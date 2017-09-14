/**
 * @file InputSet.hpp
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
 * @brief InputSet.h
 * @author jiri
 * @date 06.04.2011
 */
/**
 * @file InputSet.h
 * @author jiri
 * Created on: 04.05.2010
 */
#pragma once

#include <map>
#include <string>

#include <scai/common/config.hpp>

namespace scai
{

namespace bf
{

class COMMON_DLL_IMPORTEXPORT InputSet
{
public:
    typedef std::map<std::string,unsigned long> ComplexityMap;
    /**
     * @brief Creates this object with the given id.
     * @param[in] id The unique id of this InputSet.
     */
    InputSet( const std::string& id );
    /**
     * @brief Creates this object with the given id and the given name.
     * @param[in] id    The unique id of this InputSet.
     * @param[in] name  The name of this InputSet.
     */
    InputSet( const std::string& id, const std::string& name );
    /**
     * @brief Default Destructor. Frees or destroys all inner ressources, if
     *        necessary.
     */
    virtual ~InputSet();

    /**
     * @brief Returns the id of the InputSet.
     * @return The unique id of the InputSet.
     */
    const std::string& getId() const;

    /**
     * @brief Returns the name of the InputSet.
     * @return The name of the InputSet.
     */
    const std::string& getName() const;

    /**
     * @brief Returns the number of floating point operations of the group with
     *        the given id.
     * @param[in] gid The id of the group for which the number of floating point
     *            operations will be returned.
     * @return The number of floating point operations of the group with the
     *         given id.
     */
    virtual unsigned long getNumFloatingPointOperations( const std::string& gid ) const;
    /**
     * @brief Sets the number of floating point operations for a group with the
     *        given id.
     * @param[in] gid       The id of the group for which the number will be set.
     * @param[in] numFlops  The number of floating point operations.
     */
    virtual void setNumFloatingPointOperations( const std::string& gid, const unsigned long numFlops );

    /**
     * @brief Returns the number of processed bytes for a given id.
     * @param[in] gid               The group id.
     * @param[in] sizeOfValueType   The size of the type of values.
     * @return The number of processed bytes.
     */
    virtual unsigned long getProcessedBytes( const std::string& gid, const unsigned short sizeOfValueType ) const;
    /**
     * @brief Sets the number of processed bytes for a group id and a certain
     *        type of values.
     * @param[in] gid               The id of the group.
     * @param[in] sizeOfValueType   The size of the type of the values.
     * @param[in] numBytes          The number of processed bytes.
     */
    virtual void setProcessedBytes(
        const std::string& gid,
        const unsigned short sizeOfValueType,
        const unsigned long numBytes );

private:
    InputSet( const InputSet& other );
    InputSet& operator=( const InputSet& other );

    const std::string mId;
    const std::string mName;

    ComplexityMap mFlopMap;
    std::map<unsigned short,ComplexityMap> mBWMap;
};

} // namespace bf

} // namespace scai
