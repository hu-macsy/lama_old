/**
 * @file InputSet.h
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
 * @brief InputSet.h
 * @author jiri
 * @date 06.04.2011
 * $Id$
 */
/**
 * @file InputSet.h
 * @author jiri
 * Created on: 04.05.2010
 */
#ifndef LAMA_INPUTSET_HPP_
#define LAMA_INPUTSET_HPP_

#include <map>
#include <string>

#include <framework/src/config_framework.hpp>

/**
 * @brief The namespace sblas holds everything of the Template-Library sblas++
 */
namespace bf
{

class LAMABENCHFRAME_DLL_IMPORTEXPORT InputSet
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

#endif // LAMA_INPUTSET_HPP_
