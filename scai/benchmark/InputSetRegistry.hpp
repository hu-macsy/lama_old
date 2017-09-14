/**
 * @file InputSetRegistry.hpp
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
 * @brief Template class for registration of input sets.
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

#include <string>
#include <memory>
#include <map>
#include <sstream>
#include <cstdio>

#include <scai/benchmark/macros.hpp>
#include <scai/benchmark/BFException.hpp>
#include <scai/benchmark/InputSet.hpp>
#include <scai/benchmark/BaseInputSetRegistry.hpp>
#include <scai/benchmark/string_helper.hpp>

#include <scai/common/exception/Exception.hpp>

namespace scai
{

namespace bf
{

template<typename InputSetT>
class InputSetCreator
{
public:
    /**
     * @brief Default destructor. Frees all inner ressources, if necessary.
     */
    virtual ~InputSetCreator();

    typedef InputSetT InputSetType;
    /**
     * @brief Creates an object of InputSet and returns it within an auto_ptr.
     * @return An object of InputSet within an auto_ptr.
     */
    virtual InputSetType* create() const = 0;

    /**
     * @brief Creates an object of InputSet from the given params and returns
     *        it within an auto_ptr.
     * @param[in] params    parameters controlling the creation of the input set
     * @return              An object of InputSet within an auto_ptr.
     */
    virtual InputSetType* create( const std::string& params ) const;
};

template<typename InputSetT>
InputSetCreator<InputSetT>::~InputSetCreator()
{
}

template<typename InputSetT>
InputSetT* InputSetCreator<InputSetT>::create( const std::string& ) const
{
    throw BFException( "InputSetCreator::create( args ): This InputSet does not allow arguments." );
}

template<typename InputSetT>
class InputSetRegistry: public BaseInputSetRegistry
{
public:
    typedef InputSetT InputSetType;
    typedef typename std::map<std::string,InputSetCreator<InputSetType>*>::const_iterator const_iterator;

    /**
     * @brief Destructor. Frees all inner ressources, if nessecary.
     */
    virtual ~InputSetRegistry();

    /**
     * @brief Returns the InputSetRegistry. If not yet returned, the
     *        InputSetRegistry is created.
     * @return The InputSetRegistry.
     */
    static InputSetRegistry& getRegistry();

    /**
     * @brief Releases all resource of the registry. Once this method is called it is no longer possible to
     *        restore the orignial state.
     */
    static void freeRegistry();

    /**
     * @brief Returns the InputSetType that fits the id.
     *
     * Returns the InputSetType that fits the id. If the InputSetType has not
     * been returned, yet, it will be created.
     *
     * @param[in] id The id of the InputSetType that will be returned. The id
     *               may holds arguments for the InputSet. So id is either
     *               id=name(args) or id=name, where 'name' is the name of the
     *               InputSet and args are the arguments for the InputSet. This
     *               also means, that a name must not contain an open or closed
     *               braked ['(' or ')'].
     * @return The InputSetType fitting to the id.
     * @throws Exception, if the InputSetType could not be created.
     */
    InputSetType& get( const std::string& id ) const;

    /**
     * @brief Adds the InputSetCreator with the given id to the registry.
     * @param[in] id                The id of the InputSetCreator.
     * @param[in] inputSetCreator   The InputSetCreator.
     */
    void add( const std::string& id, InputSetCreator<InputSetType>* inputSetCreator );

    bool has( const std::string& id ) const;

    /**
     * @brief Returns iterator to beginning.
     * @return Iterator to beginning.
     */
    const_iterator begin() const;
    /**
     * @brief Returns iterator to end.
     * @return Iterator to end.
     */
    const_iterator end() const;

    void clear();

    /**
     * @brief Returns the Id used for the Special InputSetCreator that takes a
     *        filename.
     *
     * @return The Special InputSetCreator Id.
     */
    static const std::string& getFileInputSetCreatorId();

    /**
     * @brief Returns the Id used for the Special InputSetCreator that creates a
     *        random input set.
     *
     * @return The Special InputSetCreator Id.
     */
    static const std::string& getRandomInputSetCreatorId();

    virtual void getInputSetMap( std::map<std::string,std::string>& isetIds ) const;

private:
    InputSetRegistry();
    InputSetRegistry( const InputSetRegistry& other );
    static InputSetRegistry* m_instance;
    std::map<std::string,InputSetCreator<InputSetType>*> m_creatorMap;
    mutable std::map<std::string,InputSetType*> m_inputSetMap;

    static bool isFreed;

    class CGuard
    {
    public:
        /**
         * @brief Constructor.
         */
        CGuard();
        /**
         * @brief Deletes all InputSetRegistries.
         */
        ~CGuard();
    };
    friend class CGuard;
};

template<typename InputSetT>
bool InputSetRegistry<InputSetT>::isFreed = false;

template<typename InputSetT>
const std::string& InputSetRegistry<InputSetT>::getFileInputSetCreatorId()
{
    static const std::string fileInputSetCreatorId = "File";
    return fileInputSetCreatorId;
}

template<typename InputSetT>
const std::string& InputSetRegistry<InputSetT>::getRandomInputSetCreatorId()
{
    static const std::string randomInputSetCreatorId = "Random";
    return randomInputSetCreatorId;
}

/** Array of InputSetRegistries. */
template<typename InputSetT>
InputSetRegistry<InputSetT>* InputSetRegistry<InputSetT>::m_instance = 0;

template<typename InputSetT>
InputSetRegistry<InputSetT>::InputSetRegistry()
{
}

template<typename InputSetT>
InputSetRegistry<InputSetT>::~InputSetRegistry()
{
    while( !m_inputSetMap.empty() )
    {
        typename std::map<std::string,InputSetType*>::iterator begin = m_inputSetMap.begin();
        InputSetType* ptr = begin->second;
        m_inputSetMap.erase( begin );
        delete ptr;
    }

    while( !m_creatorMap.empty() )
    {
        typename std::map<std::string,InputSetCreator<InputSetType>*>::iterator begin = m_creatorMap.begin();
        InputSetCreator<InputSetType>* ptr = begin->second;
        m_creatorMap.erase( begin );
        delete ptr;
    }
}

template<typename InputSetT>
InputSetRegistry<InputSetT>& InputSetRegistry<InputSetT>::getRegistry()
{
    static CGuard g;

    if( isFreed )
    {
        COMMON_THROWEXCEPTION( "InputSetRegistry is already freed. You can't get an instance." )
    }

    if( m_instance == 0 )
    {
        m_instance = new InputSetRegistry<InputSetType>();
    }

    return *m_instance;
}

template<typename InputSetT>
void InputSetRegistry<InputSetT>::freeRegistry()
{
    if( m_instance != 0 )
    {
        delete m_instance;
        m_instance = 0;
        isFreed = true;
    }
}

template<typename InputSetT>
typename InputSetRegistry<InputSetT>::InputSetType&
InputSetRegistry<InputSetT>::get( const std::string& id ) const
{
    typename std::map<std::string,InputSetType*>::const_iterator loc = m_inputSetMap.find( id );

    if( loc == m_inputSetMap.end() )
    {
        std::string name = "";
        std::string args = "";
        std::string::size_type openBraked = id.find( '(', 0 );

        if( openBraked != std::string::npos )
        {
            ++openBraked;
            std::string::size_type closedBraked = id.find( ')', openBraked );

            if( closedBraked != std::string::npos )
            {
                name = id.substr( 0, openBraked - 1 );
                args = id.substr( openBraked, closedBraked - openBraked );
            }
            else
            {
                std::stringstream message;
                message << "InputSetRegistry<InputSetT>::get: missing ')' in " << "InputSet name " << id;
                throw BFException( message.str() );
            }

            trimm( name );
            trimm( args );

            if( name.empty() )
            {
                std::stringstream message;
                message << "InputSetRegistry<InputSetT>::get: missing name in " << "InputSet name " << id;
                throw BFException( message.str() );
            }

            if( args.empty() )
            {
                std::stringstream message;
                message << "InputSetRegistry<InputSetT>::get: missing arguments " << "in InputSet name " << id;
                throw BFException( message.str() );
            }
        }
        else if( id.find( ')', 0 ) != std::string::npos )
        {
            std::stringstream message;
            message << "InputSetRegistry<InputSetT>::get: missing '(' in InputSet " << "name " << id;
            throw BFException( message.str() );
        }
        else
        {
            name = id;
        }

        typename std::map<std::string,InputSetCreator<InputSetType>*>::const_iterator cloc = m_creatorMap.find( name );

        if( cloc == m_creatorMap.end() )
        {
            throw BFException(
                "InputSetRegistry<InputSetT>::get: No creator Function for InputSet " + name
                + " available" );
        }

        if( args.empty() )
        {
            m_inputSetMap.insert( std::make_pair( id, cloc->second->create() ) );
        }
        else
        {
            m_inputSetMap.insert( std::make_pair( id, cloc->second->create( args ) ) );
        }
    }

    return *m_inputSetMap[id];
}

template<typename InputSetT>
void InputSetRegistry<InputSetT>::add( const std::string& id, InputSetCreator<InputSetType>* inputSetCreator )
{
    m_creatorMap.insert( std::make_pair( id, inputSetCreator ) );
}

template<typename InputSetT>
typename InputSetRegistry<InputSetT>::const_iterator InputSetRegistry<InputSetT>::begin() const
{
    return m_creatorMap.begin();
}

template<typename InputSetT>
typename InputSetRegistry<InputSetT>::const_iterator InputSetRegistry<InputSetT>::end() const
{
    return m_creatorMap.end();
}

template<typename InputSetT>
void InputSetRegistry<InputSetT>::clear()
{
    m_inputSetMap.clear();
}

template<typename InputSetT>
void InputSetRegistry<InputSetT>::getInputSetMap( std::map<std::string,std::string>& isetIds ) const
{
    for( const_iterator it = m_creatorMap.begin(); it != m_creatorMap.end(); ++it )
    {
        isetIds.insert( std::make_pair( it->first, "name" ) );
    }
}

template<typename InputSetT>
InputSetRegistry<InputSetT>::CGuard::CGuard()
{
}

template<typename InputSetT>
InputSetRegistry<InputSetT>::CGuard::~CGuard()
{
    if( InputSetRegistry<InputSetT>::m_instance != 0 )
    {
        delete InputSetRegistry<InputSetT>::m_instance;
        InputSetRegistry<InputSetT>::m_instance = 0;
    }
}

template<typename T>
class InputSetRegistration
{
public:
    /**
     * @brief Registers an InputSetRegistration-object with the given id.
     * @param[in] id The id of the InputSetRegistration to be registered.
     */
    InputSetRegistration( const std::string& id )
    {
        InputSetRegistry<typename T::InputSetType>& reg = InputSetRegistry<typename T::InputSetType>::getRegistry();
        reg.add( id, new T() );
    }
};

/**
 * @brief Registrates the given type.
 */
#define LAMA_INPUTSET_REGISTRATION(creatorType)                                     \
    static bf::InputSetRegistration<creatorType>                            \
    LAMA_UNIQUE_NAME( iSetRegObj, creatorType )(creatorType::id())

} // namespace bf

} // namespace scai
