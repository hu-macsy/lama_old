/**
 * @file Factory.hpp
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
 * @brief Template class for Factory
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#pragma once

// local library
#include <scai/common/config.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

// std
#include <typeinfo>
#include <memory>
#include <map>
#include <vector>
#include <iostream>
#include <utility>

template<typename T1, typename T2>
inline std::ostream& operator<<( std::ostream& stream, const std::pair<T1, T2>& object )
{
    stream << object.first << object.second;
    return stream;
}

namespace scai
{

namespace common
{


/** @brief Templatate class for a Factory where objects of OutputType
 *         are created by a specific creator routine for different values
 *         of "InputType".
 *
 *  A factory is especially useful for dynamic extensions of libraries
 *  as creator functions might be registered later.
 */
template<typename InputType, typename OutputType>
class COMMON_DLL_IMPORTEXPORT Factory
{
public:

    /** This method creates an object by specifying an input value.
     *
     *  @param[in] type is the input value
     *  @returns value of output type
     *
     *  @throws Exception if for type no create function has been registered.
     */
    static OutputType create( const InputType type );

    /** This template class can be used as base class for derived classes
     *  to force registration.
     *
     *  Derived class must provide create function and createValue
     */
    template<class Derived>
    class COMMON_DLL_IMPORTEXPORT Register
    {
    public:

        /** Constructor forces use of initialization. */

        Register();

        /** Guard class that registers/unregisters the creator.   */

        class COMMON_DLL_IMPORTEXPORT RegisterGuard
        {
        public:
            RegisterGuard();
            ~RegisterGuard();
            bool initialized;
        };

        /** Static guard variable takes care of registration */

        static RegisterGuard registerGuard;
    };

    /** @brief check if create is supported for a given input value */

    static bool canCreate( const InputType type );

    /** @brief Method to get all registered values.
     *
     *  @param[out] values contains all registered values of InputType for
     *              which a CreateFn has been registered.
     */

    static void getCreateValues( std::vector<InputType>& values );

    static std::vector<InputType> getCreateValues();

protected:

    /** Create function that must be provided by derived classes */

    typedef OutputType ( *CreateFn ) ();

    /** Routine to register the create routine in the factory. */

    static void addCreator( const InputType type, CreateFn create );

    static void removeCreator( const InputType type );

private:

    typedef std::map<InputType, CreateFn> CreatorMap;

    /** Map container to get for the key the create function.
     *  (if the pointer is NULL, destructor has been called)
     */

    static CreatorMap& getFactory();
};

/* -----------------------------------------------------------------------------*/
/*  Implementation of methods for class Register of template class              */
/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType>
template<class Derived>
Factory<InputType, OutputType>::Register<Derived>::Register()
{
    // just some trick stuff to cheat most compilers so they instantiate the static variable
    if ( !registerGuard.initialized )
    {
        COMMON_THROWEXCEPTION( "Register without Guard" )
    }
}

/* -----------------------------------------------------------------------------*/
/*  Constructor/Destructor  RegisterGuard                                       */
/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType>
template<class Derived>
Factory<InputType, OutputType>::Register<Derived>::RegisterGuard::RegisterGuard()
{
    Derived::addCreator( Derived::createValue(), Derived::create );
    initialized = true;
}

template<typename InputType, typename OutputType>
template<class Derived>
Factory<InputType, OutputType>::Register<Derived>::RegisterGuard::~RegisterGuard()
{
    Derived::removeCreator( Derived::createValue() );
}

// ATTENTION: this instantiation is not sufficient for some compilers,
//            so it must be done explicitly later for some compilers

template<typename InputType, typename OutputType>
template<class Derived>
typename Factory<InputType, OutputType>::template Register<Derived>::RegisterGuard
Factory<InputType, OutputType>::Register<Derived>::registerGuard;

/* -----------------------------------------------------------------------------*/
/*  Implementation of methods for template class                                */
/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType>
OutputType Factory<InputType, OutputType>::create( const InputType type )
{
    using ::operator<<;
    OutputType value;
    const CreatorMap& factory = getFactory();
    typename CreatorMap::const_iterator fn = factory.find( type );

    if ( fn  != factory.end() )
    {
        value = fn->second();
    }
    else
    {
        // Be careful: operator<< for InputType must be available
        COMMON_THROWEXCEPTION( "Factory: no creator for " << type << " available" )
    }

    return value;
}

template<typename InputType, typename OutputType>
std::map<InputType, OutputType(* )() >& Factory<InputType, OutputType>::getFactory()
{
    // Factory will be created during static initialization when it is needed for the first time
    // A destructor is never called as registered objects still might remove themselves later
    static CreatorMap* factory = NULL;

    if ( factory == NULL )
    {
        factory = new CreatorMap();
    }

    return *factory;
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType>
void Factory<InputType, OutputType>::addCreator( const InputType type, CreateFn create )
{
    CreatorMap& factory = getFactory();
    // checks for multiple entries is not really necessary here, so just add entry in map container.
    factory.insert( std::pair<InputType, CreateFn>( type, create ) );
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType>
void Factory<InputType, OutputType>::removeCreator( const InputType type )
{
    // Note: factory is never deleted so this call is safe at program exit
    CreatorMap& factory = getFactory();
    // the following call is also safe if factory does not contain the creator for type
    factory.erase( type );
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType>
bool Factory<InputType, OutputType>::canCreate( InputType value )
{
    CreatorMap& factory = getFactory();
    typename CreatorMap::const_iterator it = factory.find( value );
    return it != factory.end();
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType>
void Factory<InputType, OutputType>::getCreateValues( std::vector<InputType>& values )
{
    CreatorMap& factory = getFactory();
    values.clear();
    values.reserve( factory.size() );
    typename CreatorMap::const_iterator it;

    for ( it = factory.begin(); it != factory.end(); ++it )
    {
        values.push_back( it->first );
    }
}

template<typename InputType, typename OutputType>
std::vector<InputType> Factory<InputType, OutputType>::getCreateValues()
{
    std::vector<InputType> values;
    getCreateValues( values );
    return values;
}

} /* end namespace common */

} /* end namespace scai */
