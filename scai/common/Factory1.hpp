/**
 * @file Factory1.hpp
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
 * @brief Template class for Factory where create routine has additional argument
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#pragma once

// local library
#include <scai/common/config.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

// std
#include <memory>
#include <map>
#include <vector>

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
template<typename InputType, typename ValueType, typename OutputType>
class COMMON_DLL_IMPORTEXPORT Factory1
{
public:

    /** This method creates an object by specifying two input values
     *
     *  @param[in] type is the first input value
     *  @param[in] val is the second input value
     *  @returns value of output type
     *
     *  @throws Exception if for type no create function has been registered.
     */
    static OutputType create( const InputType type, const ValueType val );

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

    /** @brief Method to get all registered values. */

    static void getCreateValues( std::vector<InputType>& values );

    static std::vector<InputType> getCreateValues();

protected:

    /** Create function that must be provided by derived classes */

    typedef OutputType ( *CreateFn ) ( ValueType val );

    /** Routine to register the create routine in the factory. */

    static void addCreator( const InputType type, CreateFn create );

    /** Routine to unregister the create routine in the factory. */

    static void removeCreator( const InputType type );

private:

    typedef std::map<InputType, CreateFn> CreatorMap;

    /** Map container to get for the key the create function. */

    static CreatorMap& getFactory();
};

template<typename InputType, typename ValueType, typename OutputType>
template<class Derived>
Factory1<InputType, ValueType, OutputType>::Register<Derived>::Register()
{
    // just some trick stuff to cheat most compilers so they instantiate the static variable
    if ( !registerGuard.initialized )
    {
        // Attention: caused problems on Intel MIC with dlopen
        // COMMON_THROWEXCEPTION( "Register without Guard" )
    }
}

/* -----------------------------------------------------------------------------*/
/*  Constructor/Destructor + static variable for registerGuard                  */
/* -----------------------------------------------------------------------------*/

// Constructor of the guard, does initialization

template<typename InputType, typename ValueType, typename OutputType>
template<class Derived>
Factory1<InputType, ValueType, OutputType>::Register<Derived>::RegisterGuard::RegisterGuard()
{
    Derived::addCreator( Derived::createValue(), Derived::create );
    initialized = true;
}

// Destructor of the guard, removes also entry from Factory

template<typename InputType, typename ValueType, typename OutputType>
template<class Derived>
Factory1<InputType, ValueType, OutputType>::Register<Derived>::RegisterGuard::~RegisterGuard()
{
    Derived::removeCreator( Derived::createValue() );
}

// ATTENTION: instantiation must be done explicitly for some compilers

template<typename InputType, typename ValueType, typename OutputType>
template<class Derived>
typename Factory1<InputType, ValueType, OutputType>::template Register<Derived>::RegisterGuard
Factory1<InputType, ValueType, OutputType>::Register<Derived>::registerGuard;

/* -----------------------------------------------------------------------------*/
/*  Implementation of methods for template class                                */
/* -----------------------------------------------------------------------------*/

template<typename InputType, typename ValueType, typename OutputType>
OutputType Factory1<InputType, ValueType, OutputType>::create( const InputType type, const ValueType val )
{
    OutputType value;
    const CreatorMap& factory = getFactory();
    typename CreatorMap::const_iterator fn = factory.find( type );

    if ( fn != factory.end() )
    {
        value = fn->second( val );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Factory: no creator for " << type << " available" )
    }

    return value;
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename ValueType, typename OutputType>
std::map<InputType, OutputType(* )( ValueType ) >& Factory1<InputType, ValueType, OutputType>::getFactory()
{
    // Factory might be already used during static initialization, so dynamic allocation is needed
    // Factory might be used at program exit, so it is never deleted
    static CreatorMap* factory = NULL;

    if ( factory == NULL )
    {
        factory =  new CreatorMap();
    }

    return *factory;
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename ValueType, typename OutputType>
void Factory1<InputType, ValueType, OutputType>::addCreator( const InputType type, CreateFn create )
{
    CreatorMap& factory = getFactory();
    // checks for multiple entries is not really necessary here, so just add entry in map container.
    factory.insert( std::pair<InputType, CreateFn>( type, create ) );
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename ValueType, typename OutputType>
void Factory1<InputType, ValueType, OutputType>::removeCreator( const InputType type )
{
    CreatorMap& factory = getFactory();
    // checks for multiple entries is not really necessary here, so just add entry in map container.
    factory.erase( type );
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename ValueType, typename OutputType>
bool Factory1<InputType, ValueType, OutputType>::canCreate( InputType value )
{
    CreatorMap& factory = getFactory();
    typename CreatorMap::const_iterator it = factory.find( value );
    return it != factory.end();
}

/* -----------------------------------------------------------------------------*/

template<typename InputType, typename ValueType, typename OutputType>
void Factory1<InputType, ValueType, OutputType>::getCreateValues( std::vector<InputType>& values )
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


template<typename InputType, typename ValueType, typename OutputType>
std::vector<InputType> Factory1<InputType, ValueType, OutputType>::getCreateValues()
{
    std::vector<InputType> values;
    getCreateValues( values );
    return values;
}

} /* end namespace common */

} /* end namespace scai */
