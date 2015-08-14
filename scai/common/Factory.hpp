/**
 * @file Factory.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Template class for Factory 
 *
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#pragma once

#include <scai/common/config.hpp>
#include <scai/common/Exception.hpp>

#include <memory>
#include <map>
#include <vector>

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
class Factory
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
    class Register 
    {
    public:

        /** Constructor forces use of initialization. */

        Register();

        /** Method that does the registration. */

        inline static bool init();

        /** Boolean variable whose initialization does the registration. */

        static bool initialized;
    };

    /** @brief check if create is supported for a given input value */

    static bool canCreate( const InputType type );

    /** @brief Method to get all registered values. 
     *
     *  @param[out] values contains all registered values of InputType for 
     *              which a CreateFn has been registered.
     */

    static void getCreateValues( std::vector<InputType>& values );

protected:

    /** Create function that must be provided by derived classes */

    typedef OutputType ( *CreateFn ) ();

    /** Routine to register the create routine in the factory. */

    static void addCreator( const InputType type, CreateFn create );

private:

    typedef std::map<InputType, CreateFn> CreatorMap;

    /** Map container to get for the key the create function. */

    static CreatorMap& getFactory();
};

/* -----------------------------------------------------------------------------*/
/*  Implementation of methods for class Register of template class              */
/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType> 
template<class Derived>
Factory<InputType, OutputType>::Register<Derived>::Register()
{
    // Use of initialized forces its initialization and therefore registration

    if ( !initialized )
    {
        init();
    }
}

template<typename InputType, typename OutputType> 
template<class Derived>
inline bool Factory<InputType, OutputType>::Register<Derived>::init()
{
   Derived::addCreator( Derived::createValue(), Derived::create );
   return true;
}

template<typename InputType, typename OutputType> 
template<class Derived>
bool Factory<InputType, OutputType>::Register<Derived>::initialized = 
    Factory<InputType, OutputType>::Register<Derived>::init();

/* -----------------------------------------------------------------------------*/
/*  Implementation of methods for template class                                */
/* -----------------------------------------------------------------------------*/

template<typename InputType, typename OutputType> 
OutputType Factory<InputType, OutputType>::create( const InputType type )
{
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
    static std::auto_ptr<CreatorMap> factory;

    if ( !factory.get() )
    {
        factory = std::auto_ptr<CreatorMap>( new CreatorMap() );
    }

    return *factory;
}

template<typename InputType, typename OutputType>
void Factory<InputType, OutputType>::addCreator( const InputType type, CreateFn create )
{
    CreatorMap& factory = getFactory();

    // checks for multiple entries is not really necessary here, so just add entry in map container.

    factory.insert( std::pair<InputType, CreateFn>( type, create ) );
}

template<typename InputType, typename OutputType>
bool Factory<InputType, OutputType>::canCreate( InputType value )
{
    CreatorMap& factory = getFactory();

    typename CreatorMap::const_iterator it = factory.find( value );

    return it != factory.end();
}

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

}
