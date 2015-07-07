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
 * @author Thomas Brandes, Jiri Krause
 * @date 14.03.2011
 * @revised 03.07.2015
 */

#pragma once

#include <common/config.hpp>
#include <common/Exception.hpp>

#include <memory>
#include <map>

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

    template<class classname, InputType value>
    class Register 
    {
    public:

        Register()
        {
            if ( !initialized )
            {
                init();
            }
        }

        inline static bool init()
        {
            classname::addCreator( value, classname::create );
            return true;
        }

       static bool initialized;
    };

protected:

    typedef OutputType ( *CreateFn ) ();

    static void addCreator( const InputType type, CreateFn create );

private:

    typedef std::map<InputType, CreateFn> CreatorMap;

    /** Map container to get for the key the create function. */

    static CreatorMap& getFactory();
};

/*
template<typename InputType, typename OutputType> 
template<class classname, int value>
inline bool Factory<InputType, OutputType>::Register<classname,value>::init()
{
   classname::addCreator( value, classname::create )
   return true;
}
*/

template<typename InputType, typename OutputType> 
template<class classname, InputType value>
bool Factory<InputType, OutputType>::Register<classname,value>::initialized = 
    Factory<InputType, OutputType>::Register<classname,value>::init();

template<typename InputType, typename OutputType> 
OutputType Factory<InputType, OutputType>::create( const InputType type )
{
    OutputType value = 0;

    const CreatorMap& factory = getFactory();

    typename CreatorMap::const_iterator fn = factory.find( type );

    if ( fn  != factory.end() )
    {
        value = fn->second();
    }
    else
    {
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
    factory[type] = create;
}

}
