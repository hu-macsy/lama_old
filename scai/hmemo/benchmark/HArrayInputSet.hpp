/**
 * @file HArrayInputSet.hpp
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
 * @brief Definition of input set class to be used for HArray benchmarks.
 * @author Thomas Brandes
 * @date 17.09.2017
 */

#pragma once

#include <scai/benchmark.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/Context.hpp>

namespace scai
{

namespace hmemo
{

/** HArray input set, specified by context + size, type is always "double" 
 *
 *  - The first argument specifies the context on which the data is allocated
 *  - The second argument specifies the size of the HArray
 *
 *  \code
 *     HArrayInputSet inputData( "Host, 100000 );
 *     std::unique_ptr<HArrayInputSet> inputData1( InputSet::create( "HArrayInputSet", "CUDA, 2000" ) );
 *     std::unique_ptr<HArrayInputSet> inputData2( InputSet::createWithArg( "HArrayInputSet( CUDA, 20000 )" ) );
 *  \endcode
 */
class COMMON_DLL_IMPORTEXPORT HArrayInputSet: 

    public benchmark::InputSet,
    private benchmark::InputSet::Register<HArrayInputSet>

{
public:

    /** Constructor of a random LAMA input set, specified by size, fillRate */

    HArrayInputSet( const std::string& argument );

    virtual ~HArrayInputSet();

    /** Override implementation InputSet::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /** Implement method createValue that returns the key used for registration in factory. */

    static std::string createValue()
    {
        return "HArrayInputSet";
    }

    /** Static method required for create, called by InputSet::create( HArrayInputSet::createValue(), argument ) */

    static InputSet* create( const std::string argument )
    {
        return new HArrayInputSet( argument );
    }

    /** Implementation of pure method InputSet::getCreateId()  */

    virtual const std::string& getCreateId() const
    {
        static std::string id = createValue();
        return id;
    }

    /** Implementation of pure method InputSet::getArgument()  */

    virtual const std::string& getArgument() const
    {
        return mArgument;
    }

    /** Return a reference to the input set array. */

    HArray<double>& getArray();

private:

    HArrayInputSet();

    SCAI_LOG_DECL_STATIC_LOGGER( logger );

    std::string mArgument;

    ContextPtr mContext;

    IndexType mSize;

    std::unique_ptr<HArray<double> > mArray;
};

}

}
