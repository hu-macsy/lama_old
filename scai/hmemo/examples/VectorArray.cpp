/**
 * @file hmemo/examples/VectorArray.cpp
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
 * @brief Use of heterogeneous arrays in C++ container classes
 * @author Thomas Brandes, Andreas Borgen Longva
 * @date 07.12.2017
 */

#include <vector>
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>

#include <iostream>
#include <algorithm>

using namespace scai;
using namespace hmemo;

struct less_than_key
{
    inline bool operator() ( const _HArray& array1, const _HArray& array2 )
    {
        return array1.size() < array2.size();
    }
};

int main()
{
    typedef DefaultReal ValueType;

    const IndexType N = 20;   // number of arrays
   
    const IndexType N_VALUES = 10000;

    const ValueType initVal = 20;

    const ValueType* ptrArray[N];

    std::vector<HArray<ValueType> > myData;

    for ( IndexType i = 0; i < N; ++i ) 
    {
        myData.push_back( HArray<ValueType>( IndexType( N_VALUES - i ), initVal + i ) );
    }

    for ( IndexType i = 0; i < N; ++i ) 
    {
        std::cout << "array " << i << " = " << myData[i] << std::endl;
        ReadAccess<ValueType> rA( myData[i] );
        ptrArray[i] = rA.get();
    }

    std::cout << "Sort the arrays by size" << std::endl;

    std::sort( myData.begin(), myData.end(), less_than_key() );

    std::cout << "End sorting" << std::endl;

    for ( IndexType i = 0; i < N; ++i ) 
    {
        std::cout << "array " << i << " = " << myData[i] << std::endl;

        ReadAccess<ValueType> rA( myData[i] );

        if ( rA.get() != ptrArray[N - 1 - i] )
        {
            std::cout << "SERIOUS problem: there was reallocation during sorting." << std::endl;
        }
    }
}
