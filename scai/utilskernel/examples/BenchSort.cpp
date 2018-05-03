/**
 * @file BenchSort.cpp
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
 * @brief Benchmarking of sort
 * @author Thomas Brandes
 * @date 21.10.2016
 */

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

using namespace std;
using namespace scai;
using namespace hmemo;
using namespace utilskernel;

void bucketSort( const IndexType N )
{
    IndexType nBuckets = 2048;

    HArray<IndexType> values( 2400 );

    HArrayUtils::setRandom( values, nBuckets - 1 );

    // use setScalar on IndexType, binaryOp is only for numeric types
    // HArrayUtils::setScalar( values, nBuckets, common::BinaryOp::MODULO );

    HArray<IndexType> offsets;
    HArray<IndexType> perm;

    double start = common::Walltime::get();

    HArrayUtils::bucketSort( offsets, perm, values, nBuckets );

    HArray<IndexType> sortedValues;
    HArrayUtils::gather( sortedValues, values, perm, common::BinaryOp::COPY );

    double time1 = common::Walltime::get() - start;

    cout << "Bucket sort of " << N << " values took " << time1 << " seconds." << endl;

    bool isSorted = HArrayUtils::isSorted( sortedValues, common::CompareOp::LE );

    if ( isSorted )
    {
        cout << "BucketSort: Values are well sorted" << endl;
    }
    else
    {
        cout << "Attention: values are not sorted by bucket sort" << endl;
    }
}

template<typename ValueType>
void sort( const IndexType N )
{
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    HArray<ValueType> values( ctx );
    HArray<IndexType> perm;

    values.resize( N );  
    HArrayUtils::setRandom( values, 1 );

    bool ascending = true;

    double start = common::Walltime::get();

    HArrayUtils::sort( values, perm, ascending, ctx );

    double time1 = common::Walltime::get() - start;

    cout << "Sort of " << N << " values took " << time1 << " seconds." << endl;

    bool isSorted = HArrayUtils::isSorted( values, ascending );

    if ( isSorted )
    {
        cout << "Values are well sorted" << endl;
    }
    else
    {
        cout << "Attention: values are not sorted" << endl;
    }
}

int main( int argc, const char* argv[] )
{
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType N = 10 * 1000 * 1000;

    // sort<double>( N );
    // sort<IndexType>( N );
    bucketSort( N );
}
