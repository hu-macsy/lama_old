/**
 * @file Token.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Demo program of using Task and SyncToken
 * @author Thomas Brandes
 * @date 02.02.2012
 */

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/logging.hpp>

#include <memory>
#include <functional>

using std::shared_ptr;

SCAI_LOG_DEF_LOGGER( logger, "Token" )

/* --------------------------------------------------------------------- */

void task( int a[], const int b[], const int c[], int N )
{
    SCAI_LOG_INFO( logger, "Run task with N = " << N )
    SCAI_LOG_INFO( logger, "a = " << a << ", b = " << b << ", c = " << c )

    for ( int i = 0; i < N; ++i )
    {
        a[i] = b[i] + c[i];
    }

    SCAI_LOG_INFO( logger, "task done" )
    // Note: COMMON_THROWEXCEPTION( "I hate you" ) would be handled by Task
}

/* --------------------------------------------------------------------- */

using namespace scai::tasking;

void simple()
{
    const int N = 10;
    int* a = new int[N];
    int* b = new int[N];
    int* c = new int[N];

    for ( int i = 0; i < N; ++i )
    {
        b[i] = i;
        c[i] = 3;
    }

    // using shared pointer will delete token automatically
    shared_ptr<SyncToken> t ( new TaskSyncToken( std::bind( task, a, b, c , N ) ) );
    t->wait();

    for ( int i = 0; i < N; ++i )
    {
        SCAI_ASSERT_EQUAL( a[i], i + 3, "wrong value" )
    }

    delete [] c;
    delete [] b;
    delete [] a;
}

struct Data : private scai::common::NonCopyable, public SyncTokenMember
{
    Data( int N )
    {
        SCAI_LOG_INFO( logger, "construct Data( N = " << N << " )" )
        mA = new int[N];
        mB = new int[N];
        mC = new int[N];
        SCAI_LOG_INFO( logger, "a = " << mA << ", b = " << mB << ", c = " << mC )

        for ( int i = 0; i < N; ++i )
        {
            mB[i] = i;
            mC[i] = 3;
        }
    }

    ~Data()
    {
        SCAI_LOG_INFO( logger, "~Data, a = " << mA << ", b = " << mB << ", c = " << mC )
        delete [] mC;
        delete [] mB;
        delete [] mA;
    }

    int* mA;
    int* mB;
    int* mC;
};

shared_ptr<SyncToken> run( int N )
{
    shared_ptr<Data> data( new Data( N ) );
    // call task synchronously
    task( data->mA, data->mB, data->mC, N );
    // call task asynchronously
    shared_ptr<SyncToken> t ( new TaskSyncToken( std::bind( task, data->mA, data->mB, data->mC , N ) ) );
    // give ownership of data to the sync token
    t->pushToken( data );
    return t;
}

int main()
{
    // simple();
    shared_ptr<SyncToken> t = run( 100000 );
    SCAI_LOG_INFO( logger, "run started" )
    t->wait();
    SCAI_LOG_INFO( logger, "wait done" )
    // destructor of token t waits for synchronization and frees data
}
