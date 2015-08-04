/**
 * @file Token.cpp
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
 * @brief Demo program of using Task and SyncToken
 * @author Thomas Brandes
 * @date 02.02.2012
 * @since 1.0.0
 */

#include <tasking/TaskSyncToken.hpp>

#include <common/Exception.hpp>
#include <logging.hpp>

#include <common/bind.hpp>
#include <common/shared_ptr.hpp>

LAMA_LOG_DEF_LOGGER( logger, "Token" )

/* --------------------------------------------------------------------- */

void task( int a[], const int b[], const int c[], int N )
{
    LAMA_LOG_INFO( logger, "Run task with N = " << N )
    LAMA_LOG_INFO( logger, "a = " << a << ", b = " << b << ", c = " << c )

    for ( int i = 0; i < N; ++i )
    {
        a[i] = b[i] + c[i];
    }

    LAMA_LOG_INFO( logger, "task done" )

    // Note: COMMON_THROWEXCEPTION( "I hate you" ) would be handled by Task
}

/* --------------------------------------------------------------------- */

using namespace tasking;

using common::shared_ptr;

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

    shared_ptr<SyncToken> t ( new TaskSyncToken( common::bind( task, a, b, c , N ) ) );

    t->wait();

    for ( int i = 0; i < N; ++i )
    {
        COMMON_ASSERT_EQUAL( a[i], i + 3, "wrong value" )
    }

    delete [] c;
    delete [] b;
    delete [] a;

}

struct Data : private common::NonCopyable, public SyncTokenMember
{
    Data( int N ) 
    {
        LAMA_LOG_INFO( logger, "construct Data( N = " << N << " )" )

        mA = new int[N];
        mB = new int[N];
        mC = new int[N];

        LAMA_LOG_INFO( logger, "a = " << mA << ", b = " << mB << ", c = " << mC )

        for ( int i = 0; i < N; ++i )
        {
            mB[i] = i;
            mC[i] = 3;
        }
    }

    ~Data()
    {
        LAMA_LOG_INFO( logger, "~Data, a = " << mA << ", b = " << mB << ", c = " << mC )
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
 
    shared_ptr<SyncToken> t ( new TaskSyncToken( common::bind( task, data->mA, data->mB, data->mC , N ) ) );

    // give ownership of data to the sync token 

    t->pushToken( data );

    return t;
}

int main()
{
    // simple();

    shared_ptr<SyncToken> t = run( 100000 );

    LAMA_LOG_INFO( logger, "run started" )

    t->wait();

    LAMA_LOG_INFO( logger, "wait done" )

    // destructor of token t waits for synchronization and frees data

}
