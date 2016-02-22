/**
 * @file common/examples/BenchPointers.cpp
 *
 * @brief Compare overhead of unique and shared pointers
 *
 * @author Thomas Brandes
 * @date 11.02.2016
 */

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/shared_ptr.hpp>
#include <scai/common/macros/assert.hpp>

class Data
{
public:

    Data()
    {
        m1 = 2;
        m2 = 1;
    }

    int m1, m2;
};

#include <iostream>

using scai::common::Walltime;
using scai::common::shared_ptr;
using scai::common::unique_ptr;

void usualPointer( int& dummy )
{
    Data* d1 = new Data;
    Data* d2 = new Data;

    dummy = dummy + d1->m1 - d2->m2;

    delete d1;
    delete d2;
}

void sharedPointer( int& dummy )
{
    shared_ptr<Data> d1 ( new Data );
    shared_ptr<Data> d2 ( new Data );

    dummy = dummy + d1->m1 - d2->m2;
}

void uniquePointer( int& dummy )
{
    unique_ptr<Data> d1 ( new Data );
    unique_ptr<Data> d2 ( new Data );

    dummy = dummy + d1->m1 - d2->m2;
}

void sub1( shared_ptr<Data> d )
{
    d->m2++;
}

void sub2( shared_ptr<Data>& d )
{
    d->m2++;
}

int main()
{
    const int NITER = 1000 * 1000;

    int dummy = 0;

    double t = Walltime::get();

    for ( int i = 0; i < NITER; ++i )
    {
         usualPointer( dummy );
    }

    double t1 = Walltime::get() - t;

    t = Walltime::get();

    for ( int i = 0; i < NITER; ++i )
    {
         sharedPointer( dummy);
    }

    double t2 = Walltime::get() - t;

    t = Walltime::get();

    for ( int i = 0; i < NITER; ++i )
    {
         uniquePointer( dummy);
    }

    double t3 = Walltime::get() - t;

    std::cout << "Time allocate/free usual  pointer = " << t1 << std::endl;
    std::cout << "Time allocate/free shared pointer = " << t2 << std::endl;
    std::cout << "Time allocate/free unique pointer = " << t3 << std::endl;

    // Measure overhead of shared pointer copy

    shared_ptr<Data> data1( new Data() );
    shared_ptr<Data> data2( new Data() );

    t = Walltime::get();

    for ( int i = 0; i < NITER; ++i )
    {
         sub1( data1 );
    }

    double t4 = Walltime::get() - t;

    t = Walltime::get();

    for ( int i = 0; i < NITER; ++i )
    {
         sub2( data2 );
    }

    double t5 = Walltime::get() - t;

    SCAI_ASSERT_EQUAL( data1->m2, data2->m2, "routines should have delivered same results" )

    std::cout << "Timing shared pointer arguments" << std::endl;
    std::cout << "===============================" << std::endl;

    std::cout << "Time shared pointer (copy) = " << t4 << std::endl;
    std::cout << "Time shared pointer (ref)  = " << t5 << std::endl;

}
