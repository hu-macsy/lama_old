/**
 * @file ConvertMatrixStorageBenchmark.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Benchmark to measure the storage conversion times between different sparse matrix storages.
 * @author Jiri Kraus and Bea Hornef
 * @date 02.12.2011
 * $Id$
 */

#pragma once

#include <scai/benchmark.hpp>

#include <scai/lama/benchmark/LAMAInputSetComplexityVisitor.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/lama/storage/CSRStorage.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;

template<typename StorageType1,typename StorageType2>
class ConvertMatrixStorageBenchmark: public bf::Benchmark
{
public:

    typedef StorageType1 SourceStorageType;
    typedef StorageType2 DestinationStorageType;
    typedef typename StorageType2::StorageValueType ValueType;

    static const std::string& id();

    ConvertMatrixStorageBenchmark();

    ConvertMatrixStorageBenchmark( const std::string& arguments );

    ConvertMatrixStorageBenchmark(
        const ConvertMatrixStorageBenchmark<SourceStorageType,DestinationStorageType>& other );

    virtual ~ConvertMatrixStorageBenchmark();

    /** Implementation of Benchmark::copy with covariant return type */

    virtual ConvertMatrixStorageBenchmark<SourceStorageType, DestinationStorageType>* copy() const;

    virtual short getValueTypeSize() const;

    virtual bool isThreadded() const;

    virtual const std::string& getId() const;

protected:
    virtual void initialize();
    virtual void setUp();
    virtual void execute();
    virtual void tearDown();
    virtual void shutdown();

    virtual CounterType getNumFloatingPointOperations() const;
    virtual CounterType getProcessedBytes() const;

private:
    SourceStorageType mSourceStorage;
    DestinationStorageType mDestinationStorage;
    int mDevice;

    HArray<IndexType> mCsrIA;
    HArray<IndexType> mCsrJA;
    HArray<ValueType> mCsrValues;

    static const LAMAInputSetComplexityVisitor::Group& group();

    static const std::string& sid();

    SCAI_LOG_DECL_STATIC_LOGGER(logger);
};

#ifndef SCAI_LOG_LEVEL_OFF

#define LAMA_KOMMA ,
template<typename StorageType1,typename StorageType2>
SCAI_LOG_DEF_LOGGER(ConvertMatrixStorageBenchmark<StorageType1 LAMA_KOMMA StorageType2>::logger, "Benchmark.ConvertMatrixStorageBenchmark");
#undef LAMA_KOMMA

#endif

template<typename StorageType1,typename StorageType2>
const std::string& ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename StorageType1,typename StorageType2>
const std::string& ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::id()
{

    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename StorageType1,typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::ConvertMatrixStorageBenchmark() : 

    bf::Benchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) )

{
    SCAI_LOG_ERROR( logger, "ConvertMatrixStorageBenchmark created" );
    // means HOST
    mDevice = -1;
}

template<typename StorageType1,typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::ConvertMatrixStorageBenchmark( const std::string& arguments ) : 

    bf::Benchmark( sid(), 
    LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments )

{
    std::istringstream arg( arguments );
    arg >> mDevice;
    SCAI_LOG_ERROR( logger, "ConvertMatrixStorageBenchmark created for CUDA (device " << mDevice << ")" );
}

template<typename StorageType1,typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::ConvertMatrixStorageBenchmark(
    const ConvertMatrixStorageBenchmark<SourceStorageType,DestinationStorageType>& other )
    : bf::Benchmark( other ), mDevice( other.mDevice )
{
    SCAI_LOG_ERROR( logger, "copy constructror" )
}

template<typename StorageType1,typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::~ConvertMatrixStorageBenchmark()
{
    SCAI_LOG_ERROR( logger, "~ConvertMatrixStorageBenchmark" )
}

template<typename StorageType1,typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>* ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::copy() const
{
    SCAI_LOG_ERROR( logger, "copy ConvertMatrixStorageBenchmark" )

    return  new ConvertMatrixStorageBenchmark<SourceStorageType,DestinationStorageType>( *this );
}

template<typename StorageType1,typename StorageType2>
short ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getValueTypeSize() const
{
    return sizeof( ValueType );
}

template<typename StorageType1,typename StorageType2>
bool ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::isThreadded() const
{
    return true;
}

template<typename StorageType1,typename StorageType2>
const std::string& ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getId() const
{
    return id();
}

template<typename StorageType1,typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::initialize()
{
    SCAI_LOG_ERROR( logger, "initialize, mInputSetId = " << mInputSetId )

    const LAMAInputSet& inputSet = bf::InputSetRegistry<LAMAInputSet>::getRegistry().get( mInputSetId );

    SCAI_LOG_ERROR( logger, "inputSet" )

    //TODO: assign should do the value type conversion

    const CSRStorage<double>& csrStorage = inputSet.getA().getLocalStorage();

    SCAI_LOG_ERROR( logger, "csrStorage = " << csrStorage)

    mSourceStorage.setCSRData( csrStorage.getNumRows(), csrStorage.getNumColumns(), csrStorage.getNumValues(),
                               csrStorage.getIA(), csrStorage.getJA(), csrStorage.getValues() );

    ContextPtr context;

    if( mDevice < 0 )
    {
        context = Context::getContextPtr( Context::Host );
    }
    else
    {
        context = Context::getContextPtr( Context::CUDA, mDevice );
    }

    mSourceStorage.setContextPtr( context );
    mDestinationStorage.setContextPtr( context );
}

template<typename StorageType1,typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::setUp()
{
    mSourceStorage.prefetch( mSourceStorage.getContextPtr() );
    mDestinationStorage.prefetch( mDestinationStorage.getContextPtr() );
    mSourceStorage.wait();
    mDestinationStorage.wait();

//    mDestinationStorage._MatrixStorage::assign( mSourceStorage );
    mSourceStorage.buildCSRData( mCsrIA, mCsrJA, mCsrValues );

}

template<typename StorageType1,typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::execute()
{
    SCAI_LOG_ERROR( logger, "execute conversion" )

    mDestinationStorage.setContextPtr( mSourceStorage.getContextPtr() );

    // Note: conversion will be done at context of target matrix mDestinationStorage

    mDestinationStorage.setCSRData( mSourceStorage.getNumRows(), mSourceStorage.getNumColumns(), mCsrJA.size(), mCsrIA,
                                    mCsrJA, mCsrValues );

    SCAI_LOG_ERROR( logger, "execute conversion done" )
}

template<typename StorageType1,typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::tearDown()
{
}

template<typename StorageType1,typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::shutdown()
{
}

template<typename StorageType1,typename StorageType2>
CounterType ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getNumFloatingPointOperations() const
{
    return 0;
}

template<typename StorageType1,typename StorageType2>
CounterType ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getProcessedBytes() const
{
    return 0;
}

