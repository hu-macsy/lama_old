/**
 * @file ConvertMatrixStorageBenchmark.hpp
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
 * @brief Benchmark to measure the storage conversion times between different sparse matrix storages.
 * @author Jiri Kraus and Bea Hornef
 * @date 02.12.2011
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
using namespace bf;

template<typename StorageType1, typename StorageType2>
class ConvertMatrixStorageBenchmark : 

    public Benchmark, 
    private Benchmark::Register<ConvertMatrixStorageBenchmark< StorageType1, StorageType2> >

{
public:

    typedef StorageType1 SourceStorageType;
    typedef StorageType2 DestinationStorageType;
    typedef typename StorageType2::StorageValueType ValueType;

    static const std::string& id();

    ConvertMatrixStorageBenchmark();

    ConvertMatrixStorageBenchmark( const std::string& arguments );

    ConvertMatrixStorageBenchmark(
        const ConvertMatrixStorageBenchmark<SourceStorageType, DestinationStorageType>& other );

    virtual ~ConvertMatrixStorageBenchmark();

    /** Implementation of Benchmark::copy with covariant return type */

    virtual ConvertMatrixStorageBenchmark<SourceStorageType, DestinationStorageType>* copy() const;

    virtual short getValueTypeSize() const;

    virtual bool isThreadded() const;

    virtual const std::string& getId() const;

    static std::string createValue();

    static Benchmark* create()
    {
        return new ConvertMatrixStorageBenchmark();
    }

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

    SCAI_LOG_DECL_STATIC_LOGGER( logger );
};

#ifndef SCAI_LOG_LEVEL_OFF

#define LAMA_KOMMA ,
template<typename StorageType1, typename StorageType2>
SCAI_LOG_DEF_LOGGER( ConvertMatrixStorageBenchmark<StorageType1 LAMA_KOMMA StorageType2>::logger, "Benchmark.ConvertMatrixStorageBenchmark" );
#undef LAMA_KOMMA

#endif

template<typename StorageType1, typename StorageType2>
const std::string& ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename StorageType1, typename StorageType2>
const std::string& ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename StorageType1, typename StorageType2>
std::string ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::createValue()
{
    std::ostringstream value;
  
    value << "Convert_" << StorageType1::typeName() << "2" << StorageType2::typeName();

    return value.str();
}

template<typename StorageType1, typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::ConvertMatrixStorageBenchmark() :

    Benchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) )

{
    SCAI_LOG_DEBUG( logger, "ConvertMatrixStorageBenchmark created" );
    // means HOST
    mDevice = -1;
}

template<typename StorageType1, typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::ConvertMatrixStorageBenchmark(
    const ConvertMatrixStorageBenchmark<SourceStorageType, DestinationStorageType>& other )
    : Benchmark( other ), mDevice( other.mDevice )
{
    SCAI_LOG_DEBUG( logger, "copy constructror" )
}

template<typename StorageType1, typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::~ConvertMatrixStorageBenchmark()
{
    SCAI_LOG_DEBUG( logger, "~ConvertMatrixStorageBenchmark" )
}

template<typename StorageType1, typename StorageType2>
ConvertMatrixStorageBenchmark<StorageType1, StorageType2>* ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::copy() const
{
    SCAI_LOG_DEBUG( logger, "copy ConvertMatrixStorageBenchmark" )

    return  new ConvertMatrixStorageBenchmark<SourceStorageType, DestinationStorageType>( *this );
}

template<typename StorageType1, typename StorageType2>
short ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getValueTypeSize() const
{
    return sizeof( ValueType );
}

template<typename StorageType1, typename StorageType2>
bool ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::isThreadded() const
{
    return true;
}

template<typename StorageType1, typename StorageType2>
const std::string& ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getId() const
{
    static std::string id = createValue();
    return id;
}

template<typename StorageType1, typename StorageType2>
const LAMAInputSetComplexityVisitor::Group& ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::group()
{
    static const LAMAInputSetComplexityVisitor::Group group = LAMAInputSetComplexityVisitor::ConvertCSR2JDS;
    return group;
}

template<typename StorageType1, typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::initialize()
{
    SCAI_LOG_DEBUG( logger, "initialize, mInputSetId = " << mInputSetId )

    const LAMAInputSet& inputSet = InputSetRegistry<LAMAInputSet>::getRegistry().get( mInputSetId );

    SCAI_LOG_DEBUG( logger, "inputSet" )

    //TODO: assign should do the value type conversion

    const CSRStorage<double>& csrStorage = inputSet.getA().getLocalStorage();

    SCAI_LOG_DEBUG( logger, "csrStorage = " << csrStorage )

    mSourceStorage.setCSRData( csrStorage.getNumRows(), csrStorage.getNumColumns(), csrStorage.getNumValues(),
                               csrStorage.getIA(), csrStorage.getJA(), csrStorage.getValues() );

    ContextPtr context;

    if ( mDevice < 0 )
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

template<typename StorageType1, typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::setUp()
{
    mSourceStorage.prefetch( mSourceStorage.getContextPtr() );
    mDestinationStorage.prefetch( mDestinationStorage.getContextPtr() );
    mSourceStorage.wait();
    mDestinationStorage.wait();

//    mDestinationStorage._MatrixStorage::assign( mSourceStorage );
    mSourceStorage.buildCSRData( mCsrIA, mCsrJA, mCsrValues );

}

template<typename StorageType1, typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::execute()
{
    SCAI_LOG_DEBUG( logger, "execute conversion" )

    mDestinationStorage.setContextPtr( mSourceStorage.getContextPtr() );

    // Note: conversion will be done at context of target matrix mDestinationStorage

    mDestinationStorage.setCSRData( mSourceStorage.getNumRows(), mSourceStorage.getNumColumns(), mCsrJA.size(), mCsrIA,
                                    mCsrJA, mCsrValues );

    SCAI_LOG_DEBUG( logger, "execute conversion done" )
}

template<typename StorageType1, typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::tearDown()
{
}

template<typename StorageType1, typename StorageType2>
void ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::shutdown()
{
}

template<typename StorageType1, typename StorageType2>
CounterType ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getNumFloatingPointOperations() const
{
    return 0;
}

template<typename StorageType1, typename StorageType2>
CounterType ConvertMatrixStorageBenchmark<StorageType1, StorageType2>::getProcessedBytes() const
{
    return 0;
}

