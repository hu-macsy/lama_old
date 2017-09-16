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
#include <scai/benchmark/Parser.hpp>
#include <scai/lama/benchmark/LAMAInputSet.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/lama/storage/CSRStorage.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;
using namespace benchmark;

template<typename ValueType>
class ConvertMatrixStorageBenchmark : 

    public Benchmark, 
    private Benchmark::Register<ConvertMatrixStorageBenchmark<ValueType> >

{
public:

    ConvertMatrixStorageBenchmark( const std::string& arguments );

    virtual ~ConvertMatrixStorageBenchmark();

    virtual short getValueTypeSize() const;

    virtual bool isThreadded() const;

    virtual const std::string& getCreateId() const;

    virtual const std::string& getArguments() const;

    static std::string createValue();

    static Benchmark* create( std::string arguments )
    {
        return new ConvertMatrixStorageBenchmark( arguments );
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

    ConvertMatrixStorageBenchmark();

    common::unique_ptr<_MatrixStorage> mSourceStorage;
    common::unique_ptr<_MatrixStorage> mTargetStorage;

    std::string mArguments;

    SCAI_LOG_DECL_STATIC_LOGGER( logger );
};

template<typename ValueType>
SCAI_LOG_DEF_LOGGER( ConvertMatrixStorageBenchmark<ValueType>::logger, "Benchmark.ConvertMatrixStorageBenchmark" );

template<typename ValueType>
std::string ConvertMatrixStorageBenchmark<ValueType>::createValue()
{
    std::ostringstream value;
  
    value << "Convert_" << common::TypeTraits<ValueType>::id();

    return value.str();
}

template<typename ValueType>
ConvertMatrixStorageBenchmark<ValueType>::ConvertMatrixStorageBenchmark( const std::string& arg ) :

    Benchmark( getCreateId() + "( " + arg + " )" , "ConvertStorage" )

{
    SCAI_LOG_ERROR( logger, "ConvertMatrixStorageBenchmark( arg = " << arg << " )" );

    std::vector<std::string> storageTypes;

    if ( arg == "" )
    {
        tokenize( storageTypes, "CSR, JDS", ", " );
    }
    else
    {
        tokenize( storageTypes, arg, ", " );
    }

    SCAI_ASSERT_EQ_ERROR( storageTypes.size(), 2, "two storage types expected as arguments" )

    _MatrixStorage::MatrixStorageFormat sourceFormat = str2Format( storageTypes[0].c_str() );
    _MatrixStorage::MatrixStorageFormat targetFormat = str2Format( storageTypes[1].c_str() );

    common::scalar::ScalarType type = common::TypeTraits<ValueType>::stype;

    // allocate source and target storage of the required type

    mSourceStorage.reset( _MatrixStorage::create( MatrixStorageCreateKeyType( sourceFormat, type ) ) );
    mTargetStorage.reset( _MatrixStorage::create( MatrixStorageCreateKeyType( targetFormat, type ) ) );

    mArguments = storageTypes[0] + ", " + storageTypes[1];
}

template<typename ValueType>
ConvertMatrixStorageBenchmark<ValueType>::~ConvertMatrixStorageBenchmark()
{
    SCAI_LOG_DEBUG( logger, "~ConvertMatrixStorageBenchmark" )
}

template<typename ValueType>
short ConvertMatrixStorageBenchmark<ValueType>::getValueTypeSize() const
{
    return sizeof( ValueType );
}

template<typename ValueType>
bool ConvertMatrixStorageBenchmark<ValueType>::isThreadded() const
{
    return true;
}

template<typename ValueType>
const std::string& ConvertMatrixStorageBenchmark<ValueType>::getCreateId() const
{
    static std::string id = createValue();
    return id;
}

template<typename ValueType>
const std::string& ConvertMatrixStorageBenchmark<ValueType>::getArguments() const
{
    return mArguments;
}

template<typename ValueType>
void ConvertMatrixStorageBenchmark<ValueType>::initialize()
{
    SCAI_LOG_DEBUG( logger, "initialize, mInputSetId = " << mInputSetId )

    common::unique_ptr<InputSet> mInputSet;

    mInputSet.reset( benchmark::InputSet::parseAndCreate( mInputSetId ) );

    SCAI_LOG_ERROR( logger, "input set: " << *mInputSet )

    SCAI_ASSERT_EQ_ERROR( mInputSet->getId(), "LAMAInputSet", "Illegal LAMAInputSet: " << *mInputSet )

    // Now it is safe to cast

    const LAMAInputSet* mLAMAInputSet = reinterpret_cast<lama::LAMAInputSet*>( mInputSet.get() );

    //TODO: assign should do the value type conversion

    const CSRStorage<double>& csrStorage = mLAMAInputSet->getA().getLocalStorage();

    SCAI_LOG_DEBUG( logger, "csrStorage = " << csrStorage )

    mSourceStorage->setCSRData( csrStorage.getNumRows(), csrStorage.getNumColumns(), csrStorage.getNumValues(),
                                csrStorage.getIA(), csrStorage.getJA(), csrStorage.getValues() );

    ContextPtr context = Context::getContextPtr();

    mSourceStorage->setContextPtr( context );
    mTargetStorage->setContextPtr( context );
}

template<typename ValueType>
void ConvertMatrixStorageBenchmark<ValueType>::setUp()
{
    mSourceStorage->prefetch( mSourceStorage->getContextPtr() );
    mSourceStorage->wait();
 
    // no prefetch on target storage required
}

template<typename ValueType>
void ConvertMatrixStorageBenchmark<ValueType>::execute()
{
    *mTargetStorage = *mSourceStorage;
}

template<typename ValueType>
void ConvertMatrixStorageBenchmark<ValueType>::tearDown()
{
}

template<typename ValueType>
void ConvertMatrixStorageBenchmark<ValueType>::shutdown()
{
}

template<typename ValueType>
CounterType ConvertMatrixStorageBenchmark<ValueType>::getNumFloatingPointOperations() const
{
    return 0;
}

template<typename ValueType>
CounterType ConvertMatrixStorageBenchmark<ValueType>::getProcessedBytes() const
{
    return 0;
}

